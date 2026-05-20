// sinGRAAL — 6-automorphism Kangaroo CUDA kernel  (affine walk edition)
//
// KEY DESIGN DECISION — affine walk (normalize every step):
//   Old Jacobian walk: 11 field-muls/step, but DP check only every 512 steps
//   New affine walk  : ~395 field-muls/step, DP check every single step
//
//   Net gain: 512 / (395/11) ≈ 14× more effective DPs per wall-second.
//   ALSO CORRECT: jump_idx uses canonical affine x so two animals at the
//   same position always take the same jump (required for the Kangaroo
//   birthday argument to hold).
//
// Compile:
//   nvcc -O3 -arch=sm_80 --compiler-options -fPIC -c kangaroo.cu -o kangaroo.o
//   ar rcs libkangaroo_cuda.a kangaroo.o

#include "secp256k1.cuh"
#include <stdint.h>
#include <string.h>
#include <stdio.h>

// ─── Tuning constants ─────────────────────────────────────────────────────────

#define NUM_JUMPS  128
#define MAX_DPS    (1u << 22)   // 4M ring buffer
#define BLOCK_SIZE 256

// ─── Jump table in constant memory ───────────────────────────────────────────

__constant__ JumpPoint g_jumps[NUM_JUMPS];

// ─── Affine Kangaroo kernel ───────────────────────────────────────────────────
//
// Every step:
//  1. canonical_x_affine(ax) → cx       [2M]
//  2. jump_idx from cx[0] % NUM_JUMPS   [deterministic = "returning kangaroo"]
//  3. DP check: cx[3] < threshold       [if yes → ring-buffer write]
//  4. affine_add(ax,ay, jp.x,jp.y)      [1 fp_inv + 4M + 2S]
//  5. sc_add(scalar, jp.s)              [mod-n add]

__global__ void kangaroo_walk(
    Animal*  __restrict__ animals,
    DPEntry* __restrict__ dp_buf,
    u32*     __restrict__ dp_count,
    u32  num_animals,
    u32  steps_per_launch,
    u64  dp_threshold,
    u32  max_dps
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_animals) return;

    Animal a = animals[tid];

    for (u32 s = 0; s < steps_per_launch; s++) {
        // ── canonical x (affine, no inversion needed — ax is already affine) ──
        u64 cx[4];
        canonical_x_affine(a.ax, cx);

        // ── DP check BEFORE advancing (record current exact position) ─────────
        if (cx[3] < dp_threshold) {
            u32 slot = atomicAdd(dp_count, 1u) % max_dps;
            DPEntry dp;
            for (int i = 0; i < 4; i++) {
                dp.canon_x[i] = cx[i];
                dp.scalar[i]  = a.scalar[i];
            }
            dp.is_wild = a.is_wild;
            dp.pad[0] = dp.pad[1] = dp.pad[2] = 0;
            dp_buf[slot] = dp;
        }

        // ── deterministic jump selection from canonical x ─────────────────────
        u32 ji = (u32)(cx[0] % (u64)NUM_JUMPS);
        const JumpPoint jp = g_jumps[ji];

        // ── affine step (1 fp_inv + 4M + 2S) ─────────────────────────────────
        u64 nx[4], ny[4];
        affine_add(a.ax, a.ay, jp.x, jp.y, nx, ny);
        for (int i = 0; i < 4; i++) { a.ax[i] = nx[i]; a.ay[i] = ny[i]; }

        // ── scalar accumulation ───────────────────────────────────────────────
        u64 ns[4];
        sc_add(a.scalar, jp.s, ns);
        for (int i = 0; i < 4; i++) a.scalar[i] = ns[i];
    }

    animals[tid] = a;
}

// ─── Extern C API ─────────────────────────────────────────────────────────────

extern "C" {

struct KangarooCtx {
    Animal*  d_animals;
    DPEntry* d_dp_buf;
    u32*     d_dp_count;
    u32      num_animals;
    u32      max_dps;
    u64      dp_threshold;
    u32      grid;
    u32      steps_per_launch;
};

int kangaroo_set_jumps(const JumpPoint* jumps, int n) {
    if (n > NUM_JUMPS) return -1;
    return cudaMemcpyToSymbol(g_jumps, jumps, n * sizeof(JumpPoint))
           == cudaSuccess ? 0 : -1;
}

KangarooCtx* kangaroo_init(
    const Animal* host_animals,
    u32   num_animals,
    u32   dp_bits,
    u32   steps_per_launch
) {
    KangarooCtx* ctx = new KangarooCtx();
    ctx->num_animals     = num_animals;
    ctx->max_dps         = MAX_DPS;
    ctx->dp_threshold    = (dp_bits >= 64) ? 0ULL : (1ULL << (64 - dp_bits));
    ctx->steps_per_launch = steps_per_launch;
    ctx->grid            = (num_animals + BLOCK_SIZE - 1) / BLOCK_SIZE;

    cudaMalloc(&ctx->d_animals,  num_animals * sizeof(Animal));
    cudaMalloc(&ctx->d_dp_buf,   MAX_DPS     * sizeof(DPEntry));
    cudaMalloc(&ctx->d_dp_count, sizeof(u32));
    cudaMemcpy(ctx->d_animals, host_animals,
               num_animals * sizeof(Animal), cudaMemcpyHostToDevice);
    cudaMemset(ctx->d_dp_count, 0, sizeof(u32));
    return ctx;
}

u32 kangaroo_step(KangarooCtx* ctx) {
    kangaroo_walk<<<ctx->grid, BLOCK_SIZE>>>(
        ctx->d_animals, ctx->d_dp_buf, ctx->d_dp_count,
        ctx->num_animals, ctx->steps_per_launch,
        ctx->dp_threshold, ctx->max_dps
    );
    cudaDeviceSynchronize();
    u32 count;
    cudaMemcpy(&count, ctx->d_dp_count, sizeof(u32), cudaMemcpyDeviceToHost);
    return count;
}

u32 kangaroo_read_dps(KangarooCtx* ctx, DPEntry* host_buf, u32 max) {
    u32 count;
    cudaMemcpy(&count, ctx->d_dp_count, sizeof(u32), cudaMemcpyDeviceToHost);
    u32 to_read = count < max ? count : max;
    to_read = to_read < MAX_DPS ? to_read : MAX_DPS;
    if (to_read > 0)
        cudaMemcpy(host_buf, ctx->d_dp_buf,
                   to_read * sizeof(DPEntry), cudaMemcpyDeviceToHost);
    cudaMemset(ctx->d_dp_count, 0, sizeof(u32));
    return to_read;
}

void kangaroo_read_animals(KangarooCtx* ctx, Animal* host_buf) {
    cudaMemcpy(host_buf, ctx->d_animals,
               ctx->num_animals * sizeof(Animal), cudaMemcpyDeviceToHost);
}

void kangaroo_write_animals(KangarooCtx* ctx, const Animal* host_buf) {
    cudaMemcpy(ctx->d_animals, host_buf,
               ctx->num_animals * sizeof(Animal), cudaMemcpyHostToDevice);
}

void kangaroo_free(KangarooCtx* ctx) {
    if (!ctx) return;
    cudaFree(ctx->d_animals);
    cudaFree(ctx->d_dp_buf);
    cudaFree(ctx->d_dp_count);
    delete ctx;
}

u32 kangaroo_num_jumps() { return NUM_JUMPS; }

int  cuda_device_count()       { int n=0; cudaGetDeviceCount(&n); return n; }
void cuda_set_device(int dev)  { cudaSetDevice(dev); }

void cuda_device_name(int dev, char* buf, int len) {
    cudaDeviceProp p; cudaGetDeviceProperties(&p, dev);
    strncpy(buf, p.name, len-1); buf[len-1] = '\0';
}

u64 cuda_device_memory(int dev) {
    cudaDeviceProp p; cudaGetDeviceProperties(&p, dev);
    return (u64)p.totalGlobalMem;
}

} // extern "C"
