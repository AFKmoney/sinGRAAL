// sinGRAAL — 6-automorphism Kangaroo CUDA kernel for secp256k1 ECDLP
//
// Each GPU thread is one independent Kangaroo animal.
// Jump table (128 entries) lives in constant memory (L1 broadcast).
// Every NORM_INTERVAL steps the animal is normalized to affine so the
// DP check uses the EXACT canonical x  (not a corrupt raw Jacobian limb).
// DPs written to a global ring buffer via atomicAdd.
//
// Compile:
//   nvcc -O3 -arch=sm_80 --compiler-options -fPIC -c kangaroo.cu -o kangaroo.o
//   ar rcs libkangaroo_cuda.a kangaroo.o

#include "secp256k1.cuh"
#include <stdint.h>
#include <string.h>
#include <stdio.h>

// ─── Tuning constants ─────────────────────────────────────────────────────────

#define NUM_JUMPS      128          // jump table size
#define MAX_DPS        (1u << 21)   // 2M DP ring buffer entries
#define BLOCK_SIZE     256
// Normalize every N steps: pays 1 fp_inv per N steps per animal
// N=512 → ~7% overhead; corrects all DP detection
#define NORM_INTERVAL  512

// ─── Jump table in constant memory ───────────────────────────────────────────

__constant__ JumpPoint g_jumps[NUM_JUMPS];

// ─── Device helpers ───────────────────────────────────────────────────────────

__device__ __forceinline__
u32 jump_idx(const u64 x[4]) { return (u32)(x[0] % (u64)NUM_JUMPS); }

// ─── Main Kangaroo kernel ─────────────────────────────────────────────────────
//
// Walk animals in Jacobian space for NORM_INTERVAL steps, then normalize to
// affine and check for a distinguished point using the exact canonical x.
// After the DP check the animal is reset to Z=1 (affine) so Z-coordinate
// drift never accumulates beyond NORM_INTERVAL multiplications.

__global__ void kangaroo_walk(
    Animal*  __restrict__ animals,
    DPEntry* __restrict__ dp_buf,
    u32*     __restrict__ dp_count,
    u32  num_animals,
    u32  steps_per_launch,   // rounded to NORM_INTERVAL multiple
    u64  dp_threshold,       // DP iff canon_x[3] < dp_threshold
    u32  max_dps
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_animals) return;

    Animal a = animals[tid];

    u32 steps_done = 0;
    while (steps_done < steps_per_launch) {
        // ── NORM_INTERVAL Jacobian steps ─────────────────────────────────────
        u32 batch = min((u32)NORM_INTERVAL, steps_per_launch - steps_done);
        for (u32 s = 0; s < batch; s++) {
            u32 ji = jump_idx(a.x);
            const JumpPoint jp = g_jumps[ji];

            u64 nx[4], ny[4], nz[4];
            pt_add_mixed(a.x, a.y, a.z, jp.x, jp.y, nx, ny, nz);
            for (int i = 0; i < 4; i++) { a.x[i]=nx[i]; a.y[i]=ny[i]; a.z[i]=nz[i]; }

            u64 ns[4];
            sc_add(a.scalar, jp.s, ns);
            for (int i = 0; i < 4; i++) a.scalar[i] = ns[i];
        }
        steps_done += batch;

        // ── Normalize to affine and exact canonical DP check ─────────────────
        u64 ax[4], ay[4], cx[4];
        pt_normalize(a.x, a.y, a.z, ax, ay);
        canonical_x_affine(ax, cx);

        if (cx[3] < dp_threshold) {
            u32 slot = atomicAdd(dp_count, 1u) % max_dps;
            DPEntry dp;
            for (int i = 0; i < 4; i++) {
                dp.canon_x[i] = cx[i];
                dp.scalar[i]  = a.scalar[i];
            }
            dp.is_wild  = a.is_wild;
            dp.pad[0] = dp.pad[1] = dp.pad[2] = 0;
            dp_buf[slot] = dp;
        }

        // Reset Z=1 (affine) — prevents Z-coordinate drift accumulation
        for (int i = 0; i < 4; i++) { a.x[i] = ax[i]; a.y[i] = ay[i]; }
        a.z[0] = 1; a.z[1] = a.z[2] = a.z[3] = 0;
    }

    animals[tid] = a;
}

// ─── Extern C API (called from Rust FFI) ─────────────────────────────────────

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
    cudaError_t e = cudaMemcpyToSymbol(g_jumps, jumps, n * sizeof(JumpPoint));
    return e == cudaSuccess ? 0 : -1;
}

KangarooCtx* kangaroo_init(
    const Animal* host_animals,
    u32   num_animals,
    u32   dp_bits,
    u32   steps_per_launch
) {
    KangarooCtx* ctx = new KangarooCtx();
    ctx->num_animals  = num_animals;
    ctx->max_dps      = MAX_DPS;
    ctx->dp_threshold = (dp_bits >= 64) ? 0ULL : (1ULL << (64 - dp_bits));

    // Round up to NORM_INTERVAL
    u32 r = ((steps_per_launch + NORM_INTERVAL - 1) / NORM_INTERVAL) * NORM_INTERVAL;
    ctx->steps_per_launch = r;
    ctx->grid = (num_animals + BLOCK_SIZE - 1) / BLOCK_SIZE;

    cudaMalloc(&ctx->d_animals,  num_animals * sizeof(Animal));
    cudaMalloc(&ctx->d_dp_buf,   MAX_DPS     * sizeof(DPEntry));
    cudaMalloc(&ctx->d_dp_count, sizeof(u32));
    cudaMemcpy(ctx->d_animals, host_animals, num_animals * sizeof(Animal),
               cudaMemcpyHostToDevice);
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
    to_read     = to_read < MAX_DPS ? to_read : MAX_DPS;
    if (to_read > 0)
        cudaMemcpy(host_buf, ctx->d_dp_buf, to_read * sizeof(DPEntry),
                   cudaMemcpyDeviceToHost);
    cudaMemset(ctx->d_dp_count, 0, sizeof(u32));
    return to_read;
}

void kangaroo_read_animals(KangarooCtx* ctx, Animal* host_buf) {
    cudaMemcpy(host_buf, ctx->d_animals, ctx->num_animals * sizeof(Animal),
               cudaMemcpyDeviceToHost);
}

void kangaroo_write_animals(KangarooCtx* ctx, const Animal* host_buf) {
    cudaMemcpy(ctx->d_animals, host_buf, ctx->num_animals * sizeof(Animal),
               cudaMemcpyHostToDevice);
}

void kangaroo_free(KangarooCtx* ctx) {
    if (!ctx) return;
    cudaFree(ctx->d_animals);
    cudaFree(ctx->d_dp_buf);
    cudaFree(ctx->d_dp_count);
    delete ctx;
}

u32 kangaroo_num_jumps()      { return NUM_JUMPS; }
u32 kangaroo_norm_interval()  { return NORM_INTERVAL; }

int  cuda_device_count() { int n=0; cudaGetDeviceCount(&n); return n; }
void cuda_set_device(int dev) { cudaSetDevice(dev); }

void cuda_device_name(int dev, char* buf, int len) {
    cudaDeviceProp p; cudaGetDeviceProperties(&p, dev);
    strncpy(buf, p.name, len-1); buf[len-1]='\0';
}

u64 cuda_device_memory(int dev) {
    cudaDeviceProp p; cudaGetDeviceProperties(&p, dev);
    return (u64)p.totalGlobalMem;
}

} // extern "C"
