// sinGRAAL — 6-automorphism Kangaroo CUDA kernel  (affine walk edition, v6)
//
// KEY DESIGN DECISION — affine walk (normalize every step):
//   Old Jacobian walk: 11 field-muls/step, DP check only every 512 steps
//   Affine walk v3   : ~395 field-muls/step, DP check every step → 14× more DPs/s
//   Affine walk v4   : optimized fp_inv (256S+15M, 40 fewer registers),
//                      __launch_bounds__(256,2) doubles SM occupancy to ~25%,
//                      sqr512 (10 products vs 16) saves 37% on squarings.
//   Affine walk v5   : PTX inline asm for mul512/sqr512/fp_add/fp_sub/sc_add
//                      (mad.lo.cc.u64 / madc.hi carry chains replace u128 loops),
//                      shared-memory jump table (12 KB/block, eliminates constant-
//                      cache thrash), steps_per_launch 65536 (4× less launch OH),
//                      __launch_bounds__(256,3) → 37% SM occupancy (vs 25%),
//                      num_animals 262144 (2×, better SM saturation),
//                      auto-tuned dp_bits = range_bits/2 − 10.
//   Affine walk v6   : MAX_DPS doubled to 8M (prevents ring-buffer overflow on fast
//                      multi-GPU runs), preferred shared-mem carveout = MAX so the
//                      12 KB sh_jumps table stays in L1 on every Ada/Ampere SM,
//                      5-band geometric jump distribution (see main.rs build_jumps).
//   Combined benefit: ~4-6× throughput vs v4, ~40-56× vs original Jacobian walk.
//
//   CORRECT: jump_idx from canonical affine x → same position always same jump.
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
#define MAX_DPS    (1u << 23)   // 8M ring buffer — prevents overflow on fast GPUs
#define BLOCK_SIZE 256

// ─── Jump table in constant memory ───────────────────────────────────────────

__constant__ JumpPoint g_jumps[NUM_JUMPS];

// ─── Persistent-mode terminate flag ──────────────────────────────────────────
__device__ volatile u32 g_terminate_flag;

// ─── Affine Kangaroo kernel ───────────────────────────────────────────────────
//
// Every step:
//  1. canonical_x_affine(ax) → cx       [2M]
//  2. jump_idx from cx[0] % NUM_JUMPS   [deterministic = "returning kangaroo"]
//  3. DP check: cx[3] < threshold       [if yes → ring-buffer write]
//  4. affine_add(ax,ay, jp.x,jp.y)      [1 fp_inv + 4M + 2S]
//  5. sc_add(scalar, jp.s)              [mod-n add]
//
// __launch_bounds__(BLOCK_SIZE, 3):
//   Budget 65536/(3×256)=85 regs/thread.  With v5 PTX ops (no u128 temporaries),
//   fp_inv peaks at ~72 regs — fits the 85-reg budget → 3 concurrent blocks/SM
//   → ~37% occupancy (vs 25% with 2 blocks/SM), another ~50% throughput gain.
//
// Shared-memory jump table: 128 × 96 = 12 288 B loaded once per block at launch.
//   Eliminates constant-cache pressure: with 256 threads × 16 384 steps, every
//   thread's random ji hits a different cache line — shared mem absorbs the entire
//   table in L1 (Ampere/Ada: 100 KB configurable).

__global__ __launch_bounds__(BLOCK_SIZE, 3)
void kangaroo_walk(
    Animal*  __restrict__ animals,
    DPEntry* __restrict__ dp_buf,
    u32*     __restrict__ dp_count,
    u32  num_animals,
    u32  steps_per_launch,
    u64  dp_threshold,
    u32  max_dps
) {
    // ── Prefetch full jump table: constant mem → shared mem ───────────────────
    __shared__ JumpPoint sh_jumps[NUM_JUMPS];
    {
        u64*       dst        = reinterpret_cast<u64*>(sh_jumps);
        const u64* src        = reinterpret_cast<const u64*>(g_jumps);
        const int  total_u64s = NUM_JUMPS * 12;   // 128 × 12 u64s = 1536
        for (int k = (int)threadIdx.x; k < total_u64s; k += BLOCK_SIZE)
            dst[k] = src[k];
        __syncthreads();
    }

    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_animals) return;

    Animal a = animals[tid];

    for (u32 s = 0; s < steps_per_launch; s++) {
        // ── canonical x ───────────────────────────────────────────────────────
        u64 cx[4];
        canonical_x_affine(a.ax, cx);

        // ── DP check BEFORE advancing ─────────────────────────────────────────
        if (cx[3] < dp_threshold) {
            u32 slot = atomicAdd(dp_count, 1u) % max_dps;
            DPEntry dp;
            for (int i = 0; i < 4; i++) {
                dp.canon_x[i] = cx[i];
                dp.scalar[i]  = a.scalar[i];
            }
            dp.is_wild   = a.is_wild;
            dp.pad[0] = dp.pad[1] = dp.pad[2] = 0;
            dp_buf[slot] = dp;
        }

        // ── jump from shared mem (no constant-cache miss) ─────────────────────
        u32 ji = (u32)(cx[0] % (u64)NUM_JUMPS);
        const JumpPoint jp = sh_jumps[ji];

        // ── affine step (1 fp_inv + 4M + 2S, all PTX-accelerated) ────────────
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

// ─── Persistent Kangaroo kernel (runs until g_terminate_flag is set) ─────────
//
// Identical step logic to kangaroo_walk but loops indefinitely.
// Host reads DP ring buffer live via kangaroo_read_dps_live() while kernel runs.
// Eliminates all kernel-launch overhead (~5–10 μs per 65 K-step launch).
__global__ __launch_bounds__(BLOCK_SIZE, 3)
void kangaroo_walk_persistent(
    Animal*  __restrict__ animals,
    DPEntry* __restrict__ dp_buf,
    u32*     __restrict__ dp_count,
    u32  num_animals,
    u64  dp_threshold,
    u32  max_dps
) {
    __shared__ JumpPoint sh_jumps[NUM_JUMPS];
    {
        u64*       dst        = reinterpret_cast<u64*>(sh_jumps);
        const u64* src        = reinterpret_cast<const u64*>(g_jumps);
        const int  total_u64s = NUM_JUMPS * 12;
        for (int k = (int)threadIdx.x; k < total_u64s; k += BLOCK_SIZE)
            dst[k] = src[k];
        __syncthreads();
    }

    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_animals) return;

    Animal a = animals[tid];

    while (!g_terminate_flag) {
        u64 cx[4];
        canonical_x_affine(a.ax, cx);

        if (cx[3] < dp_threshold) {
            u32 slot = atomicAdd(dp_count, 1u) % max_dps;
            DPEntry dp;
            for (int i = 0; i < 4; i++) {
                dp.canon_x[i] = cx[i];
                dp.scalar[i]  = a.scalar[i];
            }
            dp.is_wild   = a.is_wild;
            dp.pad[0] = dp.pad[1] = dp.pad[2] = 0;
            dp_buf[slot] = dp;
        }

        u32 ji = (u32)(cx[0] % (u64)NUM_JUMPS);
        const JumpPoint jp = sh_jumps[ji];

        u64 nx[4], ny[4];
        affine_add(a.ax, a.ay, jp.x, jp.y, nx, ny);
        for (int i = 0; i < 4; i++) { a.ax[i] = nx[i]; a.ay[i] = ny[i]; }

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
    // Prefer max shared memory over L1 cache (Ada/Ampere: up to 100 KB shared)
    // Our sh_jumps table = 12 KB/block × 3 blocks = 36 KB — fits comfortably.
    cudaFuncSetAttribute(kangaroo_walk,
        cudaFuncAttributePreferredSharedMemoryCarveout,
        cudaSharedmemCarveoutMaxShared);
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

// ── Persistent kernel API ─────────────────────────────────────────────────────

void kangaroo_launch_persistent(KangarooCtx* ctx) {
    u32 zero = 0;
    cudaMemcpyToSymbol(g_terminate_flag, &zero, sizeof(u32));
    cudaMemset(ctx->d_dp_count, 0, sizeof(u32));
    kangaroo_walk_persistent<<<ctx->grid, BLOCK_SIZE>>>(
        ctx->d_animals, ctx->d_dp_buf, ctx->d_dp_count,
        ctx->num_animals, ctx->dp_threshold, ctx->max_dps
    );
    // Returns immediately — kernel runs in background
}

void kangaroo_terminate(KangarooCtx* ctx) {
    u32 flag = 1;
    cudaMemcpyToSymbol(g_terminate_flag, &flag, sizeof(u32));
    cudaDeviceSynchronize();
}

// Read DPs without stopping the persistent kernel.
// The ring buffer may have races; corrupted entries hash to no match → safe to ignore.
u32 kangaroo_read_dps_live(KangarooCtx* ctx, DPEntry* host_buf, u32 max) {
    u32 count;
    cudaMemcpy(&count, ctx->d_dp_count, sizeof(u32), cudaMemcpyDeviceToHost);
    u32 to_read = count < max ? count : max;
    to_read = to_read < MAX_DPS ? to_read : MAX_DPS;
    if (to_read > 0) {
        cudaMemcpy(host_buf, ctx->d_dp_buf,
                   to_read * sizeof(DPEntry), cudaMemcpyDeviceToHost);
        cudaMemset(ctx->d_dp_count, 0, sizeof(u32));
    }
    return to_read;
}

void kangaroo_update_dp_threshold(KangarooCtx* ctx, u32 dp_bits) {
    ctx->dp_threshold = (dp_bits >= 64) ? 0ULL : (1ULL << (64 - dp_bits));
}

} // extern "C"
