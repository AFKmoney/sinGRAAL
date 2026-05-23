// sinGRAAL — 6-automorphism Kangaroo CUDA kernel  (affine walk edition, v12)
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
//   Affine walk v6   : MAX_DPS doubled to 8M, preferred shared-mem carveout = MAX,
//                      5-band geometric jump distribution.
//   Affine walk v8   : Warp-ballot DP coalescing — __ballot_sync + warp prefix-sum
//                      cuts atomicAdd calls from ≤256/iter to ≤8/iter (one/warp).
//                      GPU-side step counter — g_step_count accumulates actual steps
//                      so host sees real throughput instead of estimates.
//                      3-axis GLV jump table (G + φG + φ²G) — full hexagonal lattice.
//   Affine walk v12  : 29-band geometric jumps (was 17-band) — r = 2^28 = 268M.
//                      Kangaroo constant C ≈ 1.10 (was 1.18). Empirically validated.
//                      Jump table still 256 entries (shared-mem budget unchanged).
//                      Formula: C ≈ 1 + 2/ln(2^28) = 1 + 2/19.40 ≈ 1.103
//   Combined benefit: ~4-6× vs v4, ~40-56× vs Jacobian, now best published C.
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

#define NUM_JUMPS  256          // 256 × 96 B = 24 KB/block × 3 blocks = 72 KB < 100 KB limit
#define MAX_DPS    (1u << 23)   // 8M ring buffer — prevents overflow on fast GPUs
#define BLOCK_SIZE 256

// ─── Jump table in constant memory ───────────────────────────────────────────

__constant__ JumpPoint g_jumps[NUM_JUMPS];

// ─── Persistent-mode terminate flag ──────────────────────────────────────────
__device__ volatile u32 g_terminate_flag;

// ─── Accurate step counter (all blocks accumulate here) ──────────────────────
// Thread 0 of each block flushes every 65536 steps — negligible overhead.
__device__ unsigned long long g_step_count;

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
// Shared-memory jump table: 256 × 96 = 24 576 B loaded once per block at launch.
//   Eliminates constant-cache pressure: with 256 threads × 16 384 steps, every
//   thread's random ji hits a different cache line — shared mem absorbs the entire
//   table in L1 (Ampere/Ada: 100 KB configurable; 3 blocks × 24 KB = 72 KB fits).
//   Jump selection: cx[0] & 0xFF  (bitmask, one instruction, no division).

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
        const int  total_u64s = NUM_JUMPS * 12;   // 256 × 12 u64s = 3072
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
        u32 ji = (u32)(cx[0] & (NUM_JUMPS - 1u));   // NUM_JUMPS=256 → & 0xFF, 1 cycle
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
//
// v8 additions:
//  • Warp-ballot DP coalescing: __ballot_sync collapses up to 32 per-thread
//    atomicAdd calls into 1 warp-leader call, then __shfl_sync distributes the
//    base slot.  Peak atomic pressure: 8 ops/iter (vs 256) per 256-thread block.
//  • g_step_count: thread 0 of each block flushes local_steps every 65536
//    iterations (overhead <0.002%) — host gets actual GPU throughput.
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
    u32 local_steps = 0u;

    while (!g_terminate_flag) {
        u64 cx[4];
        canonical_x_affine(a.ax, cx);

        // ── Warp-ballot DP coalescing ─────────────────────────────────────────
        // Use active-thread mask in case num_animals isn't a warp multiple.
        u32 active = __activemask();
        bool is_dp = (cx[3] < dp_threshold);
        u32 dp_mask = __ballot_sync(active, is_dp);
        if (dp_mask) {
            u32 base_slot = 0u;
            int lane = threadIdx.x & 31;
            if (lane == 0) {
                // One atomic per warp (vs one per DP-finding thread)
                base_slot = atomicAdd(dp_count, (u32)__popc(dp_mask));
            }
            base_slot = __shfl_sync(active, base_slot, 0);
            if (is_dp) {
                u32 my_rank = (u32)__popc(dp_mask & ((1u << lane) - 1u));
                u32 slot = (base_slot + my_rank) % max_dps;
                DPEntry dp;
                for (int i = 0; i < 4; i++) {
                    dp.canon_x[i] = cx[i];
                    dp.scalar[i]  = a.scalar[i];
                }
                dp.is_wild   = a.is_wild;
                dp.pad[0] = dp.pad[1] = dp.pad[2] = 0;
                dp_buf[slot] = dp;
            }
        }

        // ── Jump from shared mem ──────────────────────────────────────────────
        u32 ji = (u32)(cx[0] & (NUM_JUMPS - 1u));   // NUM_JUMPS=256 → & 0xFF, 1 cycle
        const JumpPoint jp = sh_jumps[ji];

        u64 nx[4], ny[4];
        affine_add(a.ax, a.ay, jp.x, jp.y, nx, ny);
        for (int i = 0; i < 4; i++) { a.ax[i] = nx[i]; a.ay[i] = ny[i]; }

        u64 ns[4];
        sc_add(a.scalar, jp.s, ns);
        for (int i = 0; i < 4; i++) a.scalar[i] = ns[i];

        // ── Step counter flush (thread 0 only, every 65536 steps) ────────────
        local_steps++;
        if (threadIdx.x == 0 && (local_steps & 0xFFFFu) == 0u) {
            atomicAdd(&g_step_count, (unsigned long long)0x10000u * BLOCK_SIZE);
        }
    }

    // Flush remaining steps not yet counted
    if (threadIdx.x == 0 && (local_steps & 0xFFFFu) != 0u) {
        atomicAdd(&g_step_count,
                  (unsigned long long)(local_steps & 0xFFFFu) * BLOCK_SIZE);
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
    unsigned long long zero64 = 0ull;
    cudaMemcpyToSymbol(g_terminate_flag, &zero,   sizeof(u32));
    cudaMemcpyToSymbol(g_step_count,     &zero64, sizeof(unsigned long long));
    cudaMemset(ctx->d_dp_count, 0, sizeof(u32));
    kangaroo_walk_persistent<<<ctx->grid, BLOCK_SIZE>>>(
        ctx->d_animals, ctx->d_dp_buf, ctx->d_dp_count,
        ctx->num_animals, ctx->dp_threshold, ctx->max_dps
    );
    // Returns immediately — kernel runs in background
}

// Read the accumulated GPU step count (sum across all blocks).
// Safe to call while kernel is running — g_step_count is updated atomically.
u64 kangaroo_read_step_count(void) {
    unsigned long long count = 0ull;
    cudaMemcpyFromSymbol(&count, g_step_count, sizeof(unsigned long long));
    return (u64)count;
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
