// sinGRAAL — Kangaroo CUDA kernel with Toom-Cook-6 field multiplication (v13)
//
// Drop-in replacement for kangaroo.cu with Toom-6 fp_mul:
//   Schoolbook (4×4): 16 MAD instructions per field multiplication
//   Toom-Cook-6:      11 MAD instructions + forward-difference interpolation
//   Net:              31% fewer multiplications in the hottest path
//
// Architecture: persistent kernel + warp-ballot DP coalescing (v8+)
//   Every affine_add calls fp_inv (256S + 15M) which calls fp_mul repeatedly.
//   Each fp_mul replaces 1 schoolbook mul512 with 1 toom6 mul512.
//   fp_inv: 256 sqr512 + 15 mul512_toom6 = 256*SQR + 15*11 MAD (vs 15*16 MAD)
//
// Compile:
//   nvcc -O3 -arch=sm_80 --compiler-options -fPIC \
//        -c kangaroo_toom6.cu -o kangaroo_toom6.o
//   nvcc -O3 -arch=sm_80 --compiler-options -fPIC \
//        -DUSE_SCHOOLBOOK -c kangaroo_toom6.cu -o kangaroo_schoolbook.o

#include "secp256k1_toom6.cuh"
#include <stdint.h>
#include <string.h>
#include <stdio.h>

// ─── Tuning ───────────────────────────────────────────────────────────────────

#define NUM_JUMPS  256
#define MAX_DPS    (1u << 23)
#define BLOCK_SIZE 256

// ─── Constant memory ─────────────────────────────────────────────────────────

__constant__ JumpPoint g_jumps[NUM_JUMPS];
__device__ volatile u32 g_terminate_flag;
__device__ unsigned long long g_step_count;

// ─── Persistent Kangaroo kernel ───────────────────────────────────────────────
//
// Step logic (identical to v13 kangaroo.cu, fp_mul now uses Toom-6):
//  1. canonical_x_affine(ax, cx)     [2 fp_mul = 2 × Toom-6]
//  2. DP check: cx[3] < dp_threshold
//  3. Warp-ballot DP coalescing (v8 optimization)
//  4. affine_add(ax,ay, jp.x,jp.y)  [1 fp_inv + 4 fp_mul + 2 fp_sqr]
//  5. sc_add(scalar, jp.s)
//
// fp_inv = 256 sqr512 + 15 mul512.  With Toom-6: 15 × 11 = 165 MADs
// vs schoolbook: 15 × 16 = 240 MADs.  Saving: 75 MADs per inverse.
// Total per step: ~165 + 256*10 (sqr PTX) + 2*11 (canonical) = significant.

__global__ __launch_bounds__(BLOCK_SIZE, 3)
void kangaroo_walk_persistent_toom6(
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

        // Warp-ballot DP coalescing
        u32 active  = __activemask();
        bool is_dp  = (cx[3] < dp_threshold);
        u32 dp_mask = __ballot_sync(active, is_dp);
        if (dp_mask) {
            u32 base_slot = 0u;
            int lane = threadIdx.x & 31;
            if (lane == 0)
                base_slot = atomicAdd(dp_count, (u32)__popc(dp_mask));
            base_slot = __shfl_sync(active, base_slot, 0);
            if (is_dp) {
                u32 my_rank = (u32)__popc(dp_mask & ((1u << lane) - 1u));
                u32 slot = (base_slot + my_rank) % max_dps;
                DPEntry dp;
                for (int i = 0; i < 4; i++) {
                    dp.canon_x[i] = cx[i];
                    dp.scalar[i]  = a.scalar[i];
                }
                dp.is_wild = a.is_wild;
                dp.pad[0] = dp.pad[1] = dp.pad[2] = 0;
                dp_buf[slot] = dp;
            }
        }

        u32 ji = (u32)(cx[0] & (NUM_JUMPS - 1u));
        const JumpPoint jp = sh_jumps[ji];

        u64 nx[4], ny[4];
        affine_add(a.ax, a.ay, jp.x, jp.y, nx, ny);
        for (int i = 0; i < 4; i++) { a.ax[i] = nx[i]; a.ay[i] = ny[i]; }

        u64 ns[4];
        sc_add(a.scalar, jp.s, ns);
        for (int i = 0; i < 4; i++) a.scalar[i] = ns[i];

        local_steps++;
        if (threadIdx.x == 0 && (local_steps & 0xFFFFu) == 0u)
            atomicAdd(&g_step_count, (unsigned long long)0x10000u * BLOCK_SIZE);
    }

    if (threadIdx.x == 0 && (local_steps & 0xFFFFu) != 0u)
        atomicAdd(&g_step_count,
                  (unsigned long long)(local_steps & 0xFFFFu) * BLOCK_SIZE);

    animals[tid] = a;
}

// ─── Batched (non-persistent) kernel ─────────────────────────────────────────

__global__ __launch_bounds__(BLOCK_SIZE, 3)
void kangaroo_walk_toom6(
    Animal*  __restrict__ animals,
    DPEntry* __restrict__ dp_buf,
    u32*     __restrict__ dp_count,
    u32  num_animals,
    u32  steps_per_launch,
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

    for (u32 s = 0; s < steps_per_launch; s++) {
        u64 cx[4];
        canonical_x_affine(a.ax, cx);

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

        u32 ji = (u32)(cx[0] & (NUM_JUMPS - 1u));
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
    ctx->num_animals      = num_animals;
    ctx->max_dps          = MAX_DPS;
    ctx->dp_threshold     = (dp_bits >= 64) ? 0ULL : (1ULL << (64 - dp_bits));
    ctx->steps_per_launch = steps_per_launch;
    ctx->grid             = (num_animals + BLOCK_SIZE - 1) / BLOCK_SIZE;

    cudaMalloc(&ctx->d_animals,  num_animals * sizeof(Animal));
    cudaMalloc(&ctx->d_dp_buf,   MAX_DPS     * sizeof(DPEntry));
    cudaMalloc(&ctx->d_dp_count, sizeof(u32));
    cudaMemcpy(ctx->d_animals, host_animals,
               num_animals * sizeof(Animal), cudaMemcpyHostToDevice);
    cudaMemset(ctx->d_dp_count, 0, sizeof(u32));
    cudaFuncSetAttribute(kangaroo_walk_persistent_toom6,
        cudaFuncAttributePreferredSharedMemoryCarveout,
        cudaSharedmemCarveoutMaxShared);
    return ctx;
}

u32 kangaroo_step(KangarooCtx* ctx) {
    kangaroo_walk_toom6<<<ctx->grid, BLOCK_SIZE>>>(
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

u32  kangaroo_num_jumps()         { return NUM_JUMPS; }
int  cuda_device_count()          { int n=0; cudaGetDeviceCount(&n); return n; }
void cuda_set_device(int dev)     { cudaSetDevice(dev); }

void cuda_device_name(int dev, char* buf, int len) {
    cudaDeviceProp p; cudaGetDeviceProperties(&p, dev);
    strncpy(buf, p.name, len-1); buf[len-1] = '\0';
}

u64 cuda_device_memory(int dev) {
    cudaDeviceProp p; cudaGetDeviceProperties(&p, dev);
    return (u64)p.totalGlobalMem;
}

void kangaroo_launch_persistent(KangarooCtx* ctx) {
    u32 zero = 0; unsigned long long zero64 = 0ull;
    cudaMemcpyToSymbol(g_terminate_flag, &zero,   sizeof(u32));
    cudaMemcpyToSymbol(g_step_count,     &zero64, sizeof(unsigned long long));
    cudaMemset(ctx->d_dp_count, 0, sizeof(u32));
    kangaroo_walk_persistent_toom6<<<ctx->grid, BLOCK_SIZE>>>(
        ctx->d_animals, ctx->d_dp_buf, ctx->d_dp_count,
        ctx->num_animals, ctx->dp_threshold, ctx->max_dps
    );
}

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
