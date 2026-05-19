// 6-automorphism Kangaroo CUDA kernel for secp256k1 ECDLP
//
// Each GPU thread is one independent Kangaroo animal.
// Jump table lives in constant memory (fast L1 broadcast).
// DPs written to a global ring buffer via atomicAdd.
//
// Compile:
//   nvcc -O3 -arch=sm_80 --compiler-options -fPIC -c kangaroo.cu -o kangaroo.o
//   ar rcs libkangaroo_cuda.a kangaroo.o

#include "secp256k1.cuh"
#include <stdint.h>
#include <string.h>
#include <stdio.h>

// ─── Jump table in constant memory ───────────────────────────────────────────

#define NUM_JUMPS 32
#define MAX_DPS   (1u << 20)   // 1M DP ring buffer
#define BLOCK_SIZE 256

__constant__ JumpPoint g_jumps[NUM_JUMPS];

// ─── Device helpers ───────────────────────────────────────────────────────────

__device__ __forceinline__
u32 jump_idx(const u64 x[4]) { return (u32)(x[0] % NUM_JUMPS); }

// Probabilistic DP: check raw Jacobian X top bits
// dp_mask = e.g. 0xFFFFFFF0 for 28-bit DPs
__device__ __forceinline__
bool is_dp(const u64 x[4], u64 dp_mask) { return (x[3] & dp_mask) == 0; }

// ─── Main Kangaroo kernel ─────────────────────────────────────────────────────

__global__ void kangaroo_walk(
    Animal*  __restrict__ animals,
    DPEntry* __restrict__ dp_buf,
    u32*     __restrict__ dp_count,
    u32  num_animals,
    u32  steps_per_launch,
    u64  dp_mask,
    u32  max_dps
) {
    u32 tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_animals) return;

    // Load animal into registers
    Animal a = animals[tid];

    for (u32 s = 0; s < steps_per_launch; s++) {
        u32 ji = jump_idx(a.x);
        const JumpPoint jp = g_jumps[ji];

        // Mixed Jacobian + Affine point addition
        u64 nx[4], ny[4], nz[4];
        pt_add_mixed(a.x, a.y, a.z, jp.x, jp.y, nx, ny, nz);
        for (int i = 0; i < 4; i++) { a.x[i] = nx[i]; a.y[i] = ny[i]; a.z[i] = nz[i]; }

        // Accumulate scalar mod n
        u64 ns[4];
        sc_add(a.scalar, jp.s, ns);
        for (int i = 0; i < 4; i++) a.scalar[i] = ns[i];

        // DP check (probabilistic on raw Jacobian X)
        if (is_dp(a.x, dp_mask)) {
            u32 slot = atomicAdd(dp_count, 1u) % max_dps;
            DPEntry dp;
            for (int i = 0; i < 4; i++) {
                dp.x_jac[i] = a.x[i];
                dp.z_jac[i] = a.z[i];
                dp.scalar[i] = a.scalar[i];
            }
            dp.is_wild = a.is_wild;
            dp.pad[0] = dp.pad[1] = dp.pad[2] = 0;
            dp_buf[slot] = dp;
        }
    }

    animals[tid] = a;
}

// ─── Extern C API (called from Rust) ─────────────────────────────────────────

extern "C" {

struct KangarooCtx {
    Animal*  d_animals;
    DPEntry* d_dp_buf;
    u32*     d_dp_count;
    u32      num_animals;
    u32      max_dps;
    u64      dp_mask;
    u32      grid;
};

// Upload jump table to constant memory
int kangaroo_set_jumps(const JumpPoint* jumps, int n) {
    if (n > NUM_JUMPS) return -1;
    cudaError_t e = cudaMemcpyToSymbol(g_jumps, jumps, n * sizeof(JumpPoint));
    return e == cudaSuccess ? 0 : -1;
}

// Allocate and initialize GPU buffers
KangarooCtx* kangaroo_init(
    const Animal* host_animals,
    u32   num_animals,
    u32   dp_bits
) {
    KangarooCtx* ctx = new KangarooCtx();
    ctx->num_animals = num_animals;
    ctx->max_dps     = MAX_DPS;
    // dp_mask: top dp_bits of x[3] must be zero
    // e.g. dp_bits=28: mask = 0xFFFFFFF0 (top 28 bits checked as zero in u64 top limb → 28/64 of top limb)
    // We check x[3] (bits 192-255 of x). Top 28 bits of x[3] = bits 228-255.
    // mask = 0xFFFFFFF0_00000000 → upper 28 bits of u64 = (0xFFFFFFF0ULL << 32)... hmm
    // Simpler: x[3] >> (64 - dp_bits) == 0  ↔  x[3] < (1 << (64-dp_bits))
    // We store it as a threshold: dp_threshold = 1ULL << (64 - dp_bits)
    u64 threshold = (dp_bits >= 64) ? 0 : (1ULL << (64 - dp_bits));
    ctx->dp_mask = threshold; // we'll check x[3] < dp_mask
    ctx->grid    = (num_animals + BLOCK_SIZE - 1) / BLOCK_SIZE;

    cudaMalloc(&ctx->d_animals,  num_animals * sizeof(Animal));
    cudaMalloc(&ctx->d_dp_buf,   MAX_DPS     * sizeof(DPEntry));
    cudaMalloc(&ctx->d_dp_count, sizeof(u32));

    cudaMemcpy(ctx->d_animals, host_animals, num_animals * sizeof(Animal), cudaMemcpyHostToDevice);
    cudaMemset(ctx->d_dp_count, 0, sizeof(u32));

    return ctx;
}

// Launch one kernel step, return number of DPs written
u32 kangaroo_step(KangarooCtx* ctx, u32 steps_per_launch) {
    // Reuse dp_mask as threshold (x[3] < dp_mask)
    kangaroo_walk<<<ctx->grid, BLOCK_SIZE>>>(
        ctx->d_animals,
        ctx->d_dp_buf,
        ctx->d_dp_count,
        ctx->num_animals,
        steps_per_launch,
        ctx->dp_mask,
        ctx->max_dps
    );
    cudaDeviceSynchronize();

    u32 count;
    cudaMemcpy(&count, ctx->d_dp_count, sizeof(u32), cudaMemcpyDeviceToHost);
    return count;
}

// Read DP entries from device, up to `max` entries
u32 kangaroo_read_dps(KangarooCtx* ctx, DPEntry* host_buf, u32 max) {
    u32 count;
    cudaMemcpy(&count, ctx->d_dp_count, sizeof(u32), cudaMemcpyDeviceToHost);
    u32 to_read = (count < max) ? count : max;
    to_read = (to_read < MAX_DPS) ? to_read : MAX_DPS;
    if (to_read > 0)
        cudaMemcpy(host_buf, ctx->d_dp_buf, to_read * sizeof(DPEntry), cudaMemcpyDeviceToHost);
    // Reset counter
    cudaMemset(ctx->d_dp_count, 0, sizeof(u32));
    return to_read;
}

// Get current animal positions (for monitoring)
void kangaroo_read_animals(KangarooCtx* ctx, Animal* host_buf) {
    cudaMemcpy(host_buf, ctx->d_animals, ctx->num_animals * sizeof(Animal), cudaMemcpyDeviceToHost);
}

void kangaroo_free(KangarooCtx* ctx) {
    if (!ctx) return;
    cudaFree(ctx->d_animals);
    cudaFree(ctx->d_dp_buf);
    cudaFree(ctx->d_dp_count);
    delete ctx;
}

// Query CUDA device info
int cuda_device_count() {
    int n = 0;
    cudaGetDeviceCount(&n);
    return n;
}

void cuda_device_name(int dev, char* buf, int len) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, dev);
    strncpy(buf, prop.name, len - 1);
    buf[len-1] = '\0';
}

u64 cuda_device_memory(int dev) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, dev);
    return (u64)prop.totalGlobalMem;
}

} // extern "C"
