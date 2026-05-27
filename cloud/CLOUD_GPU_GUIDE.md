# sinGRAAL v13 — Cloud GPU Deployment Guide (Puzzle #135)

Deploy sinGRAAL across cloud GPU providers in minutes using Docker.

---

## GPU Architecture — Automatic Detection (zero config)

**You no longer need to pick `CUDA_ARCH`.** With `CUDA_ARCH` unset (the default),
`build.rs` compiles a **portable multi-architecture fatbinary** covering every
common datacenter and consumer GPU, plus embedded PTX for forward-compat JIT.
The CUDA driver selects the matching cubin at load time, so a single image runs
unchanged on any of these:

| GPU | compute_cap | Covered by fatbinary |
|-----|-------------|----------------------|
| H100 | 9.0 (`sm_90`) | ✅ |
| A100 | 8.0 (`sm_80`) | ✅ |
| RTX 4090 / L4 / L40 | 8.9 (`sm_89`) | ✅ |
| RTX 3090 / 3080 Ti | 8.6 (`sm_86`) | ✅ |
| RTX 2080 Ti / T4 | 7.5 (`sm_75`) | ✅ |
| V100 | 7.0 (`sm_70`) | ✅ |
| newer (Blackwell, …) | ≥ 9.0 | ✅ via PTX JIT |

At runtime the solver also auto-detects every visible GPU (name + VRAM via
`nvidia-smi`) and, if `NUM_ANIMALS` is unset, auto-calibrates the animal count
to the card's memory. **No flags required — just provide the target.**

Optional: to shrink the image / speed up the build, pin one arch:
```bash
docker build --build-arg CUDA_ARCH=sm_89 .   # single-arch (e.g. RTX 4090)
```
If you build *on the target GPU machine* with `CUDA_ARCH` unset, `build.rs`
detects the local card's compute capability and builds natively for it.

---

## RunPod

### Quick Start (Single GPU Worker)

1. **Build and push your Docker image** (do this once from your laptop):
   ```bash
   cd singraal
   # CUDA_ARCH omitted → portable fatbinary (runs on any GPU below)
   docker build -t your-dockerhub/singraal:v13 .
   docker push your-dockerhub/singraal:v13
   ```

2. **Create a RunPod template**:
   - Go to RunPod → Templates → New Template
   - Container image: `your-dockerhub/singraal:v13`
   - Container start command: leave blank (uses `entrypoint.sh`)
   - Environment variables:
     ```
     TARGET_X=<64 hex chars>
     TARGET_Y=<64 hex chars>
     RANGE_BITS=135
     ALL_GPUS=1
     ```
   - Volume: mount `/data` for checkpoint persistence

3. **Deploy pods** — select RTX 4090 or H100, deploy N pods.

### Distributed Pool on RunPod

**Coordinator pod** (1 pod, needs stable IP or use RunPod's network volume):
```
SERVE=1
TARGET_X=<hex>
TARGET_Y=<hex>
BIND=0.0.0.0:5135
```
Expose port 5135 in the pod settings.

**Worker pods** (N pods):
```
TARGET_X=<hex>
TARGET_Y=<hex>
COORDINATOR=<coordinator_pod_ip>:5135
ALL_GPUS=1
```

---

## vast.ai

### Instance Setup

1. Search for instances: filter by GPU (RTX 4090 recommended), CUDA 12.x+
2. Select instance → "Edit instance" → Docker image: `your-dockerhub/singraal:v13`
3. Set environment variables in the "Environment" tab:
   ```bash
   TARGET_X=...
   TARGET_Y=...
   RANGE_BITS=135
   ALL_GPUS=1
   COORDINATOR=<your_coordinator_ip>:5135   # if using pool mode
   ```
4. Open port 5135 on the coordinator instance

### vast.ai CLI (scripted scaling)

```bash
# Install
pip install vastai

# Find cheapest RTX 4090 instances
vastai search offers --storage 10 --gpu-name RTX_4090 --order dph_total

# Launch 10 worker instances
for i in $(seq 1 10); do
  vastai create instance <offer_id> \
    --image your-dockerhub/singraal:v13 \
    --env '-e TARGET_X=... -e TARGET_Y=... -e COORDINATOR=<ip>:5135 -e ALL_GPUS=1' \
    --disk 10
done

# Monitor
vastai show instances
```

---

## Lambda Labs

Lambda Labs provides bare-metal A100/H100 instances — best throughput/price for serious runs.

### Setup

```bash
# SSH into your Lambda instance
ssh ubuntu@<lambda_ip>

# Install Docker + nvidia-container-toolkit
curl -fsSL https://get.docker.com | sh
distribution=$(. /etc/os-release; echo $ID$VERSION_ID)
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
curl -s -L "https://nvidia.github.io/libnvidia-container/$distribution/libnvidia-container.list" \
  | sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' \
  | sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
sudo apt-get update && sudo apt-get install -y nvidia-container-toolkit
sudo systemctl restart docker

# Build — fatbinary runs on any GPU (or pin with --build-arg CUDA_ARCH=sm_80)
git clone https://github.com/afkmoney/singraal
cd singraal
docker build -t singraal:v13 .

# Run standalone (all A100s)
docker run --gpus all \
  -e TARGET_X=<hex> -e TARGET_Y=<hex> \
  -e ALL_GPUS=1 \
  -v /data/singraal:/data \
  singraal:v13
```

### Multi-node pool on Lambda (1 coordinator + N workers)

**On coordinator node:**
```bash
docker run --gpus all -d --restart always -p 5135:5135 \
  -e SERVE=1 -e TARGET_X=<hex> -e TARGET_Y=<hex> \
  -v /data/singraal:/data \
  singraal:v13
```

**On each worker node** (replace `COORD_IP` with coordinator's private IP):
```bash
docker run --gpus all -d --restart always \
  -e TARGET_X=<hex> -e TARGET_Y=<hex> \
  -e COORDINATOR=<COORD_IP>:5135 \
  -e ALL_GPUS=1 \
  -v /data/singraal:/data \
  singraal:v13
```

---

## Local Multi-GPU (docker-compose)

For machines with 4+ GPUs, use the included `docker-compose.yml`:

```bash
cd cloud
export TARGET_X=<hex64>
export TARGET_Y=<hex64>
# CUDA_ARCH optional — leave unset for a portable fatbinary, or pin one arch:
# export CUDA_ARCH=sm_89   # RTX 4090

# Build once
docker compose build

# Start coordinator + 4 workers
docker compose up -d

# Monitor logs
docker compose logs -f coordinator
docker compose logs -f worker-0

# Scale workers (if you have 8 GPUs)
docker compose up -d --scale worker=8
```

---

## Cost Estimates (Bitcoin Puzzle #135)

sinGRAAL v13: C ≈ 0.55 (6-aut + bidir + GLV4D Halton), E[ops] ≈ 0.55 × 2^67.5 / √12 ≈ 2^65.3 steps

| Config | Gstep/s | Temps | Coût estimé |
|--------|---------|-------|-------------|
| 1× RTX 4090 | ~1.5 | ~600 ans | — |
| 8× RTX 4090 (nœud unique) | ~12 | ~75 ans | ~$30/jour |
| 32× RTX 4090 (vast.ai) | ~48 | ~19 ans | ~$120/jour |
| 100× A100 (Lambda) | ~75 | ~12 ans | ~$1,200/jour |
| 1 000× RTX 4090 (farm) | ~1 500 | ~220 jours | ~$3,000/jour |
| 10 000× RTX 4090 | ~15 000 | ~22 jours | ~$30,000/jour |

> Cible : puzzle #135 — `02145d2611c823a396ef6712ce0f712f09b9b4f3135e3e0aa3230fb9b6d08d1e16`
> Récompense ≈ 135 BTC. TARGET_X/TARGET_Y dérivés de cette clé publique compressée.

---

## Checkpoint Persistence

Always mount `/data` as a volume — the coordinator saves the DP table every 60s:

```bash
docker run --gpus all \
  -v /host/path/to/data:/data \
  -e TARGET_X=... \
  singraal:v13
```

On restart, it resumes automatically from the checkpoint.

---

## Monitoring

The coordinator prints progress every 10 seconds:
```
[coord] 12.3M DPs total | 45.2k DP/s | 32 workers | table=8192
```

Workers print:
```
[GPU 0] 4.20B steps (18.3%) | 12450 DPs | 1.03 Gstep/s | 41.2 DP/s | table=8192 | ETA~1850.2d
```
