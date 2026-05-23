# sinGRAAL v12 — Cloud GPU Deployment Guide

Deploy sinGRAAL across cloud GPU providers in minutes using Docker.

---

## GPU Architecture Selection

Build the Docker image with the correct `CUDA_ARCH` for your GPU:

| GPU | CUDA_ARCH | Cloud Providers |
|-----|-----------|-----------------|
| H100 | `sm_90` | RunPod, Lambda Labs, CoreWeave |
| A100 | `sm_80` | RunPod, Lambda Labs, vast.ai, AWS p4d |
| RTX 4090 | `sm_89` | RunPod, vast.ai, Jarvis Labs |
| RTX 3090 / 3080 Ti | `sm_86` | vast.ai, RunPod |
| RTX 2080 Ti | `sm_75` | vast.ai (budget) |
| V100 | `sm_70` | AWS p3, Lambda Labs |

---

## RunPod

### Quick Start (Single GPU Worker)

1. **Build and push your Docker image** (do this once from your laptop):
   ```bash
   cd kangaroo
   docker build -t your-dockerhub/singraal:v12-4090 --build-arg CUDA_ARCH=sm_89 .
   docker push your-dockerhub/singraal:v12-4090
   ```

2. **Create a RunPod template**:
   - Go to RunPod → Templates → New Template
   - Container image: `your-dockerhub/singraal:v12-4090`
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
2. Select instance → "Edit instance" → Docker image: `your-dockerhub/singraal:v12-4090`
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
    --image your-dockerhub/singraal:v12-4090 \
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

# Build (A100 = sm_80, H100 = sm_90)
git clone https://github.com/afkmoney/singraal
cd singraal/kangaroo
docker build -t singraal:v12 --build-arg CUDA_ARCH=sm_80 .

# Run standalone (all A100s)
docker run --gpus all \
  -e TARGET_X=<hex> -e TARGET_Y=<hex> \
  -e ALL_GPUS=1 \
  -v /data/singraal:/data \
  singraal:v12
```

### Multi-node pool on Lambda (1 coordinator + N workers)

**On coordinator node:**
```bash
docker run --gpus all -d --restart always -p 5135:5135 \
  -e SERVE=1 -e TARGET_X=<hex> -e TARGET_Y=<hex> \
  -v /data/singraal:/data \
  singraal:v12
```

**On each worker node** (replace `COORD_IP` with coordinator's private IP):
```bash
docker run --gpus all -d --restart always \
  -e TARGET_X=<hex> -e TARGET_Y=<hex> \
  -e COORDINATOR=<COORD_IP>:5135 \
  -e ALL_GPUS=1 \
  -v /data/singraal:/data \
  singraal:v12
```

---

## Local Multi-GPU (docker-compose)

For machines with 4+ GPUs, use the included `docker-compose.yml`:

```bash
cd cloud
export TARGET_X=<hex64>
export TARGET_Y=<hex64>
export CUDA_ARCH=sm_89   # RTX 4090

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

sinGRAAL v12: C ≈ 1.10, E[ops] ≈ 1.10 × 2^67.5 ≈ 1.94 × 10^20 steps

| Config | Gstep/s | Time | Est. Cost |
|--------|---------|------|-----------|
| 1 RTX 4090 | ~1.0 | ~2,050 years | — |
| 8 RTX 4090 (single node) | ~8.0 | ~256 years | ~$30/day |
| 32 RTX 4090 (vast.ai) | ~32 | ~64 years | ~$120/day |
| 100 A100 (Lambda) | ~50 | ~40 years | ~$1,200/day |
| 1,000 RTX 4090 (farm) | ~1,000 | ~2 years | ~$3,000/day |
| 10,000 RTX 4090 | ~10,000 | ~75 days | ~$30,000/day |

> Note: puzzle #135 reward ≈ 135 BTC. GPU farm economics only positive if BTC price rises enough or the pool gets lucky (expected value depends on probability of early solve).

---

## Checkpoint Persistence

Always mount `/data` as a volume — the coordinator saves the DP table every 60s:

```bash
docker run --gpus all \
  -v /host/path/to/data:/data \
  -e TARGET_X=... \
  singraal:v12
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
