# sinGRAAL Makefile — Toom-Cook-6 CUDA kangaroo solver
#
# GPU architecture is auto-detected by build.rs:
#   - CUDA_ARCH unset → portable multi-arch fatbinary (runs on any GPU)
#   - CUDA_ARCH=sm_XX → single-arch build (faster compile)
# Set CUDA_ARCH only if you want to pin one architecture.

CUDA_ARCH  ?=
DOCKER_TAG ?= singraal:latest
PORT       ?= 5135

# Puzzle target (override on the command line or via env).
TARGET_X    ?=
TARGET_Y    ?=
RANGE_BITS  ?= 135
NUM_ANIMALS ?= 262144

export CUDA_ARCH

.PHONY: all build release clean docker docker-push run-local coordinator worker \
        docker-coordinator docker-worker

all: release

# ── Rust + CUDA build ──────────────────────────────────────────────────────────

build:
	cargo build --features cuda

release:
	cargo build --release --features cuda

clean:
	cargo clean

# ── Local run (requires GPU + compiled binary) ────────────────────────────────
# Usage: make run-local TARGET_X=<hex> TARGET_Y=<hex>

run-local: release
	@test -n "$(TARGET_X)" || { echo "TARGET_X required"; exit 1; }
	@test -n "$(TARGET_Y)" || { echo "TARGET_Y required"; exit 1; }
	./target/release/singraal --target-x $(TARGET_X) --target-y $(TARGET_Y) \
		--range-bits $(RANGE_BITS) --all-gpus --num-animals $(NUM_ANIMALS)

# ── Distributed mode ──────────────────────────────────────────────────────────
# Coordinator:  make coordinator TARGET_X=<hex> TARGET_Y=<hex>
# Worker:       make worker COORD=<host:port> TARGET_X=<hex> TARGET_Y=<hex>

coordinator: release
	@test -n "$(TARGET_X)" || { echo "TARGET_X required"; exit 1; }
	@test -n "$(TARGET_Y)" || { echo "TARGET_Y required"; exit 1; }
	./target/release/singraal --target-x $(TARGET_X) --target-y $(TARGET_Y) \
		--range-bits $(RANGE_BITS) --serve --bind 0.0.0.0:$(PORT)

worker: release
ifndef COORD
	$(error COORD is not set — usage: COORD=host:port make worker TARGET_X=.. TARGET_Y=..)
endif
	@test -n "$(TARGET_X)" || { echo "TARGET_X required"; exit 1; }
	@test -n "$(TARGET_Y)" || { echo "TARGET_Y required"; exit 1; }
	./target/release/singraal --target-x $(TARGET_X) --target-y $(TARGET_Y) \
		--range-bits $(RANGE_BITS) --coordinator $(COORD) --all-gpus \
		--num-animals $(NUM_ANIMALS)

# ── Docker ────────────────────────────────────────────────────────────────────
# Default build = portable fatbinary. Pin an arch with: make docker CUDA_ARCH=sm_89

docker:
	docker build --build-arg CUDA_ARCH=$(CUDA_ARCH) -t $(DOCKER_TAG) .

docker-push:
	docker push $(DOCKER_TAG)

docker-coordinator: docker
	docker run --gpus all --rm -p $(PORT):$(PORT) \
		-e TARGET_X=$(TARGET_X) -e TARGET_Y=$(TARGET_Y) \
		-e RANGE_BITS=$(RANGE_BITS) -e SERVE=1 -e BIND=0.0.0.0:$(PORT) \
		$(DOCKER_TAG)

docker-worker: docker
ifndef COORD
	$(error COORD is not set — usage: COORD=host:port make docker-worker TARGET_X=.. TARGET_Y=..)
endif
	docker run --gpus all --rm \
		-e TARGET_X=$(TARGET_X) -e TARGET_Y=$(TARGET_Y) \
		-e RANGE_BITS=$(RANGE_BITS) -e COORDINATOR=$(COORD) \
		$(DOCKER_TAG)
