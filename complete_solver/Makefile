# sinGRAAL Makefile — Toom-Cook-6 CUDA kangaroo solver

CUDA_ARCH  ?= sm_80
DOCKER_TAG ?= singraal:latest
PORT       ?= 5135

.PHONY: all build release clean docker docker-push run-local coordinator worker

all: build

# ── Rust + CUDA build ──────────────────────────────────────────────────────────

build:
	CUDA_ARCH=$(CUDA_ARCH) cargo build

release:
	CUDA_ARCH=$(CUDA_ARCH) cargo build --release

clean:
	cargo clean

# ── Docker ────────────────────────────────────────────────────────────────────

docker:
	docker build --build-arg CUDA_ARCH=$(CUDA_ARCH) -t $(DOCKER_TAG) .

docker-push:
	docker push $(DOCKER_TAG)

# ── Local run (requires GPU + compiled binary) ────────────────────────────────

run-local: release
	./target/release/singraal --config config.toml

# ── Distributed mode ──────────────────────────────────────────────────────────
# Start coordinator:   make coordinator
# Start worker:        COORD=<host:port> make worker

coordinator: release
	./target/release/singraal --config config.toml --serve 0.0.0.0:$(PORT)

worker: release
ifndef COORD
	$(error COORD is not set — usage: COORD=host:port make worker)
endif
	./target/release/singraal --config config.toml --coordinator $(COORD)

# ── Docker distributed ────────────────────────────────────────────────────────

docker-coordinator: docker
	docker run --gpus all --rm -p $(PORT):$(PORT) $(DOCKER_TAG) \
		singraal --config /singraal/config.toml --serve 0.0.0.0:$(PORT)

docker-worker: docker
ifndef COORD
	$(error COORD is not set — usage: COORD=host:port make docker-worker)
endif
	docker run --gpus all --rm $(DOCKER_TAG) \
		singraal --config /singraal/config.toml --coordinator $(COORD)
