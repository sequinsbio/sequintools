# Use multi-stage builds to reduce the size of the final image
# https://docs.docker.com/build/building/multi-stage/#use-multi-stage-builds

# Stage 1: Build the rust app
FROM mcr.microsoft.com/devcontainers/rust:1-1-bookworm AS builder

WORKDIR /usr/sequins

# Copy the Cargo.toml and Cargo.lock files first
COPY ./Cargo.toml ./Cargo.lock ./

COPY ./src ./src

# Combine commands with && to reduce the number of layers.
RUN apt-get update && apt-get install -y --no-install-recommends \
  cmake \
  clang \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && cargo build --release

# Stage 2: Create final Docker image with debian-slim.
FROM ghcr.io/linuxcontainers/debian-slim:latest

# Copy app and test data.
COPY --from=builder /usr/sequins/target/release/sequintools /usr/local/bin/sequintools

CMD [ "sequintools" ]
