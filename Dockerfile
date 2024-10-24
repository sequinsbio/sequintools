# Use multi-stage builds to reduce the size of the final image
# https://docs.docker.com/build/building/multi-stage/#use-multi-stage-builds

# Stage 1: Build the rust app
FROM mcr.microsoft.com/devcontainers/rust:1-1-bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
  cmake \
  clang \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/sequins

# Copy the Cargo.toml and Cargo.lock files first
COPY ./Cargo.toml ./Cargo.toml
COPY ./Cargo.lock ./Cargo.lock

COPY ./src ./src

RUN cargo build --release

# Stage 2: Create final Docker image
FROM ghcr.io/linuxcontainers/debian-slim:latest

# Copy app and test data.
COPY --from=builder /usr/sequins/target/release/sequintools /usr/local/bin/sequintools

COPY ./example /usr/local/example
COPY ./testdata /usr/local/testdata

CMD [ "sequintools" ]
