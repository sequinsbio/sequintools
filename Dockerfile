# lukemathwalker/cargo-chef:0.1.73-rust-1.93.0
FROM lukemathwalker/cargo-chef@sha256:9e564338159930326c71b6c3df9ed5e5b60b0437b84861c6bb6767b1964a0fd4 AS chef
WORKDIR /app

FROM chef AS planner
COPY ./Cargo.toml ./Cargo.toml
COPY ./Cargo.lock ./Cargo.lock
COPY ./build.rs ./build.rs
COPY ./src ./src
RUN cargo chef prepare --recipe-path recipe.json

# Build the rust app
FROM chef AS builder
RUN apt-get update && apt-get install -y --no-install-recommends \
  cmake=3.31.6-2 \
  clang=1:19.0-63 \
  && rm -rf /var/lib/apt/lists/*
COPY --from=planner /app/recipe.json recipe.json
# Accept version string from build context (fallback handled in build.rs)
ARG SEQUINTOOLS_GIT_VERSION
ENV SEQUINTOOLS_GIT_VERSION=${SEQUINTOOLS_GIT_VERSION}

RUN cargo chef cook --release --recipe-path recipe.json
COPY ./Cargo.toml ./Cargo.toml
COPY ./Cargo.lock ./Cargo.lock
COPY ./build.rs ./build.rs
COPY ./src ./src
# Build with locked deps; build.rs will embed version using
# SEQUINTOOLS_GIT_VERSION if provided
RUN cargo build --locked --release

# Create final Docker image with debian-slim.
FROM debian:13.3-slim

# Augment debian-slim with tools needed to run in nextflow
RUN apt-get update && \
  apt-get --no-install-recommends install -y \
  procps=2:4.0.4-9 \
  samtools=1.21-1 && \
  apt-get clean && \
  apt-get autoremove -y && \
  rm -rf /var/lib/apt/lists/*

# Copy app and test data.
COPY --from=builder /app/target/release/sequintools /usr/local/bin/sequintools

CMD [ "/usr/local/bin/sequintools" ]
