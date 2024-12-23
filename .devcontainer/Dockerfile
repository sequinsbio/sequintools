
FROM mcr.microsoft.com/devcontainers/rust:1-1-bookworm

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  cmake=3.25.1-1 \
  clang=1:14.0-55.7~deb12u1 \
  && rm -rf /var/lib/apt/lists/*
