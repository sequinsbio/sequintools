
FROM mcr.microsoft.com/devcontainers/python:3.12-bookworm

RUN apt update \
  && apt install -y --no-install-recommends \
  cmake \
  clang \
  && rm -rf /var/lib/apt/lists/*
