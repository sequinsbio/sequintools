
FROM mcr.microsoft.com/devcontainers/python:3.12-bookworm

RUN apt update \
  && apt install -y --no-install-recommends \
  build-essential \
  cmake \
  clang \
  curl \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /code

COPY ./requirements.txt /code/requirements.txt

RUN pip install --upgrade pip \
  && pip install -r /code/requirements.txt

# Install Rust using rustup
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

# Add Rust to PATH
ENV PATH="/root/.cargo/bin:${PATH}"

RUN cargo install cargo-tarpaulin@0.31.3
