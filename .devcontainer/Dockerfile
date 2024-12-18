
FROM mcr.microsoft.com/devcontainers/python:3.12-bookworm

RUN apt update \
  && apt install -y --no-install-recommends \
  cmake \
  clang \
  curl \
  build-essential \
  && rm -rf /var/lib/apt/lists/*

# This has been moved into the Dockerfile to try and speed up testing
# in github, which uses the devcontainer to run the tests.

WORKDIR /code

COPY ./requirements.txt /code/requirements.txt

RUN pip install --upgrade pip

# Required to do manual pushes to codecov from inside our container.
RUN pip install --user -r /code/requirements.txt

# Install rust as user vscode to avoid caching directory write issues when
# running the cargo commands.
USER vscode
WORKDIR /home/vscode

# Install Rust using rustup
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | \
  sh -s -- -y --no-modify-path --default-toolchain stable

# Set environment variables for all users
ENV CARGO_HOME=/home/vscode/.cargo \
  RUSTUP_HOME=/home/vscode/.rustup \
  PATH="/home/vscode/.cargo/bin:${PATH}"

# Install Rust components
RUN cargo install cargo-tarpaulin@0.31.3
