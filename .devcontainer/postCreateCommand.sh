
#!/bin/bash

set -e -x

# Ensure we have the latest version of pip
pip install --upgrade pip

# Required to do manual pushes to codecov from inside our container.
pip install --user -r ./.devcontainer/requirements.txt

# Install container dependencies.
cargo install cargo-tarpaulin cargo-llvm-cov
