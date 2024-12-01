#!/bin/bash

set -e

git config --global --add safe.directory /workspaces/sequintools

pre-commit install --install-hooks
