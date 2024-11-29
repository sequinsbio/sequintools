#!/bin/bash

set -e

git config --global --add safe.directory /workspaces/glossary

pre-commit install --install-hooks
