#!/usr/bin/env bash
set -euo pipefail

IMAGE_NAME=${IMAGE_NAME:-cryptic-ip}
IMAGE_TAG=${IMAGE_TAG:-latest}
PLATFORMS=${PLATFORMS:-linux/amd64,linux/arm64}

if [[ "${1:-}" == "--multi-arch" ]]; then
  docker buildx build \
    --platform "${PLATFORMS}" \
    -t "${IMAGE_NAME}:${IMAGE_TAG}" \
    --push \
    .
else
  docker build \
    -t "${IMAGE_NAME}:${IMAGE_TAG}" \
    .
fi
