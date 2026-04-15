#!/usr/bin/env bash
# Run a command inside a linux/amd64 container with MEME 5.1.0 pinned.
# Used on ARM Macs where bioconda has no meme=5.1.0 build for any darwin platform.
#
# Examples:
#   scripts/run-in-linux.sh fimo --version
#   scripts/run-in-linux.sh snakemake --configfile config/default_config.yaml \
#     --until step10_promoter -j 4
#   scripts/run-in-linux.sh bash            # interactive shell
set -euo pipefail

IMAGE_TAG="${GFF_MEME_IMAGE:-genefamilyflow/meme510:latest}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOCKERFILE="${REPO_ROOT}/docker/Dockerfile.meme510"

if ! command -v docker >/dev/null 2>&1; then
    echo "error: docker not found on PATH. Install Docker Desktop and retry." >&2
    exit 127
fi

if ! docker image inspect "${IMAGE_TAG}" >/dev/null 2>&1; then
    echo "[run-in-linux] image ${IMAGE_TAG} not found; building..." >&2
    docker build \
        --platform linux/amd64 \
        -f "${DOCKERFILE}" \
        -t "${IMAGE_TAG}" \
        "${REPO_ROOT}"
fi

if [[ $# -eq 0 ]]; then
    set -- bash
fi

exec docker run --rm -it \
    --platform linux/amd64 \
    -v "${REPO_ROOT}:/work" \
    -w /work \
    "${IMAGE_TAG}" \
    "$@"
