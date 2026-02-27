#!/usr/bin/env bash
set -euo pipefail

COMPOSE_FILE=${COMPOSE_FILE:-docker-compose.yml}
RUN_ADAR2=${RUN_ADAR2:-1}

mkdir -p results

docker compose -f "${COMPOSE_FILE}" build analysis

docker compose -f "${COMPOSE_FILE}" run --rm analysis \
  "python -m pytest tests/test_validation.py tests/test_pocket_detection.py"

if [[ "${RUN_ADAR2}" == "1" ]]; then
  docker compose -f "${COMPOSE_FILE}" run --rm analysis \
    "python scripts/phase1_validate_adar2.py"
else
  echo "Skipping ADAR2 validation because RUN_ADAR2=${RUN_ADAR2}"
fi
