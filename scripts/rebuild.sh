#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p results
rebuild_start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
trap 'rebuild_exit=$?; printf "{\"started_utc\":\"%s\",\"finished_utc\":\"%s\",\"exit_code\":%s,\"procedure\":\"scripts/rebuild.sh\"}\n" "$rebuild_start" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$rebuild_exit" > results/build_status.json' EXIT
if [[ "${1:-}" == "--fresh-env" ]]; then bash scripts/setup.sh; fi
.venv/bin/python --version
make all
