#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
# --clear recreates only this project's virtual environment.
/usr/bin/python3.11 -m venv --clear .venv
.venv/bin/python --version
.venv/bin/pip install -r requirements.txt
