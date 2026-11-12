#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
rm -rf "$ROOT/demo/results"
python3 "$ROOT/demo/make_demo.py"
bash "$ROOT/demo/build_reference.sh"
bash "$ROOT/scripts/run_pipeline.sh" \
  --config "$ROOT/demo/config.demo.sh" \
  --start 1 --end 11
