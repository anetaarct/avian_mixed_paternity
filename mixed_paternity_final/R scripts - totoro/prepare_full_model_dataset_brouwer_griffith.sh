#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p logs

Rscript "R scripts - totoro/full_model_brouwer_griffith.R" --dataset-only \
  > "logs/full_model_brouwer_griffith_dataset_run.log" 2>&1
