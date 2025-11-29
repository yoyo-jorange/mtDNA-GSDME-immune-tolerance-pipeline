#!/usr/bin/env bash
set -euo pipefail

CONFIG=config/paths_example.yaml
LOGDIR=results/logs
mkdir -p "$LOGDIR"

echo "Start pipeline: $(date)" | tee "$LOGDIR/run_all.log"

# Build docker image (only if not present)
docker image inspect mtDNA_pipeline:latest >/dev/null 2>&1 || \
  docker build -t mtDNA_pipeline:latest docker/

# Run containerized pipeline (bind mount base dir)
docker run --rm -it \
  -v "$(pwd)":/workspace \
  -v /tmp:/tmp \
  -e CONFIG_PATH=/workspace/$CONFIG \
  mtDNA_pipeline:latest \
  bash -lc "cd /workspace && python scripts/00_check_environment.py --config $CONFIG 2>&1 | tee results/logs/00_check_environment.log && \
           python scripts/01_download_data.py --config $CONFIG 2>&1 | tee results/logs/01_download_data.log && \
           python scripts/02_preprocess_spatial.py --config $CONFIG 2>&1 | tee results/logs/02_preprocess_spatial.log && \
           python scripts/03_preprocess_scRNA.py --config $CONFIG 2>&1 | tee results/logs/03_preprocess_scRNA.log && \
           Rscript scripts/20_spatial_deconvolution.R $CONFIG 2>&1 | tee results/logs/20_spatial_deconvolution.log && \
           jupyter nbconvert --to html notebooks/A_quicklook_results.ipynb --output results/reports/A_quicklook_results.html"
