#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# User config
# -----------------------------
CONFIG="config/paths_example.yaml"
LOGDIR="results/logs"
mkdir -p "$LOGDIR"

echo "===== mtDNA-GSDME Pipeline Started: $(date) =====" | tee "$LOGDIR/run_all.log"

# -----------------------------
# Build image if not exists
# -----------------------------
IMAGE_NAME="mtDNA_pipeline:latest"

if ! docker image inspect $IMAGE_NAME >/dev/null 2>&1; then
    echo "[INFO] Docker image not found. Building..." | tee -a "$LOGDIR/run_all.log"
    docker build -t $IMAGE_NAME .
fi

# -----------------------------
# Run inside container
# -----------------------------
docker run --rm -it \
    -v "$(pwd)":/workspace \
    -w /workspace \
    $IMAGE_NAME \
    bash -lc "
        echo '[INFO] Running pipeline inside container...'
        CONFIG_PATH=/workspace/$CONFIG

        python scripts/00_check_environment.py --config \$CONFIG_PATH 2>&1 | tee results/logs/00_check_environment.log

        python scripts/01_download_data.py --config \$CONFIG_PATH 2>&1 | tee results/logs/01_download_data.log

        python scripts/02_preprocess_spatial.py --config \$CONFIG_PATH 2>&1 | tee results/logs/02_preprocess_spatial.log

        python scripts/03_preprocess_scRNA.py --config \$CONFIG_PATH 2>&1 | tee results/logs/03_preprocess_scRNA.log

        python scripts/10_compute_mtLeakScore.py --config \$CONFIG_PATH 2>&1 | tee results/logs/10_compute_mtLeakScore.log

        python scripts/11_compute_nonCanonicalPyroScore.py --config \$CONFIG_PATH 2>&1 | tee results/logs/11_compute_nonCanonicalPyroScore.log

        Rscript scripts/20_spatial_deconvolution.R \$CONFIG_PATH 2>&1 | tee results/logs/20_spatial_deconvolution.log

        jupyter nbconvert --to html notebooks/A_quicklook_results.ipynb \
            --output results/reports/A_quicklook_results.html 2>&1 | tee results/logs/report.log
    "

echo "===== Pipeline Finished: $(date) =====" | tee -a "$LOGDIR/run_all.log"
