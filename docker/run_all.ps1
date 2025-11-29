#!/usr/bin/env pwsh
# ================================
# mtDNA-GSDME Pipeline (PowerShell)
# ================================

# Exit on error
$ErrorActionPreference = "Stop"

# -------------------------
# User Configuration
# -------------------------
$CONFIG      = "config/paths_example.yaml"
$LOGDIR      = "results/logs"
$IMAGE_NAME  = "mtDNA_pipeline:latest"

# Create log directory
if (!(Test-Path $LOGDIR)) {
    New-Item -ItemType Directory -Path $LOGDIR | Out-Null
}

$logFile = "$LOGDIR/run_all.log"
"===== mtDNA-GSDME Pipeline Started: $(Get-Date) =====" | Out-File $logFile

function Log {
    param([string]$msg)
    $msg | Tee-Object -FilePath $logFile -Append
}

# -------------------------
# Check Docker Image
# -------------------------
Log "[INFO] Checking Docker image..."
$inspect = docker image inspect $IMAGE_NAME 2>$null

if ($LASTEXITCODE -ne 0) {
    Log "[INFO] Image not found. Building Docker image..."
    docker build -t $IMAGE_NAME .
    if ($LASTEXITCODE -ne 0) {
        Log "[ERROR] Docker build failed."
        exit 1
    }
}

# -------------------------
# Run Pipeline Inside Container
# -------------------------
Log "[INFO] Starting pipeline in Docker..."

$dockerCmd = @"
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
"@

docker run --rm -it `
    -v "${PWD}:/workspace" `
    -w /workspace `
    $IMAGE_NAME `
    bash -lc $dockerCmd

if ($LASTEXITCODE -ne 0) {
    Log "[ERROR] Pipeline execution failed."
    exit 1
}

# -------------------------
# Finish
# -------------------------
Log "===== Pipeline Finished: $(Get-Date) ====="
Write-Host "`nPipeline completed successfully!" -ForegroundColor Green
exit 0
