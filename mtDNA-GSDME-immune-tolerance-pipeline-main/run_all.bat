@echo off
setlocal enabledelayedexpansion

REM ============================================
REM User config
REM ============================================
set CONFIG=config\paths_example.yaml
set LOGDIR=results\logs
set IMAGE_NAME=mtDNA_pipeline:latest

if not exist %LOGDIR% (
    mkdir %LOGDIR%
)

echo ===== mtDNA-GSDME Pipeline Started: %date% %time% ===== > %LOGDIR%\run_all.log

REM ============================================
REM Check if Docker image exists; build if not
REM ============================================
echo [INFO] Checking Docker image... >> %LOGDIR%\run_all.log
docker image inspect %IMAGE_NAME% >NUL 2>&1

if errorlevel 1 (
    echo [INFO] Docker image not found. Building... >> %LOGDIR%\run_all.log
    docker build -t %IMAGE_NAME% .
    if errorlevel 1 (
        echo [ERROR] Docker build failed. Check Dockerfile.
        exit /b 1
    )
)

REM ============================================
REM Run inside container
REM ============================================
echo [INFO] Running pipeline in Docker container... >> %LOGDIR%\run_all.log

docker run --rm -it ^
