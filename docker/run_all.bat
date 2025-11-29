#!/usr/bin/env bash
set -e

# 1. 构建或拉取镜像
IMAGE_NAME="mtdna-gsdme-pipeline:v1"

echo "=== 构建 Docker 镜像（仅首次需要，时间稍长） ==="
docker build -t ${IMAGE_NAME} docker

# 2. 运行容器并在其中依次执行脚本
echo "=== 启动容器并运行完整流程 ==="
docker run --rm -it \
  -v "$(pwd)":/workspace \
  ${IMAGE_NAME} \
  bash -c "
    cd /workspace && \
    python scripts/00_check_environment.py && \
    python scripts/01_download_data.py && \
    python scripts/02_preprocess_spatial.py && \
    python scripts/03_preprocess_scRNA.py && \
    python scripts/04_preprocess_bulk_cfDNA.py && \
    python scripts/10_compute_mtLeakScore.py && \
    python scripts/11_compute_nonCanonicalPyroScore.py
    # 后续可以继续串联 R 脚本和 Notebook
  "
