#!/bin/bash
set -euo pipefail

echo "🚀 mtDNA-GSDME-免疫耐受pipeline启动 (2025版)"

# 加载配置
CONFIG="config/paths.yaml"
if [[ ! -f "$CONFIG" ]]; then
    echo "❌ 找不到配置文件: $CONFIG"
    exit 1
fi

# 创建目录
mkdir -p data_raw data_processed/{spatial,scrna,bulk,cfdna} results/{reports,figures,tables,logs}

# 构建镜像（仅首次）
IMAGE="mtdna-gsdme-pipeline:v1"
if [[ "$(docker images -q $IMAGE 2> /dev/null || true)" == "" ]]; then
    echo "🐳 构建Docker镜像 (首次约30-60分钟)..."
    docker build -t $IMAGE docker/
fi

# 运行完整流程
docker run --rm -it \
    --name mtdna_pipeline \
    -v $(pwd):/workspace \
    -v /tmp:/tmp \
    $IMAGE bash -c "
    cd /workspace &&

    # 00. 环境检查
    echo '✅ 00. 检查环境...'
    python scripts/00_check_environment.py > results/reports/00_environment.html 2>&1 &&

    # 01. 下载数据 (12,847 Visium HD + TCGA + GEO等)
    echo '📥 01. 下载公共数据...'
    python scripts/01_download_data.py > results/logs/01_download.log 2>&1 &&

    # 02-04. 预处理
    echo '🔧 02. 预处理Visium HD空间数据...'
    python scripts/02_preprocess_spatial.py > results/logs/02_spatial.log 2>&1 &&
    echo '🔧 03. 预处理单细胞数据...'
    Rscript scripts/03_preprocess_scrna.R > results/logs/03_scrna.log 2>&1 &&
    echo '🔧 04. 预处理TCGA/cfDNA...'
    python scripts/04_preprocess_bulk.R > results/logs/04_bulk.log 2>&1 &&

    # 10-11. 核心评分
    echo '⚡ 10. 计算mtLeakScore...'
    python scripts/10_compute_mtLeakScore.py > results/logs/10_mtLeak.log 2>&1 &&
    echo '⚡ 11. 计算nonCanonicalPyroScore...'
    python scripts/11_compute_nonCanonicalPyroScore.py > results/logs/11_pyro.log 2>&1 &&

    # 20. 四重空间解卷积
    echo '🧬 20. 空间解卷积 (Cell2location+Tangram+SpatialDWLS+RCTD)...'
    Rscript scripts/20_spatial_deconvolution.R > results/logs/20_deconv.log 2>&1 &&

    # 21. 空间模式分析
    echo '📊 21. 空间相关性分析...'
    Rscript scripts/21_spatial_correlation.R > results/logs/21_spatial.log 2>&1 &&

    # 30. 因果推断
    echo '🔗 30. CausalML + DoWhy因果推断...'
    jupyter nbconvert --to notebook --execute scripts/30_causal_inference.ipynb --output-dir=results/reports/ &&

    # 40. 泛癌验证
    echo '🌍 40. TCGA泛癌验证 + cfDNA外推...'
    Rscript scripts/40_pan_cancer_validation.R > results/logs/40_pan_cancer.log 2>&1 &&

    # 生成最终报告
    echo '📋 生成综合报告...'
    Rscript scripts/50_final_report.R &&

    echo '🎉 全部流程完成！结果在 results/ 文件夹'
    echo '主要报告: results/reports/final_summary.html'
"

echo "✅ 流程结束！请查看 results/reports/final_summary.html"
