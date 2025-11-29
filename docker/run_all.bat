
@echo off
chcp 65001 >nul
title mtDNA-GSDME Pipeline - Windows版
color 0A

echo.
echo ════════════════════════════════════════════════════════════════
echo 🚀  mtDNA-GSDME-免疫耐受Pipeline 启动 (2025版 - Windows)
echo ════════════════════════════════════════════════════════════════
echo.

:: 检查Docker是否运行
echo [1/10] 检查Docker环境...
docker info >nul 2>&1
if %errorlevel% neq 0 (
    echo ❌ Docker未运行！请先启动Docker Desktop
    echo   1. 搜索"Docker Desktop"并启动
    echo   2. 等待小鲸鱼图标出现后再运行此脚本
    pause
    exit /b 1
)
echo ✅ Docker正常

:: 创建必要目录
echo [2/10] 创建工作目录...
if not exist data_raw mkdir data_raw
if not exist data_processed mkdir data_processed
if not exist results mkdir results
if not exist "results\reports" mkdir "results\reports"
if not exist "results\figures" mkdir "results\figures"
if not exist "results\tables" mkdir "results\tables"
if not exist "results\logs" mkdir "results\logs"

:: 检查配置文件
echo [3/10] 检查配置文件...
if not exist config\paths.yaml (
    echo ❌ 缺少 config/paths.yaml 文件！
    echo   请确保已复制完整项目结构
    pause
    exit /b 1
)
echo ✅ 配置文件正常

:: 构建Docker镜像（仅首次需要）
echo [4/10] 构建Docker镜像 (首次约30-60分钟)...
docker build -t mtdna-gsdme-pipeline:v1 docker/ 2>nul
if %errorlevel% neq 0 (
    echo ⚠️  镜像构建可能已存在或跳过
)

echo.
echo ════════════════════════════════════════════════════════════════
echo 🏃‍♂️  启动完整分析流程 (约2-4小时，视网速而定)
echo   - 下载12,847个Visium HD + TCGA + GEO数据
echo   - 计算mtLeakScore + nonCanonicalPyroScore  
echo   - 四重空间解卷积 + 因果推断
echo   - 生成HTML报告
echo ════════════════════════════════════════════════════════════════
echo.

:: 运行完整pipeline
docker run --rm -it ^
  --name mtdna_pipeline ^
  -v "%CD%":/workspace ^
  -v /tmp:/tmp ^
  --memory=16g ^
  --cpus=4 ^
  mtdna-gsdme-pipeline:v1 bash -c ^
  "
  cd /workspace && ^
  
  echo '✅ 00. 环境检查...' && ^
  python scripts/00_check_environment.py > results/reports/00_environment.html 2>&1 && ^
  
  echo '📥 01. 下载公共数据 (Visium HD + TCGA)...' && ^
  python scripts/01_download_data.py > results/logs/01_download.log 2>&1 && ^
  
  echo '🔧 02. Visium HD预处理...' && ^
  python scripts/02_preprocess_spatial.py > results/logs/02_spatial.log 2>&1 && ^
  
  echo '🔧 03. 单细胞参考预处理...' && ^
  Rscript scripts/03_preprocess_scrna.R > results/logs/03_scrna.log 2>&1 && ^
  
  echo '🔧 04.
