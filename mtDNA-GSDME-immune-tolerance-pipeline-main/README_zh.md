✅ mtDNA–GSDME Immune Tolerance Pipeline

Reproducible Single-cell + Spatial Transcriptomics + Causal Inference Pipeline
for studying mtDNA–GSDME axis and immune tolerance mechanisms

🚀 简介

本仓库包含用于研究 mtDNA–GSDME–免疫耐受机制 的完整分析管线，包括：

单细胞数据（scRNA-seq）的标准化流程

空间转录组（ST）的预处理、解卷积（cell2location / RCTD / DWLS）

GSDME 相关基因集打分

mtDNA 相关代谢通路分析

免疫逃逸与免疫耐受指标构建

因果推断（DoWhy / CausalML / EconML）

可重复运行的 Docker 镜像 + Snakemake 流程

完整 notebooks 用于展示关键结果

本项目完全容器化，可复现性强，支持本地、服务器、HPC 以及 GitHub Actions 自动构建与测试。

📁 目录结构
mtDNA-GSDME-immune-tolerance-pipeline/
│
├── scripts/                 # 主分析脚本 (Python + R)
│   ├── 00_download_data.R
│   ├── 01_preprocess_scRNA.py
│   ├── 02_run_seurat.R
│   ├── 03_run_cell2location.py
│   ├── 04_run_RCTD.R
│   ├── 05_run_DWLS.py / .R
│   ├── 06_causal_inference.py
│   └── utils/
│
├── notebooks/               # 分析可视化 Jupyter notebooks
│   ├── scRNA_visualization.ipynb
│   ├── spatial_mapping.ipynb
│   └── causal_plots.ipynb
│
├── config/
│   ├── paths_example.yaml   # 示例配置文件
│   └── params.yaml          # 参数配置
│
├── docker/
│   ├── Dockerfile
│   └── run_all.sh / .bat / .ps1
│
├── results/                 # 结果输出（默认被 .gitignore 忽略）
│
├── requirements.txt         # Python 环境
├── .dockerignore            # Docker 构建优化
├── .gitignore
│
└── README.md                # ← 当前文件

🐳 使用 Docker（推荐）
1. 构建镜像

⚠ Docker tag 必须为小写。

docker build -t mtdna_pipeline:latest .

2. 运行容器

挂载当前目录为 /workspace（内部由脚本自动识别）。

docker run -it -v "$(pwd)":/workspace mtdna_pipeline:latest bash

3. 在容器内运行 Pipeline

如使用 Snakemake：

snakemake --cores 8


或者直接运行整套脚本：

bash docker/run_all.sh


Windows PowerShell：

.\docker\run_all.ps1

🛠 软件依赖
Python（自动安装自 Docker）

Python ≥ 3.10

numpy / pandas / scipy

scanpy / anndata / squidpy

torch（CPU）

pyro-ppl / cell2location

tangram-sc

causalml / dowhy / econml

R（自动安装自 Docker）

Seurat / SeuratObject / SeuratDisk

GSVA / clusterProfiler / enrichplot

SingleCellExperiment / scran / scater / batchelor

RCTD / Giotto

DWLS：使用 Python 或 R 版本（视仓库可用性）

所有 R/Python 包都在 Dockerfile 中已自动安装。

⚙ 配置文件

请将示例文件复制到你的配置目录：

cp config/paths_example.yaml config/paths.yaml


编辑关键路径：

base_dir: "/workspace"
raw_data_dir: "data_raw"
processed_data_dir: "data_processed"
results_dir: "results"
n_workers: 8
memory_gb: 64

📦 数据准备
自动下载

执行：

Rscript scripts/00_download_data.R


若需要 token，请将其放置在：

~/.mtDNA_pipeline/credentials.yaml


格式：

gdc_token: "XXXX"

🔬 分析步骤（概要）
Step 1 — 单细胞预处理
python scripts/01_preprocess_scRNA.py

Step 2 — Seurat 标准流程
Rscript scripts/02_run_seurat.R

Step 3 — 空间解卷积

cell2location（Python）

RCTD（R）

DWLS（Python 或 R）

Step 4 — mtDNA & GSDME 通路分析

GSVA、meta-score、clusterProfiler。

Step 5 — 免疫耐受评分

TIDE-like scoring, exhaustion signatures。

Step 6 — 因果推断

DoWhy + EconML + CausalML。

🧪 测试（CI 支持）

你可以运行轻量测试：

pytest -q


CI 会自动进行：

Docker 构建

基础脚本运行

轻量数据检查

📜 License

建议选择 MIT / Apache-2.0，如你决定我可以自动生成 LICENSE 文件。

🤝 贡献指南

欢迎提交 Issue 和 PR。
建议遵守：

使用 feature/* 分支提交功能

在 PR 中附测试与说明

Notebook 放在 notebooks/ 中

如需要，我可以自动生成一份 modern CONTRIBUTING.md。

📧 联系方式

如需讨论管线扩展、Docker 优化、大规模 HPC 适配，请联系仓库作者。
