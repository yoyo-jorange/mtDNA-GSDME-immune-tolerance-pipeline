# 使用最新的 Rocker 镜像作为基础，支持 R 4.4+ 和 Python 3.11+
FROM rocker/tidyverse:4.4.1

# 设置非交互式安装，避免卡住
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Shanghai

# 安装系统依赖
RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    ghostscript \
    texlive-latex-recommended \
    texlive-fonts-recommended \
    texinfo \
    lmodern \
    default-jdk \
    python3-dev \
    python3-pip \
    curl \
    wget \
    git \
    nano \
    htop \
    && rm -rf /var/lib/apt/lists/*

# 设置 Python 别名
RUN ln -sf /usr/bin/python3 /usr/bin/python

# 安装核心 Python 包（单细胞 + 空间组学 + 因果推断）
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    numpy==1.21.6 \
    pandas==1.4.2 \
    scipy==1.8.0 \
    matplotlib==3.5.3 \
    seaborn==0.11.2

# 🔥 第二层：单细胞核心（官方兼容组合）
RUN pip install --no-cache-dir \
    scanpy==1.9.1 \
    anndata==0.8.0

# 🔥 第三层：空间分析（Tangram官方推荐）
RUN pip install --no-cache-dir \
    squidpy==1.2.2 \
    tangram-sc==1.0.2

# 🔥 第四层：其他分析包
RUN pip install --no-cache-dir \
    harmonypy==0.0.10 \
    gseapy==1.1.3 \
    GEOparse==2.0.4 \
    requests==2.32.3 \
    beautifulsoup4==4.12.3 \
    plotly==5.22.0 \
    jupyterlab==4.2.5 \
    dowhy==0.11.1 \
    lifelines==0.27.8

# 安装核心 R 包（Seurat + 空间解卷积 + GSVA）
RUN R -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages(c('remotes', 'BiocManager', 'devtools'), ask = FALSE)" && \
    R -e "BiocManager::install(c('SingleCellExperiment', 'SingleR', 'scater', \
                               'scran', 'batchelor', 'DoubletFinder'), \
                               ask = FALSE, update = FALSE)" && \
    R -e "install.packages(c('Seurat', 'SeuratObject', 'SeuratDisk', \
                            'ggplot2', 'ggpubr', 'patchwork', 'harmony', \
                            'DoubletFinder', 'celldex', 'scDblFinder', \
                            'data.table', 'Matrix', 'reticulate', 'future', \
                            'future.apply'), repos = 'https://cloud.r-project.org', \
                            dependencies = TRUE)" && \
    # 空间解卷积四重验证
    R -e "remotes::install_github(c('bayesLab/Cell2location@v2.0.3', \
                                   'lzheimer/Giotto@3.1.0', \
                                   'dmcable/spatialDWLS', \
                                   'dmcable/RCTD'), \
                                   dependencies = TRUE, upgrade = FALSE)" && \
    # GSVA 和基因集分析
    R -e "BiocManager::install(c('GSVA', 'clusterProfiler', 'enrichplot', \
                                'org.Hs.eg.db'), ask = FALSE)" && \
    # 因果推断 R 接口
    R -e "install.packages(c('doWhy', 'CausalML'), repos = 'https://cloud.r-project.org')"

# 验证关键包安装
RUN R -e "library(Seurat); library(Cell2location); library(tangram); cat('Seurat:', packageVersion('Seurat'), '\n')" && \
    python3 -c "import scanpy, anndata, dowhy, causalml; print('Scanpy:', scanpy.__version__)"

# 设置工作目录
WORKDIR /workspace

# 复制项目文件
COPY . /workspace

# 默认命令
CMD ["/bin/bash"]
