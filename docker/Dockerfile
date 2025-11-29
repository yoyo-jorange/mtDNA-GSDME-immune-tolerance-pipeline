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
RUN pip3 install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir \
    numpy==1.26.4 \
    pandas==2.2.2 \
    scipy==1.13.1 \
    scikit-learn==1.5.1 \
    matplotlib==3.9.2 \
    seaborn==0.13.2 \
    plotly==5.22.0 \
    scanpy==1.10.2 \
    anndata==0.11.3 \
    squidpy==1.5.0 \
    harmonypy==0.0.10 \
    gseapy==1.1.3 \
    lifelines==0.27.8 \
    statsmodels==0.14.2 \
    jupyterlab==4.2.5 \
    nbconvert==7.16.3 \
    # 因果推断
    dowhy==0.11.1 \
    causalml==0.15.0 \
    econml==0.15.1 \
    # 🔥 GDC API直接访问（替换gdc-client）
    requests==2.32.3 \
    beautifulsoup4==4.12.3 \
    lxml==5.3.0 \
    # TCGA/GEO 下载
    GEOparse==1.4.0 \
    # 空间解卷积 Python 接口
    tangram-py==0.2.0 \
    && pip3 cache purge

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
