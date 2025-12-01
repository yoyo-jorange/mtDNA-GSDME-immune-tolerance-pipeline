############ BASE IMAGE ——— rocker with full R toolchain ############
FROM ghcr.io/rocker-org/r-ver:4.3.1

# Fix Ubuntu source issues
RUN sed -i '/jammy-backports/d' /etc/apt/sources.list && \
    sed -i 's|http://archive.ubuntu.com/ubuntu|http://mirrors.aliyun.com/ubuntu|g' /etc/apt/sources.list && \
    sed -i 's|http://security.ubuntu.com/ubuntu|http://mirrors.aliyun.com/ubuntu|g' /etc/apt/sources.list

# Avoid tzdata prompt issue
ENV DEBIAN_FRONTEND=noninteractive

###########################
# System Dependencies
###########################
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        cmake \
        libxml2-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libjpeg-dev \
        libtiff-dev \
        libcairo2-dev \
        libxt-dev \
        libx11-dev \
        libgsl-dev \
        libgeos-dev \
        libglpk-dev \
        libgdal-dev \
        libhdf5-dev \
        python3 python3-dev python3-pip \
    || { echo "APT INSTALL FAILED — retrying"; apt-get update --fix-missing; apt-get install -y --no-install-recommends libpng-dev; } && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*



###########################
# R — Core Packages
###########################
RUN R -e "install.packages(c('BiocManager'), repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(version='3.18', ask=FALSE)"

# Single-cell & core bioinformatics
RUN R -e "BiocManager::install(c( \
    'SingleCellExperiment', 'scater', 'scran', 'batchelor', \
    'GSVA', 'clusterProfiler', 'enrichplot', 'org.Hs.eg.db' \
), ask=FALSE, update=FALSE)"

# Seurat ecosystem
RUN R -e "install.packages(c( \
    'Seurat', 'SeuratObject', 'SeuratDisk', \
    'ggplot2', 'ggpubr', 'patchwork', \
    'harmony', 'Matrix', 'data.table', \
    'future', 'future.apply', 'reticulate' \
), repos='https://cloud.r-project.org', dependencies=TRUE)"


###########################
# 🔥 R — Spatial Tools (Patch Included)
###########################

# Patch 0 — fix SSL issues for remotes / GitHub
RUN apt-get update && apt-get install -y --no-install-recommends git libgit2-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN git config --global http.sslverify false

# Patch 1 — lock Matrix to safe version (critical for spatialDWLS / Giotto / RCTD)
RUN R -e "install.packages('Matrix', version='1.6-1', repos='https://cloud.r-project.org')"

# Patch 2 — install Eigen toolchain (prevent spatialDWLS compilation crash)
RUN R -e "install.packages(c('RcppEigen','RcppArmadillo'), repos='https://cloud.r-project.org')"

# Patch 3 — ensure BLAS/LAPACK fully available
RUN apt-get update && apt-get install -y --no-install-recommends \
    liblapack-dev libblas-dev pkg-config \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Patch for RCTD — fix BLAS/LAPACK linking issues
RUN echo "FLIBS=-lblas -llapack" >> /usr/local/lib/R/etc/Makeconf

# Install spam and glmnet manually (RCTD critical deps)
RUN R -e "install.packages('spam', repos='https://cloud.r-project.org')"
RUN R -e "install.packages('glmnet', repos='https://cloud.r-project.org')"

# Ensure remotes is installed (critical)
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')"

# Ensure devtools exists (some packages need it)
RUN R -e "install.packages('devtools', repos='https://cloud.r-project.org')"

# Ensure git & libgit2 are available (remotes backend)
RUN apt-get update && apt-get install -y --no-install-recommends \
    liblapack-dev libblas-dev pkg-config \
    wget \
    && apt-get clean && rm -rf /var/lib/apt/lists/*



# RCTD
RUN R -e "remotes::install_github('dmcable/spacexr', dependencies=TRUE, upgrade=FALSE)"


# spatialDWLS (now safe)
# spatialDWLS replacement — DWLS (active repo)
RUN git config --global url."https://ghproxy.com/https://github.com/".insteadOf "https://github.com/"


RUN git config --global url."https://hub.fastgit.xyz/".insteadOf "https://github.com/"


# Giotto (3.1.0)
RUN R -e "remotes::install_github('drieslab/Giotto', dependencies=TRUE, upgrade=FALSE)"


##############################
# Python — Scientific Stack
##############################
RUN pip3 install --upgrade pip

# PyTorch (CPU)
RUN pip3 install torch --index-url https://download.pytorch.org/whl/cpu

# core python scientific stack
RUN pip3 install \
    numpy scipy pandas scikit-learn matplotlib seaborn \
    jupyter jupyterlab ipykernel

# spatial + single-cell python
RUN pip3 install anndata scanpy squidpy

# cell2location + pyro
RUN pip3 install pyro-ppl cell2location

# tangram
RUN pip3 install tangram-sc

# causal inference
RUN pip3 install causalml dowhy econml

##############################
# Snakemake + Workflow Tools
##############################
# 1. 安装 Miniconda（轻量、不会和系统 Python 冲突）
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh && \
    $CONDA_DIR/bin/conda clean -tipsy
ENV PATH="$CONDA_DIR/bin:$PATH"

# 2. 创建专用的 snakemake conda 环境（推荐做法）
RUN conda create -y -n snakemake python=3.11 snakemake=8.20.0 \
    # 常用插件与依赖
    drmaa \
    graphviz \
    pydot \
    jinja2 \
    requests \
    pandas \
    tabulate \
    # 绘图支持（snakemake --report）
    matplotlib \
    seaborn \
    # 集群提交常用
    clusterstatus \
    && conda clean -afy

# 3. 把 conda 环境加入 PATH（以后直接 snakemake 就能用）
ENV PATH="$CONDA_DIR/envs/snakemake/bin:$PATH"

# 4. 可选：安装 singularity/apptainer（很多空间流程用 singularity 拉镜像）
#    如果你服务器没有 fuse-overlayfs 权限，可以只装 apptainer
RUN apt-get update && apt-get install -y --no-install-recommends \
    singularity-container \
    apptainer \
    parallel \
    && rm -rf /var/lib/apt/lists/*

# 5. 让普通用户也能跑 singularity（生产环境常见需求）
RUN echo "user.maxUserNamespaces=10000" >> /etc/singularity/singularity.conf || true

# 6. 启动时默认激活 snakemake 环境（方便直接 docker run -it yourimage）
RUN echo "conda activate snakemake" >> ~/.bashrc
SHELL ["/bin/bash", "-c"]

###############################################
# 🔥 Automatically copy your whole repository
###############################################
COPY . /workspace
WORKDIR /workspace

###############################################
# Default command
###############################################
CMD ["/bin/bash"]
