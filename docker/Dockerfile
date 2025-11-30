############ BASE IMAGE ——— rocker with full R toolchain ############
FROM ghcr.io/rocker-org/r-ver:4.3.1

# Avoid tzdata prompt issue
ENV DEBIAN_FRONTEND=noninteractive

###########################
# System Dependencies
###########################
RUN apt-get update && apt-get install -y --no-install-recommends \
    ## basic compilers
    build-essential \
    gfortran \
    cmake \
    ## libraries for R / Bioconductor / Seurat / GSVA
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
    ## python basics
    python3 python3-dev python3-pip \
    ## cleanup
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


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
RUN apt-get update && apt-get install -y --no-install-recommends git libgit2-dev && apt-get clean && rm -rf /var/lib/apt/lists/*


# RCTD
RUN R -e "remotes::install_github('dmcable/RCTD@e2e898a04dbbfbbb0f1f28fef7840ccf44a8a893', dependencies=TRUE, upgrade=FALSE)"

# spatialDWLS (now safe)
RUN R -e "remotes::install_github('dmcable/spatialDWLS', dependencies=TRUE, upgrade=FALSE)"

# Giotto (3.1.0)
RUN R -e "remotes::install_github('lzheimer/Giotto', ref='3.1.0', dependencies=TRUE, upgrade=FALSE)"


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


###############################################
# 🔥 Automatically copy your whole repository
###############################################
COPY . /workspace
WORKDIR /workspace

###############################################
# Default command
###############################################
CMD ["/bin/bash"]
