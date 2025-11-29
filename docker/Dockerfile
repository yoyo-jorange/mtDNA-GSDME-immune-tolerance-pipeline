FROM rocker/tidyverse:4.4.0

# 基础工具
RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-dev \
    git wget curl nano \
    && rm -rf /var/lib/apt/lists/*

# 设置 Python 为 python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1

# 安装 Python 包
RUN pip install --no-cache-dir \
    numpy pandas scipy scikit-learn \
    matplotlib seaborn \
    scanpy anndata \
    statsmodels lifelines \
    jupyter jupyterlab \
    # 因果推断
    causalml dowhy \
    # 可按需添加更多包
    ;

# 安装 R 包（Seurat 等）
RUN R -e "install.packages(c('BiocManager','remotes'), repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('Seurat','SingleCellExperiment','SingleR','edgeR','limma','GSVA'), ask=FALSE, update=FALSE)" && \
    R -e "install.packages(c('reticulate','patchwork','data.table'), repos='https://cloud.r-project.org')"

# 空间解卷积相关 R 包（这里只放骨架，学生可按需补充具体包名）
RUN R -e "install.packages(c('sp','spdep'), repos='https://cloud.r-project.org')"

WORKDIR /workspace
COPY . /workspace

CMD ["/bin/bash"]
