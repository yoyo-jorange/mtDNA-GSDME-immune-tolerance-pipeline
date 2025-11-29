import yaml
import pandas as pd
from pathlib import Path
import gdc_client
import GEOparse

def download_visium_hd_12847():
    """下载10x 2025公共Visium HD数据集"""
    # 你的具体链接：https://www.10xgenomics.com/datasets
    urls = [
        "https://cf.10xgenomics.com/samples/cell-exp/3.0.2/Visium_HD_Human_Breast_Cancer/Visium_HD_Human_Breast_Cancer_raw_feature_bc_matrix.h5",
        # ... 继续你的12,847个样本链接
    ]
    return urls

def main():
    with open("config/paths.yaml") as f:
        config = yaml.safe_load(f)
    
    raw_dir = Path(config["data_raw"])
    raw_dir.mkdir(exist_ok=True)
    
    # Visium HD (核心)
    if config["download"]["visium_hd"]:
        urls = download_visium_hd_12847()
        for i, url in enumerate(urls[:5]):  # 先下载5个测试
            outpath = raw_dir / f"visium_hd/sample_{i+1}.h5"
            # 下载逻辑...
            print(f"下载: {url} -> {outpath}")
