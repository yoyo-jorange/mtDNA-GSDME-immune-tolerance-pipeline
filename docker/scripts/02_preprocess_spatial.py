#!/usr/bin/env python3
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
from pathlib import Path
import yaml
from anndata import AnnData

def load_visium_hd(path):
    """加载10x Visium HD数据（.h5格式）"""
    adata = sc.read_10x_h5(str(path))
    sq.pl.extract_table(adata, library_id="visium_hd")
    return adata

def qc_spatial(adata, min_genes=200, max_genes=6000, min_counts=500):
    """空间数据QC"""
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pl.highest_expr_genes(adata, n_top=20, save="qc_top_genes.png")
    
    # 过滤
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=3)
    return adata

def main():
    with open("config/paths.yaml") as f:
        config = yaml.safe_load(f)
    
    raw_dir = Path(config["data_raw"]) / "visium_hd"
    proc_dir = Path(config["data_processed"]) / "spatial"
    proc_dir.mkdir(parents=True, exist_ok=True)
    
    processed_files = []
    for h5_file in raw_dir.glob("*.h5"):
        print(f"处理: {h5_file.name}")
        try:
            adata = load_visium_hd(h5_file)
            adata = qc_spatial(adata)
            
            # 标准化基因名（HGNC）
            adata.var_names = adata.var_names.str.upper()
            
            out_path = proc_dir / f"{h5_file.stem}_processed.h5ad"
            adata.write_h5ad(out_path)
            processed_files.append(str(out_path))
            print(f"✅ 保存: {out_path}")
        except Exception as e:
            print(f"❌ 失败 {h5_file.name}: {e}")
    
    # 保存处理清单
    pd.DataFrame({"processed_files": processed_files}).to_csv(
        "results/tables/spatial_processed_list.csv", index=False
    )
    print("📋 处理完成清单已保存")

if __name__ == "__main__":
    main()
