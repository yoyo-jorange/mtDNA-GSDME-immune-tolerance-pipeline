import os
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd

INPUT_DIR = Path("data_processed/with_scores")
OUT_DIR = Path("data_processed/with_scores")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 目标基因名称（需与表达矩阵中的 var_names 一致）
GENES = {
    "GSDME": "GSDME",
    "CASP3": "CASP3",
    "IL1B": "IL1B",
    "IL18": "IL18",
    "GSDMD": "GSDMD",
    "CASP1": "CASP1",
}

def get_gene_expr(adata, gene_name):
    if gene_name not in adata.var_names:
        print(f"[警告] 基因 {gene_name} 不在 var_names 中，视为 0")
        return np.zeros(adata.n_obs)
    x = adata[:, gene_name].X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x).flatten()

def zscore(x):
    m = np.nanmean(x)
    s = np.nanstd(x)
    if s == 0:
        return np.zeros_like(x)
    return (x - m) / s

def process_one_file(path):
    print(f"处理: {path}")
    adata = ad.read_h5ad(path)

    expr = {}
    for key, g in GENES.items():
        expr[key] = get_gene_expr(adata, g)

    # Z-score 标准化
    for k in expr:
        expr[k] = zscore(expr[k])

    numerator = expr["GSDME"] * expr["CASP3"] * (expr["IL1B"] + expr["IL18"])
    denominator = expr["GSDMD"] + expr["CASP1"]
    denominator = np.where(denominator == 0, 1e-6, denominator)

    score_raw = numerator / denominator

    # 再整体做一次 Z-score，方便跨样本比较
    score = zscore(score_raw)

    adata.obs["nonCanonicalPyroScore_raw"] = score_raw
    adata.obs["nonCanonicalPyroScore"] = score

    out_path = OUT_DIR / path.name
    adata.write_h5ad(out_path)
    print(f"已保存: {out_path}")

def main():
    if not INPUT_DIR.exists():
        print(f"[错误] 找不到目录: {INPUT_DIR}")
        return

    for path in INPUT_DIR.glob("*.h5ad"):
        process_one_file(path)

if __name__ == "__main__":
    main()
