import os
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd

PROCESSED_SPATIAL_DIR = Path("data_processed/spatial")
OUT_DIR = Path("data_processed/with_scores")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 线粒体基因前缀（人类通常使用 "MT-"）
MT_PREFIX = "MT-"

# 这里需要学生根据实际情况准备一个 OXPHOS 基因列表
OXPHOS_GENES_TXT = Path("config/oxphos_genes.txt")

def load_oxphos_genes():
    if not OXPHOS_GENES_TXT.exists():
        print(f"[警告] 找不到 {OXPHOS_GENES_TXT}，请填写 OXPHOS 基因列表。")
        return []
    with open(OXPHOS_GENES_TXT, encoding="utf-8") as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes

def compute_oxphos_score(adata, oxphos_genes):
    # 简单版：每个 spot/细胞中 OXPHOS 基因的平均表达
    gene_idx = [g for g in oxphos_genes if g in adata.var_names]
    if not gene_idx:
        print("[警告] 本数据中找不到 OXPHOS 基因，得分将为 0。")
        return np.zeros(adata.n_obs)
    sub = adata[:, gene_idx].X
    if hasattr(sub, "toarray"):
        sub = sub.toarray()
    return np.asarray(sub).mean(axis=1)

def process_one_file(path, oxphos_genes):
    print(f"处理: {path}")
    adata = ad.read_h5ad(path)

    # 1) 线粒体基因表达
    mt_genes = [g for g in adata.var_names if g.startswith(MT_PREFIX)]
    if mt_genes:
        mt_expr = adata[:, mt_genes].X
        if hasattr(mt_expr, "toarray"):
            mt_expr = mt_expr.toarray()
        mt_sum = np.asarray(mt_expr).sum(axis=1)
    else:
        print("[警告] 未检测到以 'MT-' 开头的基因，将 mt 表达设为 0")
        mt_sum = np.zeros(adata.n_obs)

    # 2) OXPHOS score
    ox_score = compute_oxphos_score(adata, oxphos_genes)

    # 3) mtLeakScore = log2(mt_sum + 1) * ( - OXPHOS score )
    mt_log = np.log2(mt_sum + 1.0)
    mt_leak = mt_log * (-ox_score)

    adata.obs["mt_sum"] = mt_sum
    adata.obs["oxphos_score"] = ox_score
    adata.obs["mtLeakScore"] = mt_leak

    out_path = OUT_DIR / path.name
    adata.write_h5ad(out_path)
    print(f"已保存: {out_path}")

def main():
    oxphos_genes = load_oxphos_genes()
    if not PROCESSED_SPATIAL_DIR.exists():
        print(f"[错误] 找不到目录: {PROCESSED_SPATIAL_DIR}，请先运行预处理脚本。")
        return

    for path in PROCESSED_SPATIAL_DIR.glob("*.h5ad"):
        process_one_file(path, oxphos_genes)

if __name__ == "__main__":
    main()
