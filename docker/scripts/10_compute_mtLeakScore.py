#!/usr/bin/env python3
"""
mtDNA泄漏评分计算脚本
公式: mtLeakScore = log2(平均线粒体基因UMIs in cytoplasm + 1) × (OXPHOS signature score ↓)
支持: Visium HD空间数据、单细胞、TCGA bulk
"""

import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import yaml
import matplotlib.pyplot as plt
import seaborn as sns
from gseapy import gseapy
import warnings
warnings.filterwarnings('ignore')

def load_oxphos_geneset():
    """加载OXPHOS基因集（MSigDB HALLMARK_OXIDATIVE_PHOSPHORYLATION）"""
    oxphos_genes = [
        'MT-ATP6', 'MT-ATP8', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'MT-CYB',
        'MT-ND1', 'MT-ND2', 'MT-ND3', 'MT-ND4', 'MT-ND4L', 'MT-ND5', 'MT-ND6',
        'NDUFA1', 'NDUFA2', 'NDUFA3', 'NDUFA4', 'NDUFA5', 'NDUFA6', 'NDUFA7',
        'NDUFA8', 'NDUFA9', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13',
        'NDUFB1', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7',
        'NDUFB8', 'NDUFB9', 'NDUFB10', 'NDUFB11', 'NDUFC1', 'NDUFC2',
        'NDUFV1', 'NDUFV2', 'NDUFV3', 'COX4I1', 'COX4I2', 'COX5A', 'COX5B',
        'COX6A1', 'COX6B1', 'COX7A1', 'COX7A2', 'COX7B', 'ATP5A1', 'ATP5B',
        'ATP5C1', 'ATP5D', 'ATP5F1', 'ATP5G1', 'ATP5G2', 'ATP5G3', 'ATP5H',
        'ATP5I', 'ATP5J', 'ATP5J2', 'ATP5O', 'UQCRB', 'UQCRC1', 'UQCRC2',
        'UQCRFS1', 'UQCRQ', 'CYC1', 'CYCS'
    ]
    return [g.upper() for g in oxphos_genes]

def compute_oxphos_score(adata, oxphos_genes):
    """计算OXPHOS signature score (GSVA风格)"""
    available_genes = [g for g in oxphos_genes if g in adata.var_names]
    if len(available_genes) < 5:
        print("⚠️  OXPHOS基因过少，使用平均表达")
        return np.zeros(adata.n_obs)
    
    # 提取OXPHOS基因表达
    oxphos_expr = adata[:, available_genes].X
    if hasattr(oxphos_expr, 'toarray'):
        oxphos_expr = oxphos_expr.toarray()
    
    # GSVA-like score: 平均标准化表达
    oxphos_mean = np.mean(oxphos_expr, axis=1)
    oxphos_z = (oxphos_mean - np.nanmean(oxphos_mean)) / np.nanstd(oxphos_mean)
    return oxphos_z

def compute_mtleak_score(adata):
    """计算mtLeakScore"""
    print("🔬 计算mtLeakScore...")
    
    # 1. 识别线粒体基因 (MT-前缀)
    mt_genes = [g for g in adata.var_names if g.startswith(('MT-', 'mt-'))]
    print(f"📊 发现 {len(mt_genes)} 个线粒体基因")
    
    if len(mt_genes) == 0:
        print("⚠️  未找到MT-基因，mt表达设为0")
        mt_expr = np.zeros((adata.n_obs, 1))
    else:
        mt_expr = adata[:, mt_genes].X
        if hasattr(mt_expr, 'toarray'):
            mt_expr = mt_expr.toarray()
        mt_sum = np.sum(mt_expr, axis=1)  # 每个spot的总mt UMI/counts
    
    # 2. log2(mt + 1) 转换
    mt_log = np.log2(mt_sum + 1)
    
    # 3. OXPHOS score (负相关)
    oxphos_genes = load_oxphos_geneset()
    oxphos_score = compute_oxphos_score(adata, oxphos_genes)
    
    # 4. mtLeakScore = log2(mt+1) × (-OXPHOS)
    mt_leak
