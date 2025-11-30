#!/usr/bin/env Rscript
library(Seurat)
library(future)
library(harmony)
plan("multisession", workers = 4)

# 加载配置
config <- yaml::read_yaml("config/paths.yaml")
raw_dir <- file.path(config$data_raw, "scrna")
proc_dir <- file.path(config$data_processed, "scrna")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

# GEO单细胞数据预处理函数
process_geo_scrna <- function(geo_id) {
  cat("处理GEO:", geo_id, "\n")
  
  # 下载GEO数据（示例：GSEXXXXX）
  if (!dir.exists(file.path(raw_dir, geo_id))) {
    SeuratData::InstallData(geo_id)  # 或使用GEOquery下载
  }
  
  # 读取并标准化
  sobj <- readRDS(file.path(raw_dir, geo_id, "processed.rds"))  # 假设已下载
  sobj <- SCTransform(sobj, verbose = FALSE)
  sobj <- RunPCA(sobj, verbose = FALSE)
  sobj <- RunUMAP(sobj, dims = 1:30, verbose = FALSE)
  
  # 细胞类型注释（使用你的18种细胞类型）
  celltypes <- c("Tumor_cells", "Macrophages_GSDME_pos", "Treg", "MDSC", 
                 "CD8_T_exhausted", "CD8_T_active")  # 你的完整列表
  
  # 保存参考
  saveRDS(sobj, file.path(proc_dir, paste0(geo_id, "_ref.rds")))
  cat("✅ 保存:", geo_id, "\n")
}

# 处理所有单细胞数据集
geo_ids <- c("GSE211732", "GSE236263")  # 你的cfDNA+scRNA队列
sapply(geo_ids, process_geo_scrna)

cat("🎉 单细胞参考数据预处理完成\n")
