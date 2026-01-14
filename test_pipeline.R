#!/usr/bin/env Rscript

# ==============================================================================
# 功能：使用 Smart Wrapper 进行全流程测试 (支持自动断点续传)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(HMMcopy)
  library(MatrixEQTL)
  library(qqman)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(methods)
})

# 加载你的函数库 (确保 R 文件夹下包含 step1-4, utils, wrapper 等所有文件)
script_dir <- "./R"
sapply(list.files(script_dir, pattern = "\\.R$", full.names = TRUE), source)

# ==============================================================================
# 1. 核心路径配置 (120k Run)
# ==============================================================================

# 输出目录
work_dir    <- "./output"
file_prefix <- "Cohort"

# 输入目录 (共同父目录，涵盖 Baoan 和 Longgang)
# 脚本会自动递归搜索该目录下所有的 .wig 文件
raw_wig_dir <- "projects/hmmcopy/samples/readcounts_wigs"

# 引用文件
gc_wig  <- "Homo_sapiens_assembly38_20000bp.gc.wig"
map_wig <- "k100.Umap.MultiTrackMappability_20000bp.wig"

# 表型数据 (假设列名为 FID, IID, GDM_Status)
pheno_file <- "pheno_with_header.table"
# 协变量 (PCs, Age, Batch 等)
covar_file <- "XXXX_NIPT_plink.PC1-10.68632.maternalAge.bmi.OGTTwk.xCov.txt"

# ==============================================================================
# 2. 参数配置
# ==============================================================================

pipeline_args <- list(
  # Step 1: 矩阵构建
  step1 = list(
    input_dir    = raw_wig_dir,   # 这里传入父目录
    gc_file      = gc_wig,
    map_file     = map_wig,
    batch_size   = 2000,          # 120k样本建议加大Batch，减少碎片文件
    n_cores      = 30,            # 建议申请 30-40 核心
    qc_mapd_th   = 0.35,
    qc_median_th = 0.2,
    file_pattern = "_20000bp\\.readcounts\\.wig$"     # 正则解释：测试文件包含 _20000bp 且以 .readcounts.wig 结尾
  ),

  # Step 2: INT 标准化
  step2 = list(k = 0.375),

  # Step 3: GWAS (1/2 编码无影响，线性模型自动适应)
  step3 = list(
    pheno_file    = pheno_file,
    covar_file    = covar_file,
    phenotype_col = "GDM_Status", # 替换为你表型文件中实际的列名 (1=Con, 2=Case)
    chunk_size    = 20000         # 每次处理2万个Bin，根据内存调整
  ),

  # Step 4: 可视化
  step4 = list(
    title = "CNV Risk GWAS (N=XXXX, 20kb bins)"
  )
)

# ==============================================================================
# 3. 启动全流程
# ==============================================================================

message(">> [Cohort] Starting Pipeline...")
message(sprintf("   Input Root: %s", raw_wig_dir))

if (!dir.exists(work_dir)) {
    message(sprintf("Creating output directory: %s", work_dir))
    dir.create(work_dir, recursive = TRUE)
}

run_pipeline_smart(
  work_dir        = work_dir,
  prefix          = file_prefix,
  input_args      = pipeline_args,
  force_overwrite = FALSE # 断点续传开启：如果某一步中断，下次直接从该步开始
)

message(">> [120k Cohort] Analysis Complete.")