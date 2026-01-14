#!/usr/bin/env Rscript

# ==============================================================================
# 脚本名称：step4_visualize.R
# 功能：解析 GWAS 结果，绘制曼哈顿图和 QQ 图
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qqman)
  library(stringr)
  library(dplyr)
})

# ==================== 配置 ====================
# 输入文件 (Step 3 的输出)
gwas_file <- "./test_run/output_120k_smart/test_cohort_smart.gwas_results.txt"
# 输出前缀
out_prefix <- "./test_run/output_120k_smart/test_cohort_smart.plot"

# 显著性阈值 (画线用)
suggestive_line <- -log10(1e-5)  # 蓝色虚线
genome_wide_line <- -log10(5e-8) # 红色实线 (GWAS标准)
# 或者使用 Bonferroni 矫正阈值: -log10(0.05 / 总Bin数)

# ==================== 1. 读取数据 ====================
message(">> Loading GWAS results...")
gwas <- fread(gwas_file)

if (nrow(gwas) == 0) stop("Result file is empty!")

# ==================== 2. 解析坐标 ====================
message(">> Parsing genomic coordinates...")

# 我们的 SNP 列格式为: chr1_10000_30000
# 需要拆分为 CHR 和 BP (Start)
# 使用 stringr 或 tidyr 拆分
# 注意：有些 chr 可能是 "chrX", "chrY", 需要特殊处理

plot_data <- gwas %>%
  # 拆分 SNP 字符串
  mutate(
    temp_id = str_remove(SNP, "chr"), # 去掉 'chr' 前缀: 1_10000_30000
  ) %>%
  tidyr::separate(temp_id, into = c("CHR_STR", "BP", "END"), sep = "_", convert = TRUE) %>%
  mutate(
    # 将 X, Y, M 转为数字，以便 qqman 排序
    CHR = case_when(
      CHR_STR == "X" ~ 23,
      CHR_STR == "Y" ~ 24,
      CHR_STR == "M" | CHR_STR == "MT" ~ 25,
      TRUE ~ as.numeric(CHR_STR)
    ),
    P = `p-value`
  ) %>%
  filter(!is.na(CHR)) %>% # 过滤掉无法识别的染色体
  select(SNP, CHR, BP, P, beta)

message(sprintf(">> Ready to plot %d variants.", nrow(plot_data)))

# ==================== 3. 绘制 QQ 图 ====================
# 检查是否存在系统性膨胀 (Inflation)
message(">> Generating QQ Plot...")

png(paste0(out_prefix, "_QQ.png"), width = 1200, height = 1200, res = 150)
qq(plot_data$P, main = "Q-Q Plot of GDM CNV-GWAS")
dev.off()

# 计算 Lambda GC (Inflation Factor)
# 理想值应接近 1.0。如果 > 1.1 说明可能有未校正的群体分层或批次效应
chisq <- qchisq(1 - plot_data$P, 1)
lambda_gc <- median(chisq) / qchisq(0.5, 1)
message(sprintf(">> Lambda GC (Inflation Factor): %.4f", lambda_gc))

# ==================== 4. 绘制曼哈顿图 ====================
message(">> Generating Manhattan Plot...")

png(paste0(out_prefix, "_Manhattan.png"), width = 2000, height = 1000, res = 150)

# 颜色设置：染色体交替颜色
col_set <- c("#4197d8", "#f8c120", "#413496", "#495226", "#d60c00")

manhattan(plot_data,
          chr = "CHR",
          bp = "BP",
          p = "P",
          snp = "SNP",
          col = col_set,
          suggestiveline = suggestive_line,
          genomewideline = genome_wide_line,
          annotatePval = 0.01, # 如果有极显著的点，标注出来
          main = sprintf("Manhattan Plot (Lambda GC = %.3f)", lambda_gc))

dev.off()

message(">> Plots saved.")