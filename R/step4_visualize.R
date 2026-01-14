#' Visualize GWAS Results (Manhattan & QQ)
#'
#' @param gwas_file Path to Step 3 output file.
#' @param output_prefix Prefix for output plots.
#' @param title Plot title.
#' @export
visualize_gwas <- function(gwas_file, output_prefix, title = "CNV-GWAS") {

  requireNamespace("data.table")
  requireNamespace("qqman")
  requireNamespace("dplyr")
  requireNamespace("tidyr")
  requireNamespace("stringr")

  message(">> Step 4: Visualizing Results...")

  # 1. 读取数据
  gwas <- data.table::fread(gwas_file)
  if (nrow(gwas) == 0) {
    warning("GWAS result file is empty. Skipping plots.")
    return(NULL)
  }

  # 2. 解析坐标 (兼容 Step 3 输出的 SNP 列格式)
  # 假设 SNP 列是: chr1_10000_30000 或 1_10000_30000
  plot_data <- gwas %>%
    dplyr::mutate(
      temp_id = stringr::str_remove(SNP, "^chr"), # 去掉可能的 chr 前缀
    ) %>%
    tidyr::separate(temp_id, into = c("CHR_STR", "BP", "END"), sep = "_", convert = TRUE, extra = "drop") %>%
    dplyr::mutate(
      CHR = dplyr::case_when(
        CHR_STR == "X" ~ 23,
        CHR_STR == "Y" ~ 24,
        CHR_STR %in% c("M", "MT") ~ 25,
        TRUE ~ as.numeric(CHR_STR)
      ),
      P = `p-value`
    ) %>%
    dplyr::filter(!is.na(CHR), !is.na(P)) %>%
    dplyr::select(SNP, CHR, BP, P)

  if(nrow(plot_data) == 0) {
      warning("No valid variants for plotting after parsing.")
      return(NULL)
  }

  # 3. 绘制 QQ 图
  png(paste0(output_prefix, "_QQ.png"), width = 1200, height = 1200, res = 150)
  qqman::qq(plot_data$P, main = paste0("Q-Q Plot: ", title))
  dev.off()

  # 计算 Lambda GC
  chisq <- qchisq(1 - plot_data$P, 1)
  lambda_gc <- median(chisq) / qchisq(0.5, 1)
  message(sprintf("   Lambda GC: %.4f", lambda_gc))

  # 4. 绘制曼哈顿图
  png(paste0(output_prefix, "_Manhattan.png"), width = 2000, height = 1000, res = 150)

  # 显著线设置
  suggestive <- -log10(1e-5)
  genomewide <- -log10(5e-8)

  col_set <- c("#4197d8", "#f8c120", "#413496", "#495226", "#d60c00")

  tryCatch({
      qqman::manhattan(plot_data,
               chr = "CHR", bp = "BP", p = "P", snp = "SNP",
               col = col_set,
               suggestiveline = suggestive,
               genomewideline = genomewide,
               main = paste0(title, " (Lambda GC = ", round(lambda_gc, 3), ")"))
  }, error = function(e) {
      message("Manhattan plot failed (possibly due to data limits): ", e$message)
  })

  dev.off()
  message(">> Step 4 Done.")
}