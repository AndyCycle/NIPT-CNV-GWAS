#' Run GWAS using MatrixEQTL with Robust Data Handling
#'
#' Automatically handles file slicing, ID alignment, and large-scale regression.
#'
#' @param int_mat_file Path to INT normalized matrix (Step 2 output).
#' @param sample_list_file Path to sample ID list (Step 1 output).
#' @param pheno_file Path to Phenotype file (PLINK format or CSV).
#' @param covar_file Path to Covariate file.
#' @param output_file Output path for results.
#' @param phenotype_col Name of the phenotype column to analyze.
#' @param tmp_dir Directory for temporary files.
#' @export
run_nipt_gwas <- function(int_mat_file, sample_list_file,
                          pheno_file, covar_file,
                          output_file, phenotype_col,
                          tmp_dir = tempdir()) {

  requireNamespace("MatrixEQTL")

  # === Fix 1: 确保临时目录存在 ===
  if (!dir.exists(tmp_dir)) {
    message(sprintf("Creating temporary directory: %s", tmp_dir))
    dir.create(tmp_dir, recursive = TRUE)
  }

  # 使用绝对路径以避免系统命令路径解析错误
  int_mat_file <- normalizePath(int_mat_file)
  sample_list_file <- normalizePath(sample_list_file)
  # 临时文件路径
  clean_data_file <- file.path(tmp_dir, "clean_matrix_for_eqtl.txt")
  bin_ids_file    <- file.path(tmp_dir, "bin_ids.txt")

  message(">> Preparing GWAS Environment...")

  # 1. Load Anchor IDs
  if (!file.exists(sample_list_file)) stop("Sample list file missing.")
  anchor_ids <- readLines(sample_list_file)
  n_expected <- length(anchor_ids)

  # 2. Detect Physical Dimensions
  # 检测纯数据列数（不含ID列）
  cmd <- sprintf("head -n 2 %s | tail -n 1 | awk -F'\t' '{print NF}'", shQuote(int_mat_file))
  n_actual_cols <- as.integer(system(cmd, intern = TRUE)) - 1

  message(sprintf("Anchor IDs: %d | Physical Data Columns: %d", n_expected, n_actual_cols))

  # 3. ID Trimming Logic (处理 Ghost IDs)
  final_ids <- anchor_ids
  cut_col_idx <- n_actual_cols + 1 # +1 是因为文件里还有第1列 BinID

  if (n_actual_cols < n_expected) {
    warning(sprintf("Ghost IDs detected! Trimming IDs from %d to %d.", n_expected, n_actual_cols))
    final_ids <- anchor_ids[1:n_actual_cols]
  } else if (n_actual_cols > n_expected) {
    stop("Data columns exceed Sample IDs. Critical structural error.")
  }

  # 4. Create Clean Data Slice (No Header, No BinID)
  message("Slicing CNV matrix (system call)...")
  # Fix: 确保使用 shQuote 处理路径，防止空格报错
  cmd_clean <- sprintf("tail -n +2 %s | cut -f 2-%d > %s",
                       shQuote(int_mat_file), cut_col_idx, shQuote(clean_data_file))
  if (system(cmd_clean) != 0) stop("Failed to slice data matrix. Check disk space or permissions.")

  # 5. Extract Bin IDs
  cmd_bins <- sprintf("cut -f 1 %s | tail -n +2 > %s", shQuote(int_mat_file), shQuote(bin_ids_file))
  system(cmd_bins)
  bin_ids <- readLines(bin_ids_file)

  # 6. Align Phenotype & Covariates (Internal Helper)
  .align_meta <- function(infile, ids, target_col = NULL, type_name = "Meta") {
    message(sprintf("Aligning %s...", type_name))
    df <- data.table::fread(infile, header = FALSE, data.table = FALSE)

    # PLINK Header Guessing
    colnames(df)[1:2] <- c("FID", "IID")
    if (ncol(df) >= 6 && is.null(target_col)) {
       # 如果没指定列，默认不做处理，后面会报错
    }
    rownames(df) <- df$IID

    # Short ID mapping
    short_ids <- sub("_[0-9]+bp$", "", ids)

    # Determine columns to keep
    vars_to_keep <- if(is.null(target_col)) setdiff(colnames(df), c("FID", "IID")) else target_col

    # Build Output Matrix (Row = Var, Col = Sample)
    mat <- matrix(NA, nrow = length(vars_to_keep), ncol = length(ids))
    colnames(mat) <- ids
    rownames(mat) <- vars_to_keep

    # Intersection
    common <- intersect(short_ids, rownames(df))
    if (length(common) == 0) stop(sprintf("No overlapping samples found in %s file.", type_name))

    # Fill Data
    map_idx <- match(common, short_ids)
    long_ids_common <- ids[map_idx]

    mat[, long_ids_common] <- t(df[common, vars_to_keep, drop=FALSE])

    # Write to temp file
    tmp_out <- file.path(tmp_dir, paste0("aligned_", type_name, ".txt"))
    out_df <- data.frame(id = rownames(mat), mat, check.names = FALSE)

    # Force NA handling
    out_df[is.na(out_df)] <- "NA"

    data.table::fwrite(out_df, tmp_out, sep = "\t", quote = FALSE)
    return(tmp_out)
  }

  # Generate Aligned Files
  pheno_ready <- .align_meta(pheno_file, final_ids, target_col = phenotype_col, "Phenotype")
  covar_ready <- .align_meta(covar_file, final_ids, target_col = NULL, "Covariates")

  # 7. MatrixEQTL Setup
  message("Initializing MatrixEQTL SlicedData...")

  # SNP Data
  snps <- MatrixEQTL::SlicedData$new()
  snps$fileDelimiter <- "\t"
  snps$fileOmitCharacters <- "NA"
  snps$fileSkipRows <- 0    # 纯数据，无表头
  snps$fileSkipColumns <- 0 # 纯数据，无ID列
  snps$fileSliceSize <- 2000
  snps$LoadFile(clean_data_file)
  rownames(snps) <- bin_ids # 恢复行名

  # Gene Data (Phenotype)
  gene <- MatrixEQTL::SlicedData$new()
  gene$fileDelimiter <- "\t"
  gene$fileOmitCharacters <- "NA"
  gene$fileSkipRows <- 1
  gene$fileSkipColumns <- 1
  gene$LoadFile(pheno_ready)

  # Covariates
  cvrt <- MatrixEQTL::SlicedData$new()
  cvrt$fileDelimiter <- "\t"
  cvrt$fileOmitCharacters <- "NA"
  cvrt$fileSkipRows <- 1
  cvrt$fileSkipColumns <- 1
  cvrt$LoadFile(covar_ready)

  # === Final Alignment Check (MatrixEQTL Best Practice) ===
  # MatrixEQTL 引擎假设列顺序完全一致。
  # 这里的检查确保我们之前的对齐逻辑没有出错。
  if (snps$nCols() != gene$nCols() || snps$nCols() != cvrt$nCols()) {
    stop(sprintf("Dimension Mismatch! SNP(%d) vs Pheno(%d) vs Covar(%d)",
                 snps$nCols(), gene$nCols(), cvrt$nCols()))
  }

  # 8. Run Engine
  message("Running Linear Regression Engine...")
  me <- MatrixEQTL::Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file,
    pvOutputThreshold = 1e-2, # 输出 P < 0.01 的结果，节省空间
    useModel = MatrixEQTL::modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )

  # 9. Cleanup
  # 清理临时文件以节省空间，特别是那个巨大的 clean_data_file
  if (file.exists(clean_data_file)) file.remove(clean_data_file)
  if (file.exists(bin_ids_file)) file.remove(bin_ids_file)
  if (file.exists(pheno_ready)) file.remove(pheno_ready)
  if (file.exists(covar_ready)) file.remove(covar_ready)

  message("GWAS Done. Temporary files cleaned.")
}