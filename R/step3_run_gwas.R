#' Run GWAS using MatrixEQTL with Robust Data Handling (Physical Chunking)
#'
#' Automatically handles file slicing, ID alignment, and large-scale regression.
#' Uses physical file splitting to cap memory usage strictly.
#'
#' @param int_mat_file Path to INT normalized matrix (Step 2 output).
#' @param sample_list_file Path to sample ID list (Step 1 output).
#' @param pheno_file Path to Phenotype file (PLINK format or CSV).
#' @param covar_file Path to Covariate file.
#' @param output_file Output path for results.
#' @param phenotype_col Name of the phenotype column to analyze.
#' @param chunk_size Number of SNP rows to process at once (Default 10,000).
#' @param tmp_dir Directory for temporary files.
#' @export
run_nipt_gwas <- function(int_mat_file, sample_list_file,
                          pheno_file, covar_file,
                          output_file, phenotype_col,
                          chunk_size = 10000, # 新增参数：每次处理1万个Bin
                          tmp_dir = tempdir()) {

  requireNamespace("MatrixEQTL")
  requireNamespace("data.table")
  report_mem("Step3: Start")

  # === Init ===
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)

  int_mat_file <- normalizePath(int_mat_file)
  sample_list_file <- normalizePath(sample_list_file)

  # 定义核心临时文件
  # 注意：我们不再生成一个巨大的 clean_matrix，而是生成一个用于取数的源文件
  clean_source_file <- file.path(tmp_dir, "source_matrix_pure.txt")
  bin_ids_file      <- file.path(tmp_dir, "bin_ids.txt")

  message(">> Preparing GWAS Environment...")

  # 1. Load Anchor IDs
  if (!file.exists(sample_list_file)) stop("Sample list file missing.")
  anchor_ids <- readLines(sample_list_file)
  n_expected <- length(anchor_ids)

  # 2. Detect Physical Dimensions
  cmd <- sprintf("head -n 2 %s | tail -n 1 | awk -F'\t' '{print NF}'", shQuote(int_mat_file))
  n_actual_cols <- as.integer(system(cmd, intern = TRUE)) - 1

  message(sprintf("Anchor IDs: %d | Physical Data Columns: %d", n_expected, n_actual_cols))

  # 3. ID Trimming Logic
  final_ids <- anchor_ids
  cut_col_idx <- n_actual_cols + 1

  if (n_actual_cols < n_expected) {
    warning(sprintf("Ghost IDs detected! Trimming IDs from %d to %d.", n_expected, n_actual_cols))
    final_ids <- anchor_ids[1:n_actual_cols]
  } else if (n_actual_cols > n_expected) {
    stop("Data columns exceed Sample IDs.")
  }

  # 4. Generate Pure Source Matrix (No Header, No BinID)
  # 这一步是为了后续切分方便。生成一个纯数字的大文件。
  message("Generating pure data source (system call)...")

  # A. 提取 Bin IDs
  system(sprintf("cut -f 1 %s | tail -n +2 > %s", shQuote(int_mat_file), shQuote(bin_ids_file)))
  bin_ids <- readLines(bin_ids_file)
  n_snps <- length(bin_ids)

  # B. 生成纯数据源 (不含 id 列，不含表头)
  cmd_clean <- sprintf("tail -n +2 %s | cut -f 2-%d > %s",
                       shQuote(int_mat_file), cut_col_idx, shQuote(clean_source_file))

  if (system(cmd_clean) != 0) stop("Failed to create source matrix.")

  report_mem("Step3: Source Prepared")

  # 5. Align Phenotype & Covariates (Ultra-Low Memory Helper)
  .align_meta <- function(infile, ids, target_col = NULL, type_name = "Meta") {
    message(sprintf("Aligning %s (Low Memory Mode)...", type_name))

    # 1. 智能预读：检测是否有 Header
    first_line <- readLines(infile, n = 1)
    has_header <- grepl("[a-zA-Z]", first_line) && !grepl("^V[0-9]+", first_line) # 简单启发式：如果有字母且不是V1,V2

    # 如果用户明确知道有无 Header，最好在 run_nipt_gwas 参数里加一个 header=TRUE/FALSE
    # 这里我们假设通常都有 Header (PLINK .pheno / .cov 通常有 FID IID)

    # 使用 fread 自动处理
    # 如果 target_col 是 "V6" 这种索引型字符串，我们需要小心

    dt <- data.table::fread(infile, header = has_header)

    # 标准化列名：强制前两列为 FID, IID
    old_names <- names(dt)
    if (length(old_names) >= 2) {
      setnames(dt, old_names[1:2], c("FID", "IID"))
    }

    # 确保 IID 是字符型 (防止前几行是纯数字被转为 numeric，导致无法匹配 ID)
    dt[, IID := as.character(IID)]

    # 确定要保留的列
    if (is.null(target_col)) {
      # 保留除 FID, IID 外的所有列 (协变量模式)
      vars_to_keep <- setdiff(names(dt), c("FID", "IID"))
    } else {
      # 表型模式：处理 target_col
      # 如果用户传的是 "V6" 且文件有 Header，这会找不到列
      # 建议用户传真实的列名（如 "GDM"），或者数字索引

      if (target_col %in% names(dt)) {
        vars_to_keep <- target_col
      } else if (grepl("^V[0-9]+", target_col) && has_header) {
        # 如果用户传了 V6 但文件有 Header，尝试转为数字索引
        col_idx <- as.integer(sub("V", "", target_col))
        if (col_idx <= ncol(dt)) {
          vars_to_keep <- names(dt)[col_idx]
          message(sprintf("Mapping %s to column name: %s", target_col, vars_to_keep))
        } else {
          stop("Column index out of bounds.")
        }
      } else {
         # 没 Header 的情况，列名就是 V1, V2...
         vars_to_keep <- target_col
      }
    }

    # 强制将目标列转为 numeric (防止因为 "NA" 字符串导致整列变成 character)
    for (v in vars_to_keep) {
      if (v %in% names(dt)) {
         # 尝试转为 numeric，非数值变为 NA
         dt[[v]] <- as.numeric(as.character(dt[[v]]))
      }
    }

    report_mem(sprintf("Step3: Aligning %s - Fread Done", type_name))

    # 2. 对齐与转置 (使用 data.table 引用语义，避免大对象复制)
    data.table::setkey(dt, IID)

    # 短 ID 映射
    short_ids <- sub("_[0-9]+bp$", "", ids)

    dt <- unique(dt, by = "IID")

    # 提取并排序 (nomatch=NA 自动填充缺失值为 NA)
    # 这一步生成一个新的 data.table，行顺序严格等于 short_ids
    # 我们只提取需要的 columns
    dt_subset <- dt[short_ids, ..vars_to_keep]

    # 显式清理大对象
    rm(dt); gc()

    # 3. 转置并写入
    # 因为 fwrite 不支持直接转置写入，我们需要先转置为 matrix 或 list
    # 对于 13 x 120k 这样形状的矩阵，转置后是 13行 x 120k列

    mat_t <- t(as.matrix(dt_subset))

    # 构造输出 data.table (第一列是 ID)
    out_dt <- data.table::data.table(id = vars_to_keep)
    out_dt <- cbind(out_dt, as.data.table(mat_t))

    # 设置列名 (id + 样本长ID)
    data.table::setnames(out_dt, c("id", ids))

    # 写入文件
    tmp_out <- file.path(tmp_dir, paste0("aligned_", type_name, ".txt"))

    # 关键：MatrixEQTL 对 NA 很敏感，确保写入 "NA" 字符串
    # data.table::fwrite 的 na 参数可以处理
    data.table::fwrite(out_dt, tmp_out, sep = "\t", quote = FALSE, na = "NA")

    rm(dt_subset, mat_t, out_dt); gc()
    report_mem(sprintf("Step3: Aligning %s - Done", type_name))

    return(tmp_out)
  }


  # Generate Aligned Files
  pheno_ready <- .align_meta(pheno_file, final_ids, target_col = phenotype_col, "Phenotype")
  covar_ready <- .align_meta(covar_file, final_ids, target_col = NULL, "Covariates")

  report_mem("Step3: Meta Aligned")

  # 6. Initialize Global MatrixEQTL Objects (Pheno/Covar)
  # 这部分只需加载一次
  message("Loading Phenotype and Covariates into RAM...")

  gene <- MatrixEQTL::SlicedData$new()
  gene$fileDelimiter <- "\t"
  gene$fileOmitCharacters <- "NA"
  gene$fileSkipRows <- 1
  gene$fileSkipColumns <- 1
  gene$LoadFile(pheno_ready)

  cvrt <- MatrixEQTL::SlicedData$new()
  cvrt$fileDelimiter <- "\t"
  cvrt$fileOmitCharacters <- "NA"
  cvrt$fileSkipRows <- 1
  cvrt$fileSkipColumns <- 1
  cvrt$LoadFile(covar_ready)

  report_mem("Step3: Meta Loaded")

  # 7. Physical Chunking Loop (The Memory Solver)

  # 初始化结果文件 (写入 Header)
  cat("SNP\tgene\tbeta\tt-stat\tp-value\tFDR\n", file = output_file)

  n_chunks <- ceiling(n_snps / chunk_size)
  message(sprintf(">> Starting GWAS in %d physical chunks (Size: %d)...", n_chunks, chunk_size))

  # 构造一个带列名的 Header 文件用于合并
  # 我们需要创建一个只有 Header 的文件，用于每次和 chunk 数据合并，骗过 MatrixEQTL
  chunk_header_file <- file.path(tmp_dir, "chunk_header.txt")
  # Header 必须不含 "id" (因为纯数据没有id列), 只有样本名
  writeLines(paste(final_ids, collapse = "\t"), chunk_header_file)

  for (i in 1:n_chunks) {
    # 计算行号范围
    start_row <- (i - 1) * chunk_size + 1
    end_row   <- min(i * chunk_size, n_snps)
    n_rows_chunk <- end_row - start_row + 1

    message(sprintf("--- Processing Chunk %d/%d (Bin %d - %d) ---", i, n_chunks, start_row, end_row))

    # A. 提取纯数据块 (System Call)
    # tail -n +K | head -n N
    raw_chunk <- file.path(tmp_dir, "raw_chunk.txt")
    cmd_extract <- sprintf("tail -n +%d %s | head -n %d > %s",
                           start_row, shQuote(clean_source_file),
                           n_rows_chunk, shQuote(raw_chunk))
    system(cmd_extract)

    # B. 合并 Header + Data
    # MatrixEQTL 需要 Header 才能正确识别列数 (Lazy Mode)
    ready_chunk <- file.path(tmp_dir, "ready_chunk.txt")
    # cat header data > ready
    system(paste("cat", shQuote(chunk_header_file), shQuote(raw_chunk), ">", shQuote(ready_chunk)))

    # C. 加载 SNP Chunk
    snps <- MatrixEQTL::SlicedData$new()
    snps$fileDelimiter <- "\t"
    snps$fileOmitCharacters <- "NA"
    snps$fileSkipRows <- 1      # 跳过 Header
    snps$fileSkipColumns <- 0   # 无行名
    snps$fileSliceSize <- 2000

    snps$LoadFile(ready_chunk)

    # 设置行名
    rownames(snps) <- bin_ids[start_row:end_row]

    # D. 运行引擎
    tmp_res <- file.path(tmp_dir, "chunk_res.txt")

    # 捕获可能的空结果错误
    tryCatch({
      me <- MatrixEQTL::Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = tmp_res,
        pvOutputThreshold = 1.0,  # <--- 修改点 1: 设为 1.0，保留所有P值
                                  # 因为 Bin 总数不多(15w)，全保留也不会很大(约10MB)
                                  # 这对后续 FDR 计算最准确。
        useModel = MatrixEQTL::modelLINEAR,
        errorCovariance = numeric(),
        verbose = FALSE,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
      )

      # E. 追加结果 (跳过 Header)
      if (file.exists(tmp_res)) {
        # 检查是否包含数据
        n_res_lines <- as.integer(strsplit(system(paste("wc -l", tmp_res), intern=TRUE), " ")[[1]][1])
        if (n_res_lines > 1) {
           system(paste("tail -n +2", shQuote(tmp_res), ">>", shQuote(output_file)))
        }
      }
    }, error = function(e) {
      message(sprintf("Warning: Chunk %d failed. Error: %s", i, e$message))
    })

    # F. 清理本次循环内存
    rm(snps, me)
    gc()

    # 每 5 个 Chunk 报一次内存
    if(i %% 5 == 0) report_mem(sprintf("Chunk %d Done", i))
  }

  # =========================================================
  # 8. [新增] Post-hoc FDR Recalculation (全局 FDR 修正)
  # =========================================================
  message(">> Recalculating Global FDR...")

  # 使用 data.table 快速读取合并后的结果
  # 15万行数据瞬间就能读完，内存占用极低
  res_dt <- data.table::fread(output_file, header = TRUE)

  if (nrow(res_dt) > 0) {
    # 核心修正：使用 Benjamini-Hochberg 方法重算 FDR
    # 这里的 n 自动等于总行数 (Total Bins)
    res_dt$FDR <- p.adjust(res_dt$`p-value`, method = "BH")

    # 覆盖原文件
    data.table::fwrite(res_dt, output_file, sep = "\t", quote = FALSE)
    message(sprintf(">> Global FDR updated. Total SNPs: %d", nrow(res_dt)))
  } else {
    warning("Result file is empty.")
  }

  # 9. Cleanup
  file.remove(clean_source_file, bin_ids_file, chunk_header_file)
  if(file.exists(file.path(tmp_dir, "raw_chunk.txt"))) file.remove(file.path(tmp_dir, "raw_chunk.txt"))
  if(file.exists(file.path(tmp_dir, "ready_chunk.txt"))) file.remove(file.path(tmp_dir, "ready_chunk.txt"))
  if(file.exists(file.path(tmp_dir, "chunk_res.txt"))) file.remove(file.path(tmp_dir, "chunk_res.txt"))
  if (file.exists(pheno_ready)) file.remove(pheno_ready)
  if (file.exists(covar_ready)) file.remove(covar_ready)

  message("GWAS Done.")
}