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
                          n_cores = 1, # <--- 新增参数：并行核数
                          tmp_dir = tempdir()) {

  requireNamespace("MatrixEQTL")
  requireNamespace("data.table")
  requireNamespace("parallel") # 引入并行包
  report_mem("Step3: Start")

  # === Init ===
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)

  int_mat_file <- normalizePath(int_mat_file)
  sample_list_file <- normalizePath(sample_list_file)

  # 定义核心临时文件
  # 注意：我们不再生成一个巨大的 clean_matrix，而是生成一个用于取数的源文件
  clean_source_file <- file.path(tmp_dir, "source_matrix_pure.txt")
  bin_ids_file      <- file.path(tmp_dir, "bin_ids.txt")

  message(sprintf(">> Preparing GWAS Environment (Cores: %d)...", n_cores))

  # 1. Load Anchor IDs
  if (!file.exists(sample_list_file)) stop("Sample list file missing.")
  anchor_ids <- readLines(sample_list_file)
  anchor_ids <- anchor_ids[anchor_ids != ""] # 去除空行
  anchor_ids <- unique(anchor_ids)           # 再次强制去重
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

  # =========================================================
  # [新增模块] 统计有效样本量 (Calculate Effective Sample Size)
  # =========================================================
  message(">> Calculating effective sample sizes...")

  # 1. 快速读取对齐后的数据 (包含 FID, IID, Data...)
  p_dt <- data.table::fread(pheno_ready)
  c_dt <- data.table::fread(covar_ready)

  # 2. 找到有效样本交集 (Pheno非NA 且 Covar完整非NA)
  # 注意：MatrixEQTL 会自动剔除任何含有 NA 的列(样本)

  # 获取表型值 (假设第3列开始是数据，前两列是FID IID)
  # 我们的 .align_meta 生成的文件格式是：ID(列名) ...
  # MatrixEQTL转置后：行是变量，列是样本。但 .align_meta 生成的是 MatrixEQTL 格式：
  # 第一行是 Header (ID)，第一列是 RowName。
  # 等等，我们需要确认 .align_meta 输出的是 MatrixEQTL 要求的 SlicedData 格式：
  # MatrixEQTL 要求：第一行 Header (SampleIDs)，第一列 RowNames (VariableNames)。

  # 修正读取逻辑：读取 .align_meta 生成的文件
  # 该文件第一行是 Header (id, Sample1, Sample2...)
  # 实际上我们在内存里统计比较麻烦，因为 MatrixEQTL 会处理 NA。
  # 最准确的方法是模拟 MatrixEQTL 的逻辑：

  # 简单策略：提取两者的共有样本列名，且值不为 NA
  # 因为我们之前的逻辑比较复杂，这里用一种更直接的方法：
  # 分析 p_dt 和 c_dt 的数值部分

  # 提取表型数值向量 (去除第一列变量名)
  # p_dt 只有一行数据 (因为只选了一个表型)
  pheno_vals <- as.numeric(p_dt[1, -1, with=FALSE]) # 第一行，排除第一列ID

  # 提取协变量完整性
  # 只要某一列(样本)在任意协变量上为NA，该样本就会被丢弃
  # c_dt 可能有多行
  covar_mat <- as.matrix(c_dt[, -1, with=FALSE])
  covar_valid <- colSums(is.na(covar_mat)) == 0

  # 表型也必须非 NA
  pheno_valid <- !is.na(pheno_vals)

  # 最终有效样本索引
  valid_idx <- pheno_valid & covar_valid
  final_y <- pheno_vals[valid_idx]

  # 3. 判断数据类型并统计
  n_stat <- list()
  unique_y <- unique(final_y)
  unique_cnt <- length(unique_y)

  # 判定逻辑：如果唯一值数量 <= 5 且包含 0/1 或 1/2，视为二分类
  is_binary <- FALSE
  if (unique_cnt == 2) {
    is_binary <- TRUE
    # 排序：假设小的是Control，大的是Case
    vals <- sort(unique_y)
    n_ctrl <- sum(final_y == vals[1])
    n_case <- sum(final_y == vals[2])

    message(sprintf("   Type: Binary | Control(%s): %d | Case(%s): %d",
                    vals[1], n_ctrl, vals[2], n_case))

    n_stat$N_Control <- n_ctrl
    n_stat$N_Case    <- n_case
  } else {
    # 连续变量
    n_total <- length(final_y)
    message(sprintf("   Type: Continuous | Total Samples: %d", n_total))
    n_stat$N_Total <- n_total
  }

  # 清理内存
  rm(p_dt, c_dt, pheno_vals, covar_mat); gc()

  # =========================================================
  # [结束新增模块]
  # =========================================================

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

  # =========================================================
  # 7. Parallel Physical Chunking Loop
  # =========================================================

  n_chunks <- ceiling(n_snps / chunk_size)
  message(sprintf(">> Starting GWAS in %d chunks (Parallel Cores: %d)...", n_chunks, n_cores))

  # 准备 Header 文件路径 (在并行中只读不写，安全)
  chunk_header_file <- file.path(tmp_dir, "chunk_header.txt")
  writeLines(paste(final_ids, collapse = "\t"), chunk_header_file)

  # --- 定义并行处理函数 ---
  process_chunk <- function(i) {
    # 每个子进程需要独立的临时文件，否则会冲突！
    # 使用 PID 或 index 区分
    my_tmp_raw   <- file.path(tmp_dir, paste0("raw_chunk_", i, ".txt"))
    my_tmp_ready <- file.path(tmp_dir, paste0("ready_chunk_", i, ".txt"))
    my_tmp_res   <- file.path(tmp_dir, paste0("chunk_res_", i, ".txt"))

    start_row <- (i - 1) * chunk_size + 1
    end_row   <- min(i * chunk_size, n_snps)
    n_rows_chunk <- end_row - start_row + 1

    # A. 提取数据 (System Call)
    cmd_extract <- sprintf("tail -n +%d %s | head -n %d > %s",
                           start_row, shQuote(clean_source_file),
                           n_rows_chunk, shQuote(my_tmp_raw))
    system(cmd_extract)

    # B. 合并 Header
    system(paste("cat", shQuote(chunk_header_file), shQuote(my_tmp_raw), ">", shQuote(my_tmp_ready)))

    # C. 加载 SNP Chunk
    # 注意：在 mclapply 中，SlicdeData 对象是 Copy-on-Write 的，读取是安全的
    snps <- MatrixEQTL::SlicedData$new()
    snps$fileDelimiter <- "\t"
    snps$fileOmitCharacters <- "NA"
    snps$fileSkipRows <- 1
    snps$fileSkipColumns <- 0
    snps$fileSliceSize <- 2000
    snps$LoadFile(my_tmp_ready)

    # 设置行名
    rownames(snps) <- bin_ids[start_row:end_row]

    # D. 运行引擎
    # 结果先写到该进程独有的临时文件
    me <- MatrixEQTL::Matrix_eQTL_engine(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name = my_tmp_res,
      pvOutputThreshold = 1.0,
      useModel = MatrixEQTL::modelLINEAR,
      errorCovariance = numeric(),
      verbose = FALSE,
      pvalue.hist = FALSE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE
    )

    # E. 读取结果并返回 Data.Table
    res_chunk <- NULL
    if (file.exists(my_tmp_res)) {
      # 检查文件是否为空 (只有Header算空)
      if (file.size(my_tmp_res) > 50) {
        res_chunk <- data.table::fread(my_tmp_res, header = TRUE)
      }
      file.remove(my_tmp_res) # 读完即删
    }

    # 清理中间文件
    file.remove(my_tmp_raw, my_tmp_ready)
    rm(snps, me); gc(verbose = FALSE)

    return(res_chunk)
  }

  # --- 执行并行 ---
  # 使用 mclapply (Linux only)
  results_list <- parallel::mclapply(1:n_chunks, process_chunk, mc.cores = n_cores)

  # --- 汇总结果 ---
  message(">> Aggregating results from all cores...")

  # 剔除 NULL (失败或空的 Chunk)
  results_list <- Filter(Negate(is.null), results_list)

  if (length(results_list) > 0) {
    final_dt <- data.table::rbindlist(results_list)

    # =========================================================
    # 8. Post-hoc FDR Recalculation (全局 FDR 修正)
    # =========================================================
    message(">> Recalculating Global FDR & Appending Sample Size...")

    # 计算 FDR
    final_dt$FDR <- p.adjust(final_dt$`p-value`, method = "BH")

    # 追加 N (逻辑同之前)
    if (exists("is_binary") && is_binary) {
      final_dt$N_Ctrl <- n_stat$N_Control
      final_dt$N_Case <- n_stat$N_Case
    } else if (exists("n_stat")) {
      final_dt$N_Total <- n_stat$N_Total
    }

    # 写入最终文件
    data.table::fwrite(final_dt, output_file, sep = "\t", quote = FALSE)
    message(sprintf(">> Output written to: %s", output_file))

  } else {
    warning("No results generated from any chunk!")
    file.create(output_file) # 创建空文件防止下游报错
  }

  # 9. Cleanup
  file.remove(clean_source_file, bin_ids_file, chunk_header_file)
  if (file.exists(pheno_ready)) file.remove(pheno_ready)
  if (file.exists(covar_ready)) file.remove(covar_ready)

  message("GWAS Done.")
}