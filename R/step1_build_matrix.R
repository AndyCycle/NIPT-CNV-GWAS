#' Build Giant CNV Matrix from WIG Files
#'
#' Stream-processes WIG files, performs HMMcopy correction, applies QC,
#' and appends to a giant matrix on disk.
#'
#' @param input_dir Directory containing WIG files.
#' @param output_prefix Prefix for output files (matrix and sample list).
#' @param gc_file Path to GC wig file.
#' @param map_file Path to Mappability wig file.
#' @param batch_size Number of samples to process in memory before flushing to disk.
#' @param n_cores Number of cores for parallel processing.
#' @param qc_mapd_th Threshold for MAPD (samples > th rejected).
#' @param qc_median_th Threshold for Abs Median (samples > th rejected).
#' @param file_pattern Regex pattern to select specific WIG files.
#' @export
build_cnv_matrix <- function(input_dir, output_prefix, gc_file, map_file,
                             batch_size = 500, n_cores = 10,
                             qc_mapd_th = 0.35, qc_median_th = 0.2,
                             file_pattern = "\\.readcounts\\.wig$") {
  report_mem("Step1: Start")
  # Check system requirements
  if (.Platform$OS.type != "unix") stop("This package requires a Unix-like OS (Linux/macOS) for efficient IO.")

  # Define outputs
  out_mat <- paste0(output_prefix, ".cnv_raw.txt")
  out_samples <- paste0(output_prefix, ".sample_ids.txt")
  out_bins <- paste0(output_prefix, ".bin_coords.txt")

  # File scanning
  # 使用参数 file_pattern 替代原来的硬编码字符串
  wig_files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE, recursive = TRUE)
  if (length(wig_files) == 0) stop("No WIG files found in input_dir.")

  message(sprintf("Found %d samples. Starting build process...", length(wig_files)))

  # Initialize Genome Skeleton (using first sample)
  message("Initializing genomic bins...")
  skeleton <- HMMcopy::wigsToRangedData(wig_files[1], gc_file, map_file, verbose = FALSE)
  bin_ids <- paste0(skeleton$chr, "_", skeleton$start, "_", skeleton$end)

  # Write Bin Info
  data.table::fwrite(data.frame(id = bin_ids), out_mat, sep = "\t", col.names = TRUE)
  data.table::fwrite(data.frame(id = bin_ids, chr = skeleton$chr, start = skeleton$start), out_bins, sep = "\t")

  # Processing Loop
  valid_samples <- character()
  temp_files <- character()

  batches <- split(wig_files, ceiling(seq_along(wig_files)/batch_size))

  report_mem("Step1: Skeleton Built")

  for (i in seq_along(batches)) {
    batch_files <- batches[[i]]
    message(sprintf("Processing Batch %d/%d (%d files)...", i, length(batches), length(batch_files)))

    # Parallel Processing
    batch_res <- parallel::mclapply(batch_files, function(f) {
      tryCatch({
        # 1. HMMcopy Workflow
        uncorr <- HMMcopy::wigsToRangedData(f, gc_file, map_file, verbose = FALSE)
        corr   <- HMMcopy::correctReadcount(uncorr, verbose = FALSE)
        vals   <- corr$copy

        # 2. Quality Control
        clean_vals <- na.omit(vals)
        if (length(clean_vals) < 1000) return(NULL)

        mapd <- .calc_mapd(clean_vals)
        med  <- median(clean_vals)

        if (mapd > qc_mapd_th || abs(med) > qc_median_th) return(NULL)

        # 修改这里：在生成 id 时直接去除后缀
        raw_id <- sub("\\.readcounts\\.wig$", "", basename(f))
        clean_id <- sub("_[0-9]+bp$", "", raw_id) # 提前清洗 _xxxxbp 的窗口后缀（这里涉及原始的wig文件命名，例如原本命名为CL100164113_L01_1_20000bp.readcounts.wig）

        return(list(id = clean_id, data = vals))
      }, error = function(e) NULL)
    }, mc.cores = n_cores)

    # Filter NULLs
    batch_res <- Filter(Negate(is.null), batch_res)

    report_mem(sprintf("Step1: Batch %d/%d Done", i, length(batches)))

    if (length(batch_res) > 0) {
      # Extract Data
      mat_chunk <- do.call(cbind, lapply(batch_res, `[[`, "data"))
      ids_chunk <- vapply(batch_res, `[[`, "id", FUN.VALUE = character(1))

      # Write Temp File
      tmp_f <- paste0(output_prefix, ".tmp_batch_", i, ".txt")
      data.table::fwrite(as.data.table(mat_chunk), tmp_f, sep = "\t", col.names = TRUE) # Header is SampleID

      # Track
      temp_files <- c(temp_files, tmp_f)
      valid_samples <- c(valid_samples, ids_chunk)
    }
  }

  # Final Merge (System Paste)
  message("Merging batches...")
  if (length(temp_files) > 0) {
    # Construct huge paste command safely
    # If too many files, we might need iterative paste, but for <5000 batches single command is usually fine
    cmd <- paste("paste", out_mat, paste(temp_files, collapse = " "), ">", paste0(out_mat, ".tmp"))
    if (system(cmd) != 0) stop("Merge failed.")

    file.rename(paste0(out_mat, ".tmp"), out_mat)
    unlink(temp_files)

    # Write Sample List
    writeLines(valid_samples, out_samples)

    message(sprintf("Success! Matrix: %s", out_mat))
    message(sprintf("Total Valid Samples: %d", length(valid_samples)))
    report_mem("Step1: Merge Done")
  } else {
    stop("No samples passed QC.")
  }
}