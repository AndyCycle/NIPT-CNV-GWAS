#' Apply Row-wise Inverse Normal Transformation (INT)
#'
#' Reads a giant matrix row by row, applies INT, and writes to a new file.
#' Extremely memory efficient.
#'
#' @param input_mat Path to the raw matrix (from step 1).
#' @param output_mat Path for the normalized matrix.
#' @param k Blom's constant (default 0.375).
#' @export
normalize_matrix <- function(input_mat, output_mat, k = 0.375) {
  report_mem("Step2: Start")
  message(sprintf("Starting INT normalization on: %s", input_mat))

  con_in <- file(input_mat, "r")
  con_out <- file(output_mat, "w")
  on.exit({ close(con_in); close(con_out) })

  # Copy Header
  header_line <- readLines(con_in, n = 1)
  if (length(header_line) == 0) stop("Empty input matrix")
  writeLines(header_line, con_out) # 原样写入 Header (保留样本ID)

  cnt <- 0

  repeat {
    lines <- readLines(con_in, n = 1000) # Read in chunks of 1000 lines for speed
    if (length(lines) == 0) break

    # Process Chunk
    processed_lines <- vapply(lines, function(line) {
      parts <- strsplit(line, "\t", fixed = TRUE)[[1]] # 使用 fixed=TRUE 加速
      bin_id <- parts[1]

      # 注意：parts[1] 是行名，parts[-1] 是数据
      # 必须确保 parts 长度一致。如果某行数据损坏（比如最后少个TAB），as.numeric会产生NA
      vals <- as.numeric(parts[-1])

      # INT Logic
      n_obs <- length(vals)
      valid_idx <- which(!is.na(vals))
      n_valid <- length(valid_idx)

      if (n_valid > 10) {
        ranks <- rank(vals[valid_idx], ties.method = "average")
        transformed <- qnorm((ranks - k) / (n_valid - 2*k + 1))
        # Rounding to save disk space and IO time
        vals[valid_idx] <- round(transformed, 6)
      }

      # Reconstruct line
      # Handle NAs explicitly for output
      vals_str <- as.character(vals)
      vals_str[is.na(vals)] <- "NA"

      # ++++++ 更稳健的拼接 ++++++
      paste(c(bin_id, vals_str), collapse = "\t")
      # ++++++++++++++++++++++++++++++++
    }, FUN.VALUE = character(1), USE.NAMES = FALSE)

    writeLines(processed_lines, con_out)
    cnt <- cnt + length(lines)
    if (cnt %% 10000 == 0) {
      report_mem("Step2: Start")
      message(sprintf("Processed %d bins...", cnt))
    }
  }

  message("Normalization complete.")
}