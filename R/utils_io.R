# R/utils_io.R

#' Detect columns in a tab-separated file using system tools
#' @noRd
.detect_cols <- function(filepath, row_idx = 1) {
  cmd <- sprintf("head -n %d %s | tail -n 1 | awk -F'\t' '{print NF}'", row_idx, shQuote(filepath))
  as.integer(system(cmd, intern = TRUE))
}

#' Fast column cutter using system 'cut'
#' @noRd
.cut_columns <- function(infile, outfile, start_col, end_col) {
  cmd <- sprintf("cut -f %d-%d %s > %s", start_col, end_col, shQuote(infile), shQuote(outfile))
  if (system(cmd) != 0) stop("System cut command failed.")
}

#' Fast row skipper using system 'tail'
#' @noRd
.skip_rows <- function(infile, outfile, skip = 1) {
  cmd <- sprintf("tail -n +%d %s > %s", skip + 1, shQuote(infile), shQuote(outfile))
  if (system(cmd) != 0) stop("System tail command failed.")
}
