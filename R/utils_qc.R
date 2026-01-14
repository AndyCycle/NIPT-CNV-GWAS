# R/utils_qc.R

#' Calculate Median Absolute Pairwise Difference
#' @noRd
.calc_mapd <- function(x) {
  # Robust implementation handling NAs
  x <- na.omit(x)
  if (length(x) < 2) return(Inf)
  median(abs(diff(x)))
}