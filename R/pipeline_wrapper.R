#' Smart Pipeline Runner with Resume Capability
#'
#' Sequentially runs the NiptGWAS pipeline. Skips steps if output files already exist.
#'
#' @param work_dir Directory for all outputs.
#' @param prefix File prefix.
#' @param step_args List of arguments for each step.
#' @param force_overwrite Logical. If TRUE, re-runs all steps regardless of existing files.
#' @export
run_pipeline_smart <- function(work_dir, prefix,
                               input_args = list(),
                               force_overwrite = FALSE) {

  full_prefix <- file.path(work_dir, prefix)

  # --- Step 1 Check ---
  raw_mat <- paste0(full_prefix, ".cnv_raw.txt")
  ids_file <- paste0(full_prefix, ".sample_ids.txt")

  if (!force_overwrite && file.exists(raw_mat) && file.exists(ids_file)) {
    message(">> Step 1 output exists. Skipping Matrix Build.")
  } else {
    message(">> Running Step 1: Matrix Build...")
    do.call(build_cnv_matrix, c(list(output_prefix = full_prefix), input_args$step1))
  }

  # --- Step 2 Check ---
  int_mat <- paste0(full_prefix, ".cnv_int.txt")

  if (!force_overwrite && file.exists(int_mat)) {
    message(">> Step 2 output exists. Skipping Normalization.")
  } else {
    message(">> Running Step 2: Normalization...")
    normalize_matrix(input_mat = raw_mat, output_mat = int_mat)
  }

  # --- Step 3 Check ---
  gwas_out <- paste0(full_prefix, ".gwas_results.txt")

  if (!force_overwrite && file.exists(gwas_out)) {
    message(">> Step 3 output exists. Skipping GWAS.")
  } else {
    message(">> Running Step 3: GWAS...")
    do.call(run_nipt_gwas, c(list(
      int_mat_file = int_mat,
      sample_list_file = ids_file,
      output_file = gwas_out,
      tmp_dir = file.path(work_dir, "tmp")
    ), input_args$step3))
  }

  # --- Step 4 Check ---
  plot_file <- paste0(full_prefix, "_Manhattan.png")

  if (!force_overwrite && file.exists(plot_file)) {
      message(">> Step 4 output exists. Skipping Visualization.")
  } else {
      message(">> Running Step 4: Visualization...")
      # 调用上面封装的函数
      visualize_gwas(
          gwas_file = gwas_out,
          output_prefix = full_prefix,
          title = input_args$step4$title %||% "CNV-GWAS Analysis" # %||% 需自定义或用 ifelse
      )
  }

  message(">> Pipeline finished.")
}

# 辅助小函数，防止 input_args$step4$title 为 NULL 时报错
%||%` <- function(a, b) if (!is.null(a)) a else b