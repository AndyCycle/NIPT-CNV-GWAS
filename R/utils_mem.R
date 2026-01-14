#' Report Current Memory Usage
#'
#' Prints current memory usage (RSS) and triggers garbage collection.
#' @param label A string to identify the checkpoint.
#' @export
report_mem <- function(label = "Checkpoint") {
  # 1. 强制垃圾回收 (虽然会稍微拖慢速度，但在调试期非常必要)
  gc(verbose = FALSE)

  # 2. 获取系统级内存使用 (RSS - Resident Set Size)
  # 使用 ps 命令获取当前进程的内存
  pid <- Sys.getpid()
  tryCatch({
    # Linux ps command to get RSS in KB
    cmd <- sprintf("ps -o rss= -p %d", pid)
    rss_kb <- as.numeric(system(cmd, intern = TRUE))
    rss_gb <- rss_kb / 1024 / 1024

    # 获取 R 内部对象占用 (Ncells + Vcells)
    mem_r <- sum(gc()[, 2]) # MB

    message(sprintf("[MEM] %s | System RSS: %.2f GB | R Objects: %.2f MB",
                    label, rss_gb, mem_r))
  }, error = function(e) {
    message(sprintf("[MEM] %s | (Failed to retrieve system RSS)", label))
  })
}