# =============================================================================
# load.R —— 载入 bioinfo-figures 的全部作图函数
# 用法（在项目脚本里，固定好 .libPaths() 之后）：
#   source("<本仓库所在目录>/load.R", encoding = "UTF-8")
# 找自己的位置：先看 source() 传进来的路径（嵌套 source 也行）；不行再看环境变量 BIOFIG_DIR；再看当前工作目录。
# 本文件不改 .libPaths()（由项目的 _common.R 负责）；只检查必需的包，可选包在用到的函数里再查。缺包时不自动安装。
# =============================================================================
FIG_SKILL_VERSION <- "0.2.1 (2026-10-12)"
local({
  is_repo <- function(d) length(d) == 1 && !is.na(d) && nzchar(d) && file.exists(file.path(d, "R", "00_core.R"))
  # source() 把文件路径放在它自己调用帧的 ofile 里；从最里层往外找，第一个是 load.R 且旁边有 R/ 的就是本仓库
  here <- NULL
  for (fr in rev(sys.frames())) {
    f <- get0("ofile", envir = fr, inherits = FALSE)
    if (is.character(f) && length(f) == 1 && basename(f) == "load.R") {
      d <- dirname(normalizePath(f, mustWork = FALSE)); if (is_repo(d)) { here <- d; break }
    }
  }
  if (!is_repo(here)) here <- Sys.getenv("BIOFIG_DIR")
  if (!is_repo(here)) here <- getwd()
  if (!is_repo(here)) stop("找不到 bioinfo-figures 的目录：请写 source(\"<仓库目录>/load.R\")，或先设环境变量 BIOFIG_DIR", call. = FALSE)
  here <- normalizePath(here, winslash = "/")
  need <- c("ggplot2", "scales", "ggrepel", "patchwork")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("bioinfo-figures 需要的包缺失：", paste(miss, collapse = ", "),
                         "。先运行 Rscript \"", file.path(here, "install_deps.R"), "\" 查看缺什么，确认后加 --install 安装；",
                         "或把装有这些包的库加进 .libPaths()。", call. = FALSE)
  files <- sort(list.files(file.path(here, "R"), pattern = "^[0-9]{2}_.*[.]R$", full.names = TRUE))
  for (f in files) sys.source(f, envir = globalenv(), keep.source = FALSE, toplevel.env = globalenv())
  assign("FIG_SKILL_DIR", here, envir = globalenv())
})
