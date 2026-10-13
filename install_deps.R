# =============================================================================
# install_deps.R —— 检查并安装 bioinfo-figures 需要的 R 包
# 默认只检查、不安装：列出缺哪些包、各自影响哪些函数。确认后加 --install 才安装。
# 运行：
#   Rscript install_deps.R                       # 只检查
#   Rscript install_deps.R --install             # 安装缺的必需包和可选包
#   Rscript install_deps.R --install --required  # 只装必需包
#   可加 --lib=<目录>（装到指定库，默认 .libPaths() 的第一个）、--repo=<CRAN 镜像>（如 https://mirrors.tuna.tsinghua.edu.cn/CRAN/）
# 包库：设了环境变量 BIOFIG_LIBS（分号分隔）就先用它；否则仓库根目录有 local_libpaths.R（本机专用，不进仓库）就执行它。
# 测试过的版本：R 4.5.3；ggplot2 4.0.3、scales 1.4.0、ggrepel 0.9.8、patchwork 1.3.2、ggrastr 1.0.2、cowplot 1.2.0、colorspace 2.1.3、
#   circlize 0.4.18、ComplexUpset 1.3.3、ggalluvial 0.12.6、ComplexHeatmap 2.26.1。
# =============================================================================
args <- commandArgs(TRUE)
opt <- function(name) { a <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(a)) sub("^[^=]*=", "", a[1]) else "" }
install <- "--install" %in% args; required_only <- "--required" %in% args

here <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
})
libs <- Sys.getenv("BIOFIG_LIBS", "")
if (nzchar(libs)) { libs <- strsplit(libs, ";")[[1]]; .libPaths(libs[dir.exists(libs)])
} else if (file.exists(file.path(here, "local_libpaths.R"))) source(file.path(here, "local_libpaths.R"), encoding = "UTF-8")

PKGS <- data.frame(
  pkg  = c("ggplot2", "scales", "ggrepel", "patchwork",
           "ggrastr", "cowplot", "colorspace", "circlize", "ComplexUpset", "ggalluvial", "ComplexHeatmap"),
  need = c(rep("必需", 4), rep("可选", 7)),
  src  = c(rep("CRAN", 10), "Bioconductor"),
  used = c("所有图", "所有图", "所有带标签的图", "拼图（森林图表格、点图细胞数条等）",
           "大数据点层栅格化（没有时点层保持矢量，PDF 会很大）", "save_plot_pub(legend_separate = TRUE)",
           "palette_check() 的色盲模拟", "heatmap_pub", "upset_pub（默认引擎）", "sankey_pub", "heatmap_pub"),
  stringsAsFactors = FALSE)
if (getRversion() < "4.1.0") stop("需要 R ≥ 4.1（测试环境是 R 4.5.3），当前是 ", getRversion(), call. = FALSE)
PKGS$have <- vapply(PKGS$pkg, function(p) nzchar(system.file(package = p)), logical(1))
PKGS$version <- vapply(PKGS$pkg, function(p) if (nzchar(system.file(package = p))) as.character(utils::packageVersion(p)) else "", "")

cat("R ", as.character(getRversion()), "；包库：", paste(.libPaths(), collapse = " ; "), "\n\n", sep = "")
for (i in seq_len(nrow(PKGS))) with(PKGS[i, ], cat(sprintf("  %-4s %-15s %-13s %-10s %s\n", need, pkg, src,
                                                        if (have) version else "【缺】", used)))
if (PKGS$have[PKGS$pkg == "ggplot2"] && utils::packageVersion("ggplot2") < "4.0.0")
  cat("\n注意：ggplot2 低于 4.0（测试用的是 4.0.3），多数图应该能画，但没有测过。\n")

todo <- PKGS[!PKGS$have & (PKGS$need == "必需" | !required_only), ]
if (!nrow(todo)) { cat("\n需要的包都在。\n"); quit(save = "no") }
# 默认装到个人库（R_LIBS_USER 的第一个）：刚装好 R 的电脑上个人库目录还不存在，.libPaths()[1] 会是 R 自带的库，通常没有写权限
lib <- opt("lib")
if (!nzchar(lib)) lib <- strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep, fixed = TRUE)[[1]][1]
if (is.na(lib) || !nzchar(lib)) lib <- .libPaths()[1]
lib <- normalizePath(path.expand(lib), winslash = "/", mustWork = FALSE)
cat("\n缺：", paste(todo$pkg, collapse = ", "), "\n将安装到：", lib, if (!dir.exists(lib)) "（目录还不存在，安装时新建）", "\n", sep = "")
if (!install) {
  cat("这次只检查，没有安装任何东西。确认后运行：Rscript install_deps.R --install",
      if (any(todo$need == "可选")) "（只装必需包再加 --required）" else "", "\n", sep = "")
  quit(save = "no")
}

dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))   # 新装的 BiocManager 和最后的核对都要能在这个库里找到包
repo <- opt("repo"); if (!nzchar(repo)) { repo <- getOption("repos")[["CRAN"]]; if (is.null(repo) || repo == "@CRAN@") repo <- "https://cloud.r-project.org" }
cat("安装到：", lib, "；CRAN 镜像：", repo, "\n", sep = "")
cran <- todo$pkg[todo$src == "CRAN"]; bioc <- todo$pkg[todo$src == "Bioconductor"]
if (length(cran)) utils::install.packages(cran, lib = lib, repos = repo)
if (length(bioc)) {
  if (!nzchar(system.file(package = "BiocManager"))) utils::install.packages("BiocManager", lib = lib, repos = repo)
  BiocManager::install(bioc, lib = lib, update = FALSE, ask = FALSE)
}
left <- todo$pkg[!vapply(todo$pkg, function(p) nzchar(system.file(package = p)), logical(1))]
if (length(left)) stop("这些包没装上：", paste(left, collapse = ", "), "（看上面的安装日志）", call. = FALSE)
cat("\n装好了。下一步：Rscript --vanilla tests/test_templates.R，应该显示 ALL PASS。\n")
