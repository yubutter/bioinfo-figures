# =============================================================================
# 00_core.R —— 共用：版面尺寸、字号、主题、颜色、栅格化、存图、p 值格式、列检查
# 来源：single-cell-best-practice v0.5 的 S10 作图函数（sc_helpers.R），2026-10-12 作图调研二期后搬到本 skill 并扩充。
# 规则见 references/conventions.md。
# =============================================================================

FIG_MM <- c(single = 89, onehalf = 120, double = 183, word = 155)  # 版面宽度（mm）：Nature 单栏/1.5 栏/双栏；Word 版心
FIG_PT <- 7                                                         # 成图字号（pt）；期刊要求 5–7 pt
FONT_PUB <- Sys.getenv("BIOFIG_FONT", "Arial")                      # 报告档/投稿档字体（无衬线，全文统一）；没有 Arial 的机器（常见于 Linux）可设环境变量 BIOFIG_FONT
OKABE_ITO <- c(orange = "#E69F00", sky_blue = "#56B4E9", bluish_green = "#009E73", yellow = "#F0E442",
               blue = "#0072B2", vermillion = "#D55E00", reddish_purple = "#CC79A7", grey = "#999999")
DIR_COLS <- c(up = "#D6604D", down = "#4393C3", ns = "grey85")      # 上调/下调/不显著（柔和红蓝，不用纯红纯蓝）

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- 检查 ----------------------------------------------------------------------------------------------------
#' 输入表必须有的列；缺了就停下并用中文列出（模板的"数据契约"）
.need_cols <- function(d, cols, fun) {
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop(sprintf("%s：输入表缺少列 %s（需要：%s）", fun, paste(miss, collapse = ", "), paste(cols, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}
#' 可选包：没有就停下，给出安装命令，不自动安装。bioc = TRUE 表示是 Bioconductor 的包
.need_pkg <- function(pkg, bioc = FALSE) {
  for (p in pkg) if (!requireNamespace(p, quietly = TRUE)) {
    how <- if (bioc) sprintf("BiocManager::install(\"%s\")", p) else sprintf("install.packages(\"%s\")", p)
    stop(sprintf("缺少 R 包 %s。本模板库不自动安装：可以运行 %s（或运行仓库里的 install_deps.R 查看全部缺包），也可以把装有它的库加进 .libPaths()",
                 p, how), call. = FALSE)
  }
  invisible(TRUE)
}

# ---- 颜色 ----------------------------------------------------------------------------------------------------
#' 把颜色和另一种颜色按比例 p 混合（p = 0 不变，p = 1 变成 with）；用来造同色相的深浅、给描边加深、把面积色调淡
mix_col <- function(col, with, p) {
  out <- vapply(seq_along(col), function(i) {
    a <- grDevices::col2rgb(col[i]) / 255; b <- grDevices::col2rgb(with) / 255
    grDevices::rgb(t(a * (1 - p[min(i, length(p))]) + b * p[min(i, length(p))]))
  }, character(1))
  stats::setNames(out, names(col))
}

#' 谱系色族：每个谱系一个主色，谱系内的细胞类型保持这个色相、按亮度由深到浅分档（HCL 空间，相邻两档亮度差 l_step）。
#' lineage：命名字符向量，names = 细胞类型（按想要的顺序），值 = 所属谱系；base_cols：每个谱系的主色
#' 依据：Zhuang 组衰老脑（大类一个色相、亚群取深浅）；Lake 2023 人肾图谱（按谱系成色族）。在 HCL 里调亮度而不是混黑/混白：混黑发闷发脏。
lineage_palette <- function(lineage, base_cols, l_step = 14, l_lim = c(33, 88)) {
  miss <- setdiff(unique(lineage), names(base_cols))
  if (length(miss)) stop("缺少谱系主色：", paste(miss, collapse = ", "))
  out <- character(0)
  for (lg in unique(lineage)) {
    types <- names(lineage)[lineage == lg]; n <- length(types)
    if (n == 1) {
      cols <- base_cols[[lg]]
    } else {
      luv <- grDevices::convertColor(t(grDevices::col2rgb(base_cols[[lg]]) / 255), from = "sRGB", to = "Luv")
      L <- luv[1]; C <- sqrt(luv[2]^2 + luv[3]^2); H <- (atan2(luv[3], luv[2]) * 180 / pi) %% 360
      span <- l_step * (n - 1)
      lo <- max(l_lim[1], min(L - span / 2, l_lim[2] - span))
      ls <- lo + l_step * (0:(n - 1))
      cs <- C * (1 - 0.25 * (ls - lo) / span)          # 变浅时略降饱和度，避免超出色域
      cols <- grDevices::hcl(H, cs, ls, fixup = TRUE)
    }
    out <- c(out, stats::setNames(cols, types))
  }
  out[names(lineage)]
}

#' 类别色：≤ 8 组用 Okabe–Ito（色盲友好）；9–32 组加入其加深 25% / 变浅 30% / 加深 45% 版本。细胞类型多于 8 种时优先用 lineage_palette()
cat_palette <- function(n, names = NULL) {
  base <- unname(OKABE_ITO)
  pool <- c(base, mix_col(base, "black", 0.25), mix_col(base, "white", 0.30), mix_col(base, "black", 0.45))
  if (n > length(pool)) stop("类别多于 ", length(pool), " 个：请按谱系分色族（lineage_palette）或合并类别")
  out <- unname(pool[seq_len(n)])
  if (!is.null(names)) out <- stats::setNames(out, names)
  out
}

#' 以 0 为中心的对称发散色标：只给有正负的量（得分、log2FC、z 分数、NES）用；超出 ±lim 压到端点色（squish）
scale_div0 <- function(lim, aesthetics = "colour", low = "#2166AC", mid = "white", high = "#B2182B", ...) {
  ggplot2::scale_colour_gradient2(low = low, mid = mid, high = high, midpoint = 0, limits = c(-lim, lim),
                                  oob = scales::squish, aesthetics = aesthetics, ...)
}
div0_limit <- function(x, cap = Inf) min(max(abs(x), na.rm = TRUE), cap)

# ---- 主题 ----------------------------------------------------------------------------------------------------
#' 报告/投稿档主题：细线（0.25）、小字号、黑色坐标文字
theme_pub <- function(base_size = FIG_PT, base_family = "") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      line = ggplot2::element_line(linewidth = 0.25),
      axis.line = ggplot2::element_line(linewidth = 0.25), axis.ticks = ggplot2::element_line(linewidth = 0.25),
      axis.ticks.length = ggplot2::unit(1, "mm"), axis.text = ggplot2::element_text(colour = "black", size = base_size),
      axis.title = ggplot2::element_text(size = base_size), legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size), legend.key.size = ggplot2::unit(3, "mm"),
      strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold"),
      plot.margin = ggplot2::margin(2, 2, 2, 2, "mm"))
}
theme_umap_pub <- function(base_size = FIG_PT, base_family = "") {
  ggplot2::theme_void(base_size = base_size, base_family = base_family) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = base_size), legend.title = ggplot2::element_text(size = base_size),
                   legend.key.size = ggplot2::unit(3, "mm"), plot.title = ggplot2::element_text(size = base_size + 1, face = "bold"),
                   strip.text = ggplot2::element_text(size = base_size, face = "bold"),
                   plot.margin = ggplot2::margin(1, 1, 1, 1, "mm"))
}

#' 只把点层栅格化（ggrastr），文字和坐标轴保持矢量；dpi = NULL 或没装 ggrastr 时不栅格化
.rast <- function(layer, dpi) {
  if (is.null(dpi) || !requireNamespace("ggrastr", quietly = TRUE)) return(layer)
  ggrastr::rasterise(layer, dpi = dpi, dev = "cairo")
}

# ---- 文字 ----------------------------------------------------------------------------------------------------
#' 斜体基因名（plotmath）。必须写成带引号的字符串：italic(Gapdh) 会把字母逐个排开
lab_italic <- function(gene, rest = NULL) {
  txt <- sprintf("italic('%s')", gene)
  if (!is.null(rest)) txt <- paste0(txt, "~", rest)
  parse(text = txt)[[1]]
}
#' 例：lab_gene_lfc("Gapdh") → Gapdh（斜体）log2 fold change (old vs young)
lab_gene_lfc <- function(gene, contrast = "old vs young") lab_italic(gene, sprintf("log[2]~'fold change (%s)'", contrast))

#' 精确 p 值的 plotmath 文本（parse = TRUE 用）：≥ 0.001 写两位有效数字，更小写成 m × 10^e。lab 是前缀，如 "italic(P)" 或 "FDR"
fmt_p_math <- function(p, lab = "italic(P)") {
  vapply(p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x >= 0.001) return(sprintf("%s == '%s'", lab, formatC(signif(x, 2), format = "fg", digits = 2)))
    e <- floor(log10(x)); m <- x / 10^e
    sprintf("%s == %s %%*%% 10^%d", lab, sprintf("%.1f", m), e)
  }, "")
}
#' 普通文本版（表格列、图注）：0.0085、2.4e-4
fmt_p_txt <- function(p) ifelse(is.na(p), "", ifelse(p >= 0.001, formatC(signif(p, 2), format = "fg", digits = 2),
                                                       sub("e-0*", "e-", formatC(p, format = "e", digits = 1))))

#' 通路名清理：去掉 HALLMARK_ / REACTOME_ 等前缀，下划线变空格，句首大写（保留常见缩写）。Reactome 优先用 msigdbr 的 gs_description（原始大小写）
clean_term <- function(x) {
  x <- sub("^(HALLMARK|REACTOME|KEGG|GOBP|GOMF|GOCC|WP|BIOCARTA|PID)_", "", x)
  x <- tolower(gsub("_", " ", x))
  keep <- c("dna", "rna", "mtorc1", "tnfa", "nfkb", "kras", "il6", "jak", "stat3", "il2", "stat5", "tgf", "e2f", "g2m", "uv",
            "p53", "myc", "pi3k", "akt", "mtor", "mrna", "trna", "rrna", "atp", "nmd", "ecm", "nad", "nadh", "gpcr", "mapk",
            "erk", "hif", "ifn", "tlr", "nk", "tcr", "bcr", "il", "cd28", "ii", "iii", "iv", "srp", "fgfr", "egfr")
  w <- strsplit(x, " ")
  x <- vapply(w, function(v) paste(ifelse(v %in% keep, toupper(v), v), collapse = " "), "")
  x <- paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  gsub("\\bRRNA\\b", "rRNA", gsub("\\bTRNA\\b", "tRNA", gsub("\\bMRNA\\b", "mRNA", x)))   # 句首大写之后再改，否则又被改回 RRNA
}

# ---- 溯源与冻结（项目用"指针"方式载入本 skill 时用；2026-10-12 用户选方案 A）--------------------------------------
#' 本 skill 的版本号和 git 提交号，写进项目的 software_versions.csv
fig_skill_provenance <- function() {
  dir <- if (exists("FIG_SKILL_DIR")) FIG_SKILL_DIR else Sys.getenv("BIOFIG_DIR", NA_character_)
  if (is.na(dir) || !nzchar(dir)) stop("fig_skill_provenance：没有通过 load.R 载入，也没设 BIOFIG_DIR，不知道模板库在哪", call. = FALSE)
  head <- file.path(dir, ".git", "HEAD"); sha <- NA_character_
  if (file.exists(head)) {
    h <- trimws(readLines(head, warn = FALSE)[1])
    if (startsWith(h, "ref: ")) {
      ref <- sub("^ref: ", "", h); f <- file.path(dir, ".git", ref)
      if (file.exists(f)) sha <- trimws(readLines(f, warn = FALSE)[1]) else {
        pk <- file.path(dir, ".git", "packed-refs")
        if (file.exists(pk)) { l <- grep(paste0(" ", ref, "$"), readLines(pk, warn = FALSE), value = TRUE); if (length(l)) sha <- sub(" .*", "", l[1]) }
      }
    } else sha <- h
  }
  data.frame(software = "bioinfo-figures", version = if (exists("FIG_SKILL_VERSION")) FIG_SKILL_VERSION else NA_character_,
             commit = sha, path = dir)
}

#' 交付或投稿前把当时的作图代码复制进项目（冻结），以后 skill 再改，这个项目的图也能原样重画
#' dest：项目里的目录（如 "R/figures_frozen"）；已存在就拒绝覆盖
fig_skill_freeze <- function(dest) {
  if (dir.exists(dest)) stop("fig_skill_freeze：", dest, " 已存在，不覆盖", call. = FALSE)
  pv <- fig_skill_provenance(); dir.create(file.path(dest, "R"), recursive = TRUE)
  file.copy(list.files(file.path(pv$path, "R"), full.names = TRUE), file.path(dest, "R"))
  file.copy(file.path(pv$path, "load.R"), dest)
  utils::write.csv(pv, file.path(dest, "FROZEN_FROM.csv"), row.names = FALSE)
  message("作图代码已冻结到 ", dest, "（", pv$version, "，", pv$commit, "）；项目改为 source(file.path(\"", dest, "\", \"load.R\"))")
  invisible(dest)
}

# ---- 存图 ----------------------------------------------------------------------------------------------------
#' 按最终尺寸（mm）存 PDF（cairo，Arial，文字保持为文字）和 PNG；p 可以是 ggplot / patchwork / ComplexHeatmap 对象
#'   data：作图数据（data.frame），给了就另存 <stem>_data.csv（UTF-8）——项目约定"每张图都有作图原始数据"
#'   legend_separate = TRUE：另存"无图例版"和"单独的图例 PDF"，方便拼版
save_plot_pub <- function(p, dir, stem, width_mm, height_mm, dpi = 600, family = FONT_PUB, legend_separate = FALSE, data = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  w <- width_mm / 25.4; h <- height_mm / 25.4
  draw <- function(obj) if (inherits(obj, c("Heatmap", "HeatmapList"))) ComplexHeatmap::draw(obj, merge_legend = TRUE) else print(obj)
  grDevices::cairo_pdf(file.path(dir, paste0(stem, ".pdf")), width = w, height = h, family = family)
  draw(p); grDevices::dev.off()
  grDevices::png(file.path(dir, paste0(stem, ".png")), width = w * dpi, height = h * dpi, res = dpi,
                 type = "cairo", bg = "white", family = family)
  draw(p); grDevices::dev.off()
  if (legend_separate && inherits(p, "ggplot")) {
    .need_pkg("cowplot")
    leg <- cowplot::get_legend(p + ggplot2::theme(legend.position = "right"))
    grDevices::cairo_pdf(file.path(dir, paste0(stem, "_legend.pdf")), width = 40 / 25.4, height = h, family = family)
    grid::grid.newpage(); grid::grid.draw(leg); grDevices::dev.off()
    grDevices::cairo_pdf(file.path(dir, paste0(stem, "_nolegend.pdf")), width = w, height = h, family = family)
    print(p + ggplot2::theme(legend.position = "none")); grDevices::dev.off()
  }
  if (!is.null(data)) {
    con <- file(file.path(dir, paste0(stem, "_data.csv")), open = "w", encoding = "UTF-8")
    utils::write.csv(as.data.frame(data), con, row.names = FALSE); close(con)
  }
  msg <- paste0("保存图 ", file.path(dir, stem), "（", width_mm, " × ", height_mm, " mm）")
  if (exists("log_msg", mode = "function")) log_msg(msg) else message(msg)
  invisible(stem)
}
