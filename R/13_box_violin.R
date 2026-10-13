# =============================================================================
# 13_box_violin.R —— 每只鼠（每个样本）一个值的组间比较：箱线版、小提琴版
# 输入 d：panel（分面，如细胞类型；可为因子定顺序）、group（组，如年龄）、value；每行一个样本（不是细胞）
#      s（可选，每个分面一行）：panel、p（按样本的模型给出的 FDR/P）、g1、g2（括号两端的组）
# 用户选择（2026-10-12）：箱线版 = v0.5 现状 box_points_pub 的样子；小提琴版 = gandallab 式；两者都加 FDR 括号和每格 n。
# 科研护栏（代码里强制）：
#   - p 只在两组都 ≥ min_n_p 个样本时显示，且只显示 s 里给的按样本模型的结果，函数不现场做检验（避免逐细胞检验的伪重复）；
#   - 某组 < min_n_box 个样本不画箱、< min_n_violin 个不画小提琴（几个点的"分布"没有意义），只画点；
#   - 每个样本的点都画出来；每格底部写 n。
# 来源：Zhuang 组、人脂肪图谱（箱线 + 每个点）；gandallab/devbrain_xqtl（Nature 2024）小提琴；silverbushlab（Cell 2026，MIT）隔组条纹
# =============================================================================

.box_prep <- function(d, cols) {
  .need_cols(d, c("panel", "group", "value"), "box_pub / violin_pub")
  d$group <- factor(d$group, levels = names(cols))
  if (anyNA(d$group)) stop("box_pub / violin_pub：group 里有不在 cols 名字中的组", call. = FALSE)
  d$panel <- droplevels(if (is.factor(d$panel)) d$panel else factor(d$panel, levels = unique(d$panel)))
  d
}
.shape_ok <- function(d, min_n) {   # 只保留 n ≥ min_n 的组（画汇总形状用）；merge 后重新设因子
  n <- stats::aggregate(value ~ panel + group, data = d, FUN = length)
  out <- merge(d, n[n$value >= min_n, c("panel", "group")])
  out$group <- factor(as.character(out$group), levels = levels(d$group))
  out$panel <- factor(as.character(out$panel), levels = levels(d$panel))
  out
}
.n_layer <- function(d, base_size) {
  n <- stats::aggregate(value ~ panel + group, data = d, FUN = length); n$lab <- paste0("n = ", n$value)
  ggplot2::geom_text(data = n, ggplot2::aes(group, -Inf, label = lab), inherit.aes = FALSE, vjust = -0.6,
                     size = (base_size - 1.5) / ggplot2::.pt, colour = "grey40")
}
.p_layer <- function(d, s, base_size, min_n = 3, lab = "FDR") {
  if (is.null(s) || !nrow(s)) return(NULL)
  .need_cols(s, c("panel", "p", "g1", "g2"), "box_pub / violin_pub 的 s")
  lv <- levels(d$group); rng <- diff(range(d$value, na.rm = TRUE)); tick <- 0.025 * rng
  br <- do.call(rbind, lapply(seq_len(nrow(s)), function(i) {
    z <- d[d$panel == s$panel[i], ]; n1 <- sum(z$group == s$g1[i]); n2 <- sum(z$group == s$g2[i])
    if (is.na(s$p[i]) || n1 < min_n || n2 < min_n) return(NULL)
    data.frame(panel = s$panel[i], x1 = match(s$g1[i], lv), x2 = match(s$g2[i], lv), y = max(z$value, na.rm = TRUE) + 0.08 * rng, p = s$p[i])
  }))
  if (is.null(br)) return(NULL)
  br$panel <- factor(br$panel, levels = levels(d$panel)); br$lab <- fmt_p_math(br$p, lab)
  list(ggplot2::geom_segment(data = br, ggplot2::aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.3),
       ggplot2::geom_segment(data = br, ggplot2::aes(x = x1, xend = x1, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = 0.3),
       ggplot2::geom_segment(data = br, ggplot2::aes(x = x2, xend = x2, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = 0.3),
       ggplot2::geom_text(data = br, ggplot2::aes(x = (x1 + x2) / 2, y = y, label = lab), inherit.aes = FALSE, parse = TRUE,
                          vjust = -0.3, size = (base_size - 1) / ggplot2::.pt))
}
.stripes <- function(lv) ggplot2::geom_rect(data = data.frame(x = seq(2, length(lv), by = 2)),
                                            ggplot2::aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey95")
.box_finish <- function(p, cols, ylab, show_n, show_p, base_size) {
  p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(if (show_n) 0.12 else 0.05, if (show_p) 0.15 else 0.05))) +
    ggplot2::facet_wrap(~panel, nrow = 1) + ggplot2::scale_x_discrete(limits = names(cols)) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") + ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::labs(x = NULL, y = ylab) + theme_pub(base_size) + ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}

#' 箱线版：浅色填充箱（不画离群点）+ 同色描边箱 + 每个样本的点 + FDR 括号 + 每格 n
box_pub <- function(d, s = NULL, cols, ylab = NULL, show_p = TRUE, show_n = TRUE, min_n_box = 3, min_n_p = 3, stripes = FALSE,
                    point_size = 1.1, base_size = FIG_PT, seed = 42L) {
  d <- .box_prep(d, cols); ok <- .shape_ok(d, min_n_box)
  p <- ggplot2::ggplot(d, ggplot2::aes(group, value))
  if (stripes) p <- p + .stripes(names(cols))
  p <- p + ggplot2::geom_boxplot(data = ok, ggplot2::aes(fill = group), alpha = 0.25, colour = NA, outlier.shape = NA, width = 0.6) +
    ggplot2::geom_boxplot(data = ok, ggplot2::aes(colour = group), fill = NA, outlier.shape = NA, width = 0.6, linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = group), size = point_size, shape = 16, position = ggplot2::position_jitter(width = 0.12, height = 0, seed = seed))
  if (show_p) p <- p + .p_layer(d, s, base_size, min_n_p)
  if (show_n) p <- p + .n_layer(d, base_size)
  .box_finish(p, cols, ylab, show_n, show_p, base_size)
}

#' 小提琴版：半透明填充 + 同色描边的小提琴、同色细箱、白色中位数点、每个样本的点；style = "framed" 为 silverbushlab 式（黑细边、深灰窄箱、条纹、四边框）
violin_pub <- function(d, s = NULL, cols, ylab = NULL, show_p = TRUE, show_n = TRUE, min_n_violin = 5, min_n_p = 3, stripes = FALSE,
                       style = c("soft", "framed"), base_size = FIG_PT, seed = 42L) {
  style <- match.arg(style)
  d <- .box_prep(d, cols); ok <- .shape_ok(d, min_n_violin)
  med <- stats::aggregate(value ~ panel + group, data = d, FUN = stats::median)
  p <- ggplot2::ggplot(d, ggplot2::aes(group, value, fill = group, colour = group))
  if (stripes || style == "framed") p <- p + .stripes(names(cols))
  if (style == "soft") {
    p <- p + ggplot2::geom_violin(data = ok, alpha = 0.3, trim = TRUE, linewidth = 0.35, scale = "width", width = 0.8) +
      ggplot2::geom_boxplot(data = ok, width = 0.05, outlier.shape = NA, linewidth = 0.35) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15, height = 0, seed = seed), size = 0.8, shape = 16, alpha = 0.9) +
      ggplot2::geom_point(data = med, shape = 21, fill = "white", size = 1.8, stroke = 0.6)
  } else {
    p <- p + ggplot2::geom_violin(data = ok, colour = "grey20", linewidth = 0.3, scale = "width", width = 0.8) +
      ggplot2::geom_boxplot(data = ok, width = 0.12, fill = "grey30", colour = "grey30", outlier.shape = NA, linewidth = 0.3) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12, height = 0, seed = seed), size = 0.7, shape = 16, colour = "black", alpha = 0.8) +
      ggplot2::geom_point(data = med, shape = 21, fill = "white", colour = "black", size = 1.4, stroke = 0.3)
  }
  if (show_p) p <- p + .p_layer(d, s, base_size, min_n_p)
  if (show_n) p <- p + .n_layer(d, base_size)
  p <- .box_finish(p, cols, ylab, show_n, show_p, base_size)
  if (style == "framed") p <- p + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.5),
                                                 axis.line = ggplot2::element_blank())
  p
}

#' 兼容 v0.5 的旧接口（single-cell-best-practice 的脚本在用）：d 的 x、y 两列给列名；不带 FDR、n。新代码请用 box_pub()
box_points_pub <- function(d, x, y, cols, point_size = 1.1, base_size = FIG_PT, seed = 42L) {
  dd <- d
  dd$x_ <- factor(d[[x]], levels = names(cols)); dd$y_ <- d[[y]]
  ggplot2::ggplot(dd, ggplot2::aes(x_, y_)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = x_), alpha = 0.25, colour = NA, outlier.shape = NA, width = 0.6) +
    ggplot2::geom_boxplot(ggplot2::aes(colour = x_), fill = NA, outlier.shape = NA, width = 0.6, linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(colour = x_), size = point_size, shape = 16, position = ggplot2::position_jitter(width = 0.12, height = 0, seed = seed)) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") + ggplot2::scale_colour_manual(values = cols, guide = "none") +
    theme_pub(base_size) + ggplot2::theme(axis.title.x = ggplot2::element_blank())
}
