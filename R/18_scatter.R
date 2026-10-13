# =============================================================================
# 18_scatter.R —— 样本层面的相关散点（每个点 = 一个样本/一只鼠，不是细胞）
# 输入 d：x、y、group（组，如年龄）；dataset（可选，多个数据集时用形状区分）
# 用户选择（2026-10-12，第二批盲选）：style = "black"——黑色回归线 + 灰色 95% 置信带、点按组浅色填充、R / P / n 写在面板上方（不压点）。
# 开关（缝合体）：style = "group"（点按组着色 + 形状 = 数据集 + 灰色虚线回归线；数据集 > 1 时建议用）、"red"（silverbushlab 式红线、统计写在面板内最空的角）
# 科研护栏：相关按样本算（cor.test），n 写进图里；多个数据集混在一起时可能有辛普森悖论，图注要说明是否分数据集看过。
# 来源（读原代码后重写）：silverbushlab ggscatter（Cell 2026，MIT）；aifimmunology
# =============================================================================

.cor_lab <- function(d, method) {
  ct <- stats::cor.test(d$x, d$y, method = method)
  sprintf("italic(%s) == %.2f*','~%s*','~italic(n) == %d", if (method == "pearson") "R" else "rho", ct$estimate, fmt_p_math(ct$p.value), nrow(d))
}
.cor_annot <- function(d, method, base_size) {   # 写在点最少的那个角
  qx <- stats::quantile(d$x, c(0.3, 0.7)); qy <- stats::quantile(d$y, c(0.3, 0.7))
  n <- c(tl = sum(d$x <= qx[2] & d$y >= qy[2]), tr = sum(d$x >= qx[1] & d$y >= qy[2]), bl = sum(d$x <= qx[2] & d$y <= qy[1]), br = sum(d$x >= qx[1] & d$y <= qy[1]))
  k <- names(which.min(n))
  ggplot2::annotate("text", x = if (k %in% c("tl", "bl")) -Inf else Inf, y = if (k %in% c("tl", "tr")) Inf else -Inf, label = .cor_lab(d, method), parse = TRUE,
                    hjust = if (k %in% c("tl", "bl")) -0.06 else 1.06, vjust = if (k %in% c("tl", "tr")) 1.4 else -0.6, size = base_size / ggplot2::.pt)
}

#' 样本层面相关散点。gcols：组 → 颜色；method = "pearson" / "spearman"（n 小或有离群点时建议同时看 Spearman）
scatter_pub <- function(d, xlab, ylab, gcols = NULL, style = c("black", "group", "red"), method = c("pearson", "spearman"), base_size = FIG_PT) {
  style <- match.arg(style); method <- match.arg(method)
  .need_cols(d, c("x", "y"), "scatter_pub")
  if (style != "red") .need_cols(d, "group", "scatter_pub")
  sub <- ggplot2::theme(plot.subtitle = ggplot2::element_text(size = base_size))
  if (style == "black") {
    return(ggplot2::ggplot(d, ggplot2::aes(x, y)) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey85", linewidth = 0.4) +
      ggplot2::geom_point(ggplot2::aes(fill = group), size = 1.6, shape = 21, colour = "grey20", stroke = 0.25) +
      ggplot2::scale_fill_manual(values = gcols, name = NULL) +
      ggplot2::labs(x = xlab, y = ylab, subtitle = parse(text = .cor_lab(d, method))[[1]]) + theme_pub(base_size) +
      ggplot2::theme(legend.position = "bottom") + sub)
  }
  if (style == "group") {
    if (!"dataset" %in% names(d)) d$dataset <- "all"
    return(ggplot2::ggplot(d, ggplot2::aes(x, y)) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = "grey45", fill = "grey90", linewidth = 0.35, linetype = "dashed") +
      ggplot2::geom_point(ggplot2::aes(colour = group, shape = dataset), size = 1.5, stroke = 0.6) +
      ggplot2::scale_colour_manual(values = gcols, name = NULL) + ggplot2::scale_shape_manual(values = c(16, 1, 17, 2, 15, 0), name = "Dataset") +
      ggplot2::labs(x = xlab, y = ylab, subtitle = parse(text = .cor_lab(d, method))[[1]]) + theme_pub(base_size) + sub)
  }
  ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = "red", fill = "lightgray", linewidth = 0.5) +
    ggplot2::geom_point(size = 1.1, shape = 16, colour = "black") + .cor_annot(d, method, base_size) +
    ggplot2::labs(x = xlab, y = ylab) + theme_pub(base_size)
}
