# =============================================================================
# 20_gsea_pca.R —— GSEA 三联图、bulk / pseudobulk PCA（规范版，不做盲选）
# 来源（读原代码后重写）：orca-hd/trop2_in_crc、enrichplot gseaplot2 的三联布局；taiji2（形状 = 来源、颜色 = 分组）
# =============================================================================

#' GSEA 三联图：富集分数曲线 + 命中竖线 + 排序指标（底部蓝白红色带）；峰值处虚线；图上写 NES 和 FDR。
#' ranks：命名的数值向量（基因 → 排序量，如 sign(logFC) × −log10 P）；set：基因集
gsea_curve_pub <- function(ranks, set, title, nes, fdr, rank_label = "Gene rank", base_size = FIG_PT, col = "#2166AC") {
  r <- sort(ranks, decreasing = TRUE); hit <- names(r) %in% set
  if (!any(hit)) stop("gsea_curve_pub：基因集里没有一个基因在 ranks 中", call. = FALSE)
  w <- abs(r) * hit; es <- cumsum(w / sum(w) - (!hit) / sum(!hit)); pk <- which.max(abs(es))
  d <- data.frame(rank = seq_along(r), es = es, hit = hit, metric = r)
  pos <- es[pk] > 0
  p1 <- ggplot2::ggplot(d, ggplot2::aes(rank, es)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = pk, colour = "grey50", linetype = "dashed", linewidth = 0.25) +
    ggplot2::geom_line(colour = col, linewidth = 0.6) +
    ggplot2::annotate("text", x = if (pos) nrow(d) * 0.97 else nrow(d) * 0.03, y = if (pos) max(es) else min(es),
                      label = sprintf("NES = %.2f\nFDR = %s", nes, fmt_p_txt(fdr)), hjust = if (pos) 1 else 0, vjust = if (pos) 1 else 0,
                      size = (base_size - 0.5) / ggplot2::.pt) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::labs(x = NULL, y = "Enrichment score", title = title) + theme_pub(base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank())
  p2 <- ggplot2::ggplot(d[d$hit, ], ggplot2::aes(x = rank, xend = rank, y = 0, yend = 1)) + ggplot2::geom_segment(linewidth = 0.2, colour = "grey15") +
    ggplot2::scale_x_continuous(limits = c(1, nrow(d)), expand = c(0, 0)) + ggplot2::theme_void()
  rng <- diff(range(d$metric)); ybar <- min(d$metric) - 0.06 * rng
  p3 <- ggplot2::ggplot(d, ggplot2::aes(rank, metric)) + ggplot2::geom_area(fill = "grey75") +
    ggplot2::geom_tile(ggplot2::aes(y = ybar, fill = metric), height = 0.08 * rng) +
    scale_div0(stats::quantile(abs(d$metric), 0.98), aesthetics = "fill", guide = "none") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::labs(x = rank_label, y = "Ranking metric") + theme_pub(base_size)
  patchwork::wrap_plots(p1, p2, p3, ncol = 1, heights = c(3, 0.5, 1.4))
}

#' PCA：logcpm（基因 × 样本，已 log 化），info（每个样本一行，顺序与列一致）；颜色 = colour 列，形状 = shape 列；轴标题写方差占比
pca_pub <- function(logcpm, info, colour, shape = NULL, cols, n_top = 2000, title = NULL, base_size = FIG_PT) {
  .need_cols(info, c(colour, shape), "pca_pub")
  if (ncol(logcpm) != nrow(info)) stop("pca_pub：logcpm 的列数和 info 的行数不一致", call. = FALSE)
  v <- apply(logcpm, 1, stats::var); x <- logcpm[order(-v)[seq_len(min(n_top, length(v)))], ]
  pc <- stats::prcomp(t(x), scale. = FALSE); ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  d <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], col = info[[colour]], shp = if (is.null(shape)) "all" else info[[shape]])
  ggplot2::ggplot(d, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_point(ggplot2::aes(colour = col, shape = shp), size = 1.8, stroke = 0.6) +
    ggplot2::scale_colour_manual(values = cols, name = colour) +
    ggplot2::scale_shape_manual(values = c(16, 1, 17, 2, 15, 0), name = shape %||% "", guide = if (is.null(shape)) "none" else "legend") +
    ggplot2::labs(x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]), title = title) + theme_pub(base_size) +
    ggplot2::theme(aspect.ratio = 1)
}
