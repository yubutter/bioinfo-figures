# =============================================================================
# 16_enrich.R —— 富集图：单个比较（enrich_pub）、多组网格（enrich_grid_pub）
# enrich_pub 输入 d：term（通路名，已清理，见 clean_term()；Reactome 建议用 msigdbr 的 gs_description）、nes（有正负）、fdr、n_le（leading edge 基因数）
# enrich_grid_pub 输入 d：row（通路）、col（细胞类型/条件）、nes、fdr；行列顺序按因子
# 用户选择（2026-10-12）：
#   单个比较 = v = "fdr_bar_dot"：x = −log10 FDR，按方向着色的浅色横条 + 末端空心点（大小 = 基因数）+ 条外写 NES；
#     点的大小按数据范围映射到 size_range（用户："最大的不过大、最小的不过小"）；
#   多组网格 = style = "dot"：圆点（大小 = −log10 FDR，填充 = NES 蓝白红）+ 显著格细黑圈（用户："边框稍微细一点"，0.35）。
# 开关（缝合体）：enrich_pub v = "nes_bar_dot" / "lollipop" / "dot"；enrich_grid_pub style = "tile"（方格 + 黑框）/ "triangle"（三角 = 方向）
# 科研要求（图注要写）：基因集来源和版本、背景基因集（ORA 时）、多重校正方法；上调和下调分开做或用 GSEA；除了 p 值同时给出效应（NES）。
# 来源（读原代码后重写）：linnarsson go_plot 横条、totalx 棒棒糖、heiniglab 心脏图谱网格黑圈（Nature 2022）、allonkleinlab 显著黑框
# =============================================================================

.enrich_prep <- function(d, wrap, dir_labels) {
  .need_cols(d, c("term", "nes", "fdr", "n_le"), "enrich_pub")
  d$dir <- factor(ifelse(d$nes > 0, dir_labels[1], dir_labels[2]), levels = dir_labels)
  d$term_w <- vapply(d$term, function(x) paste(strwrap(x, wrap), collapse = "\n"), "")
  d$lfdr <- -log10(d$fdr)
  d
}

#' 单个比较的富集图。dir_labels：NES > 0 / < 0 时的说法，如 c("Higher in old", "Lower in old")
#' ORA（超几何，如 fgsea::fora、clusterProfiler::enricher）没有 NES：传 nes = 带方向的 log2 富集倍数（下调的取负），n_le = 重叠基因数，
#'   并设 effect_label = "log2 FE"、effect_title = "log2 fold enrichment"，轴标题和条外数字就不会写成 NES
enrich_pub <- function(d, v = c("fdr_bar_dot", "nes_bar_dot", "lollipop", "dot"), dir_labels = c("Up", "Down"),
                       effect_label = "NES", effect_title = "Normalized enrichment score (NES)", size_title = "Leading-edge\ngenes",
                       dir_cols = DIR_COLS[c("up", "down")], fdr_cols = c("#FFF5EB", "#FDAE6B", "#E6550D", "#7F2704"),
                       size_range = c(1.3, 3.2), wrap = 40, base_size = FIG_PT) {
  v <- match.arg(v)
  d <- .enrich_prep(d, wrap, dir_labels)
  dir_cols <- stats::setNames(unname(dir_cols), dir_labels); edge <- mix_col(dir_cols, "black", 0.3)
  # 点的大小：数据里最小的基因数 → 直径 size_range[1]，最大的 → size_range[2]，中间按面积线性插值（不像 scale_size_area 从 0 起算）
  size_sc <- ggplot2::scale_size(range = size_range, name = size_title, breaks = scales::breaks_pretty(3))
  th <- theme_pub(base_size) + ggplot2::theme(axis.ticks.y = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank(),
                                              axis.text.y = ggplot2::element_text(lineheight = 0.85), legend.key.height = ggplot2::unit(3, "mm"),
                                              legend.box = "vertical", legend.spacing.y = ggplot2::unit(1, "mm"))
  if (v == "fdr_bar_dot") {
    d$term_w <- stats::reorder(d$term_w, d$lfdr)
    return(ggplot2::ggplot(d, ggplot2::aes(lfdr, term_w)) +
      ggplot2::geom_col(ggplot2::aes(fill = dir, colour = dir), alpha = 0.35, width = 0.62, linewidth = 0.25) +
      ggplot2::geom_vline(xintercept = -log10(0.05), linewidth = 0.3, colour = "grey40", linetype = "dashed") +
      ggplot2::geom_point(ggplot2::aes(size = n_le, colour = dir), shape = 21, fill = "white", stroke = 0.45) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.1f", nes)), hjust = -0.55, size = (base_size - 1.5) / ggplot2::.pt, colour = "grey30") +
      ggplot2::scale_fill_manual(values = dir_cols, name = NULL) + ggplot2::scale_colour_manual(values = edge, guide = "none") + size_sc +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
      ggplot2::labs(x = bquote(-log[10]~FDR~"(number ="~.(effect_label)*")"), y = NULL) + th)
  }
  if (v == "nes_bar_dot") {
    d$term_w <- stats::reorder(d$term_w, d$nes)
    d$bar_fill <- grDevices::adjustcolor(dir_cols[as.character(d$dir)], alpha.f = 0.35); d$bar_edge <- edge[as.character(d$dir)]
    xr <- range(c(0, d$nes))
    arrows <- data.frame(x = xr, lab = c(paste("←", dir_labels[2]), paste(dir_labels[1], "→")), h = c(0, 1), col = dir_cols[c(2, 1)])
    arrows <- arrows[c(xr[1] < 0, xr[2] > 0), ]
    return(ggplot2::ggplot(d, ggplot2::aes(nes, term_w)) +
      ggplot2::geom_col(ggplot2::aes(fill = I(bar_fill), colour = I(bar_edge)), width = 0.62, linewidth = 0.25) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey30") +
      ggplot2::geom_point(ggplot2::aes(size = n_le, fill = lfdr), shape = 21, colour = "grey20", stroke = 0.25) +
      ggplot2::geom_text(data = arrows, ggplot2::aes(x = x, y = Inf, label = lab, hjust = h, colour = I(col)), inherit.aes = FALSE,
                         vjust = -0.6, size = (base_size - 1) / ggplot2::.pt) +
      ggplot2::scale_fill_gradientn(colours = fdr_cols, name = quote(-log[10]~FDR)) + size_sc +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.08)) + ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(x = effect_title, y = NULL) + th + ggplot2::theme(plot.margin = ggplot2::margin(5, 2, 2, 2, "mm")))
  }
  d$term_w <- stats::reorder(d$term_w, d$nes)
  if (v == "lollipop") {
    return(ggplot2::ggplot(d, ggplot2::aes(nes, term_w, colour = dir)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = nes, yend = term_w), linewidth = 0.5) +
      ggplot2::geom_point(ggplot2::aes(size = n_le), shape = 16) + ggplot2::geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) +
      ggplot2::scale_colour_manual(values = dir_cols, name = NULL) + size_sc +
      ggplot2::labs(x = effect_title, y = NULL) + th)
  }
  ggplot2::ggplot(d, ggplot2::aes(nes, term_w)) +                    # dot
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(size = n_le, fill = lfdr), shape = 21, colour = "grey30", stroke = 0.25) +
    ggplot2::scale_fill_gradientn(colours = fdr_cols, name = quote(-log[10]~FDR)) + size_sc +
    ggplot2::labs(x = effect_title, y = NULL) + th +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = 0.2))
}

#' 多组富集网格（类型 × 通路）。lim：NES 色标截断；caption 写明黑圈/黑框的含义
enrich_grid_pub <- function(d, style = c("dot", "tile", "triangle"), lim = 2.5, max_size = 3.2, ring_stroke = 0.35, edge_stroke = 0.15,
                            fdr_cut = 0.05, base_size = FIG_PT) {
  style <- match.arg(style)
  .need_cols(d, c("row", "col", "nes", "fdr"), "enrich_grid_pub")
  d$lfdr <- -log10(pmax(d$fdr, 1e-10)); d$sig <- !is.na(d$fdr) & d$fdr < fdr_cut
  th <- theme_pub(base_size) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), axis.title = ggplot2::element_blank(),
                                              panel.grid.major = ggplot2::element_line(colour = "grey93", linewidth = 0.2),
                                              legend.key.height = ggplot2::unit(3, "mm"),
                                              plot.caption = ggplot2::element_text(size = base_size - 1, colour = "grey30"))
  if (style == "dot") {
    return(ggplot2::ggplot(d, ggplot2::aes(col, row)) +
      ggplot2::geom_point(ggplot2::aes(size = lfdr, fill = nes), shape = 21, colour = "grey55", stroke = edge_stroke) +
      ggplot2::geom_point(data = d[d$sig, ], ggplot2::aes(size = lfdr), shape = 21, fill = NA, colour = "black", stroke = ring_stroke) +
      scale_div0(lim, aesthetics = "fill", name = "NES") + ggplot2::scale_size_area(max_size = max_size, name = quote(-log[10]~FDR)) +
      ggplot2::labs(caption = sprintf("Black ring: FDR < %s", fdr_cut)) + th)
  }
  if (style == "tile") {
    return(ggplot2::ggplot(d, ggplot2::aes(col, row)) +
      ggplot2::geom_tile(ggplot2::aes(fill = nes), colour = "white", linewidth = 0.4) +
      ggplot2::geom_tile(data = d[d$sig, ], fill = NA, colour = "black", linewidth = 0.4, width = 0.9, height = 0.9) +
      scale_div0(lim, aesthetics = "fill", name = "NES") + ggplot2::coord_fixed() +
      ggplot2::labs(caption = sprintf("Box: FDR < %s", fdr_cut)) + th +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()))
  }
  d$dir <- factor(ifelse(d$nes > 0, "Up", "Down"), c("Up", "Down"))
  ggplot2::ggplot(d, ggplot2::aes(col, row)) +
    ggplot2::geom_point(ggplot2::aes(shape = dir, fill = abs(nes), size = lfdr, colour = sig), stroke = 0.35) +
    ggplot2::scale_shape_manual(values = c(Up = 24, Down = 25), name = NULL) +
    ggplot2::scale_fill_gradientn(colours = c("#FEE6CE", "#FD8D3C", "#D94801", "#7F2704"), name = "|NES|") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75"), labels = c(`TRUE` = sprintf("FDR < %s", fdr_cut), `FALSE` = "n.s."), name = NULL) +
    ggplot2::scale_size_area(max_size = max_size, name = quote(-log[10]~FDR)) +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 2, fill = "grey60"))) + th
}
