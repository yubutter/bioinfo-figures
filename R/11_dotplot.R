# =============================================================================
# 11_dotplot.R —— 标志基因点图（卫冕方 v0.5 dotplot_pub；2026-10-12 盲选后加：基因分组括号、右侧细胞数条）
# 输入 d：每行一个"基因 × 细胞类型"，列名默认按 Seurat::DotPlot()$data（features.plot、id、pct.exp、avg.exp.scaled）
# 用户选择：第二批点图选"括号 + 细胞数条"（信息多；不要时去掉 = 不传 ncell）。
# 来源：Chu 2023（shape 21 + 发散填充）、KPMP；分组括号 HLCA / heiniglab（scanpy var_group_positions，读代码后重写）；
#   细胞数条 aifimmunology（读代码后重写）
# =============================================================================

#' 点图。圆 = 表达比例（scale_size_area，0% 时面积为 0）；填充 = 平均表达。
#'   scaled = TRUE：填充是组间标准化后的 z（有正负）→ 以 0 为中心的发散色；FALSE：非负平均表达 → 顺序色。
#'   gene_group：命名向量（基因 → 组），给了就按组分面，组名上方画一条横线当括号（组里 ≤ 3 个基因时组名斜放）；
#'   ncell：命名向量（类型 → 细胞数），给了就在右侧加一列细胞数横条（对数轴），提醒哪些类型细胞少、比例不可靠。
dotplot_pub <- function(d, gene = "features.plot", group = "id", pct = "pct.exp", value = "avg.exp.scaled",
                        scaled = TRUE, max_size = 3.2, lim = 2.5, min_pct = 1, base_size = FIG_PT,
                        gene_order = NULL, group_order = NULL, gene_group = NULL, group_levels = NULL, ncell = NULL, widths = c(8, 1)) {
  .need_cols(d, c(gene, group, pct, value), "dotplot_pub")
  dd <- data.frame(gene = d[[gene]], group = d[[group]], pct = d[[pct]], value = d[[value]])
  # 表达比例 < min_pct% 的格子不画点（否则满屏都是带描边的小圈）；坐标轴仍保留全部基因和类型
  dd$value[dd$pct < min_pct] <- NA; dd$pct[dd$pct < min_pct] <- NA
  dd$gene <- factor(dd$gene, levels = gene_order %||% unique(as.character(dd$gene)))
  dd$group <- factor(dd$group, levels = rev(group_order %||% unique(as.character(dd$group))))
  fill_scale <- if (scaled) scale_div0(lim, aesthetics = "fill", name = "Scaled\nexpression") else
    ggplot2::scale_fill_gradientn(colours = c("grey95", "#FDD49E", "#EF6548", "#990000"), name = "Mean\nexpression")
  p <- ggplot2::ggplot(dd, ggplot2::aes(gene, group)) +
    ggplot2::geom_point(ggplot2::aes(size = pct, fill = value), shape = 21, stroke = 0.2, colour = "grey30", na.rm = TRUE) +
    ggplot2::scale_size_area(max_size = max_size, name = "% expressing", breaks = c(25, 50, 75, 100), limits = c(0, 100)) +
    fill_scale + theme_pub(base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
                   axis.title = ggplot2::element_blank(), legend.position = "bottom", legend.box = "horizontal",
                   panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.2))
  if (!is.null(gene_group)) {
    gl <- group_levels %||% unique(unname(gene_group[levels(dd$gene)]))
    p$data$ggroup <- factor(gene_group[as.character(p$data$gene)], levels = gl)
    narrow <- any(table(gene_group[intersect(levels(dd$gene), names(gene_group))]) <= 3)
    st <- if (narrow) ggplot2::element_text(size = base_size - 0.5, angle = 35, hjust = 0, vjust = 0) else ggplot2::element_text(size = base_size - 0.5)
    p <- p + ggplot2::facet_grid(cols = ggplot2::vars(ggroup), scales = "free_x", space = "free_x") +
      ggplot2::geom_segment(data = data.frame(ggroup = factor(gl, levels = gl), x = -Inf, xend = Inf, y = Inf),
                            ggplot2::aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = st, panel.spacing.x = ggplot2::unit(1.2, "mm"), strip.clip = "off")
  }
  if (is.null(ncell)) return(p)
  tord <- levels(dd$group)
  nb <- data.frame(group = factor(tord, levels = tord), n = as.numeric(ncell[tord]), strip = "Cells")
  bar <- ggplot2::ggplot(nb, ggplot2::aes(n, group)) +
    ggplot2::geom_col(fill = "grey70", width = 0.6, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label = format(n, big.mark = ",", trim = TRUE)), hjust = -0.15, size = (base_size - 2) / ggplot2::.pt,
                       colour = "grey30", na.rm = TRUE) +
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0, 0.05)), breaks = c(10, 1000, 1e5),
                           labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(x = "Cells (log)", y = NULL) + theme_pub(base_size) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(2, 9, 2, 1, "mm"))
  # 右侧条形图也加一个同样样式的分面标题，上下高度才和点图对齐
  if (!is.null(gene_group)) bar <- bar + ggplot2::facet_grid(cols = ggplot2::vars(strip)) +
    ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text = p$theme$strip.text)
  patchwork::wrap_plots(p, bar, widths = widths) + patchwork::plot_layout(guides = "collect") & ggplot2::theme(legend.position = "bottom")
}
