# =============================================================================
# 10_umap.R —— UMAP / 嵌入图（卫冕方：v0.5 的样式；2026-10-12 新增突出子集、斜体基因小图阵列）
# 输入 df：每行一个细胞，至少有 umap_1、umap_2（可改 x/y）和分组列 / 基因表达列
# 来源：v0.5 S10（sscVis、KPMP、Humphreys 的点大小/透明度；scanpy add_outline 的描边比例）；
#   U1 子集突出：theislab、HLCA（读代码后重写）；U4 斜体基因小图阵列：tgoliver（读代码后重写）
# =============================================================================

#' 点大小和透明度（动态，2026-10-12）：按"每平方毫米多少个细胞"连续插值，细胞越密，点越小、越透明。
#'   标定点（89 mm 方形面板）：5,000 → 0.8 / 1；3 万 → 0.5 / 0.8；10 万 → 0.35 / 0.6；27.5 万 → 0.25 / 0.4（v0.5 在 27.5 万细胞上试画定的）。
#'   面板变大（如 Word 155 mm）同样的细胞数密度变低、点相应变大；超出标定范围时取端点值。新数据仍要看一眼，可用 tune_grid() 并排比几档。
#'   panel_mm：面板边长（mm，按方形算）；不是方形时传 sqrt(宽 × 高)
umap_point_style <- function(n, panel_mm = 89) {
  ref_n <- c(5e3, 3e4, 1e5, 2.75e5); ref_size <- c(0.8, 0.5, 0.35, 0.25); ref_alpha <- c(1, 0.8, 0.6, 0.4)
  x <- log10(ref_n / 89^2); xo <- log10(max(n, 1) / panel_mm^2)
  list(size = round(stats::approx(x, ref_size, xout = xo, rule = 2)$y, 3), alpha = round(stats::approx(x, ref_alpha, xout = xo, rule = 2)$y, 3))
}

#' 试参数：把一个作图函数在一组参数上各画一张、并排成一页，看了再挑（用户：动态调整的量"我看了图再调"）。
#'   make_plot：function(参数列表) → ggplot；grid：data.frame，每行一组参数；title_fmt 用列名拼标题
tune_grid <- function(make_plot, grid, ncol = NULL, base_size = FIG_PT) {
  ps <- lapply(seq_len(nrow(grid)), function(i) {
    pr <- as.list(grid[i, , drop = FALSE])
    make_plot(pr) + ggplot2::labs(title = paste(names(pr), unlist(pr), sep = " = ", collapse = ", ")) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = base_size, face = "plain"))
  })
  patchwork::wrap_plots(ps, ncol = ncol %||% ceiling(sqrt(length(ps))))
}

#' 按行数给出图高（mm）：森林图、点图、热图这类"每行一个对象"的图（每行 per_row mm + 上下留白 base mm）
fig_height_mm <- function(n_rows, per_row = 4.2, base = 22, max_mm = 230) min(max_mm, base + per_row * n_rows)

#' 左下角的小箭头坐标轴（代替完整坐标轴）
umap_arrow_layers <- function(x, y, base_size = FIG_PT, frac = 0.12) {
  rx <- range(x); ry <- range(y)
  x0 <- rx[1] - 0.02 * diff(rx); y0 <- ry[1] - 0.02 * diff(ry)
  arr <- data.frame(x = x0, y = y0, xend = c(x0 + frac * diff(rx), x0), yend = c(y0, y0 + frac * diff(ry)))
  list(
    ggplot2::geom_segment(data = arr, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE,
                          arrow = ggplot2::arrow(length = ggplot2::unit(1.2, "mm"), type = "closed"), linewidth = 0.3, colour = "grey25"),
    ggplot2::annotate("text", x = x0, y = y0, label = "UMAP 1", hjust = 0, vjust = 1.5, size = (base_size - 1) / ggplot2::.pt, colour = "grey25"),
    ggplot2::annotate("text", x = x0, y = y0, label = "UMAP 2", hjust = 0, vjust = -0.5, angle = 90, size = (base_size - 1) / ggplot2::.pt, colour = "grey25"))
}

#' 细胞类型 UMAP：打乱绘制顺序；点大小/透明度默认取 umap_point_style()；只栅格化点；
#'   outline = TRUE 时加簇外描边（只画在细胞密集的格子里，默认关：簇松散时会出碎线）；labels = TRUE 时在各类中位点写类名（白色描边），不放图例
plot_umap_pub <- function(df, group, pal, x = "umap_1", y = "umap_2", size = NULL, alpha = NULL,
                          outline = FALSE, outline_min_bin = 4L, labels = TRUE, legend = !labels, arrows = TRUE,
                          raster_dpi = 600, base_size = FIG_PT, seed = 42L, label_size = base_size, title = NULL, panel_mm = 89) {
  .need_cols(df, c(x, y, group), "plot_umap_pub")
  sty <- umap_point_style(nrow(df), panel_mm)
  size <- size %||% sty$size; alpha <- alpha %||% sty$alpha
  set.seed(seed); df <- df[sample(nrow(df)), , drop = FALSE]
  d <- data.frame(x = df[[x]], y = df[[y]], g = factor(df[[group]], levels = names(pal)))
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y))
  if (outline) {
    key <- paste(cut(d$x, 150, labels = FALSE), cut(d$y, 150, labels = FALSE))
    dense <- d[as.vector(table(key)[key]) >= outline_min_bin, , drop = FALSE]
    p <- p + .rast(ggplot2::geom_point(data = dense, colour = "black", size = size * 1.7, shape = 16, stroke = 0), raster_dpi) +
      .rast(ggplot2::geom_point(data = dense, colour = "white", size = size * 1.1, shape = 16, stroke = 0), raster_dpi)
  }
  p <- p + .rast(ggplot2::geom_point(ggplot2::aes(colour = g), size = size, alpha = alpha, shape = 16, stroke = 0), raster_dpi) +
    ggplot2::scale_colour_manual(values = pal, name = NULL, drop = TRUE) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
    ggplot2::coord_fixed() + theme_umap_pub(base_size) + ggplot2::labs(title = title)
  if (labels) {
    lab <- stats::aggregate(cbind(x, y) ~ g, data = d, FUN = stats::median)
    p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(x, y, label = g), inherit.aes = FALSE,
                                      size = label_size / ggplot2::.pt, colour = "grey10", bg.color = "white", bg.r = 0.15,
                                      segment.size = 0.2, segment.colour = "grey30", min.segment.length = 0.2,
                                      box.padding = 0.25, max.overlaps = Inf, seed = seed)
  }
  if (arrows) p <- p + umap_arrow_layers(d$x, d$y, base_size)
  if (!legend) p <- p + ggplot2::theme(legend.position = "none")
  p
}

#' 单基因 UMAP：0 值浅灰垫底，阳性细胞按表达量从低到高画在上层并略放大；顺序色（浅橙 → 深红），各面板共用 limits。
#' 图注必须写"阳性细胞画在上层并略放大"。
plot_feature_pub <- function(df, value, x = "umap_1", y = "umap_2", size = NULL, enlarge = 1.6, limits = NULL,
                             zero_col = "#E4E4E4", cols = c("#FDBB84", "#EF6548", "#B30000", "#67000D"),
                             facet = NULL, raster_dpi = 600, base_size = FIG_PT, name = value, panel_mm = 89) {
  .need_cols(df, c(x, y, value), "plot_feature_pub")
  size <- size %||% umap_point_style(nrow(df), panel_mm)$size
  d <- data.frame(x = df[[x]], y = df[[y]], v = df[[value]])
  if (!is.null(facet)) d$facet <- df[[facet]]
  zero <- d[d$v <= 0, , drop = FALSE]; pos <- d[d$v > 0, , drop = FALSE]; pos <- pos[order(pos$v), , drop = FALSE]
  limits <- limits %||% if (nrow(pos)) range(pos$v) else c(0, 1)   # 全为 0 时给占位范围，避免 range() 返回 Inf
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    .rast(ggplot2::geom_point(data = zero, colour = zero_col, size = size, shape = 16, stroke = 0), raster_dpi) +
    .rast(ggplot2::geom_point(data = pos, ggplot2::aes(colour = v), size = size * enlarge, shape = 16, stroke = 0), raster_dpi) +
    ggplot2::scale_colour_gradientn(colours = cols, limits = limits, oob = scales::squish, name = name) +
    ggplot2::coord_fixed() + theme_umap_pub(base_size) +
    ggplot2::theme(legend.key.width = ggplot2::unit(2, "mm"), legend.key.height = ggplot2::unit(4, "mm"))
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(~facet)
  p
}

#' 只突出一个子集（某个数据集、某个细胞类型、阳性细胞），其余浅灰打底；子集点略大、画在上层；副标题写子集细胞数
umap_highlight_pub <- function(df, flag, col = "#D55E00", title = NULL, x = "umap_1", y = "umap_2", base_size = FIG_PT,
                               raster_dpi = 600, seed = 42L, arrows = TRUE, panel_mm = 89) {
  .need_cols(df, c(x, y, flag), "umap_highlight_pub")
  sty <- umap_point_style(nrow(df), panel_mm)
  d <- data.frame(umap_1 = df[[x]], umap_2 = df[[y]], f = as.logical(df[[flag]]))
  bg <- d[!d$f, ]; fg <- d[d$f, ]; set.seed(seed); fg <- fg[sample(nrow(fg)), ]
  p <- ggplot2::ggplot(d, ggplot2::aes(umap_1, umap_2)) +
    .rast(ggplot2::geom_point(data = bg, colour = "grey88", size = sty$size, shape = 16, stroke = 0), raster_dpi) +
    .rast(ggplot2::geom_point(data = fg, colour = col, size = sty$size * 1.3, shape = 16, stroke = 0, alpha = min(1, sty$alpha + 0.3)), raster_dpi) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, subtitle = sprintf("%s of %s cells", format(nrow(fg), big.mark = ","), format(nrow(d), big.mark = ","))) +
    theme_umap_pub(base_size) + ggplot2::theme(plot.subtitle = ggplot2::element_text(size = base_size - 1, colour = "grey35"))
  if (arrows) p <- p + umap_arrow_layers(d$umap_1, d$umap_2, base_size)
  p
}

#' 多个基因的表达 UMAP 小图阵列：基因名斜体作小标题，每个小图一根短色标（各基因色标独立，图注要说明）
umap_genes_pub <- function(df, genes, ncol = 2, x = "umap_1", y = "umap_2", base_size = FIG_PT) {
  .need_cols(df, c(x, y, genes), "umap_genes_pub")
  ps <- lapply(genes, function(g) plot_feature_pub(df, g, x = x, y = y, base_size = base_size, name = NULL) +
                 ggplot2::labs(title = g) +
                 ggplot2::theme(plot.title = ggplot2::element_text(face = "italic", size = base_size + 0.5),
                                legend.key.height = ggplot2::unit(3, "mm"), legend.key.width = ggplot2::unit(1.5, "mm")))
  patchwork::wrap_plots(ps, ncol = ncol)
}

#' 图注里写明 UMAP 参数（UMAP 只作展示；距离、簇大小、簇间位置不可解读）
umap_caption <- function(reduction, dims, n_neighbors, min_dist, seed, metric = "cosine") {
  sprintf("UMAP computed on %s dimensions %s (n.neighbors = %s, min.dist = %s, metric = %s, seed = %s); for display only: distances, cluster sizes and relative positions are not interpretable.",
          reduction, dims, n_neighbors, min_dist, metric, seed)
}
