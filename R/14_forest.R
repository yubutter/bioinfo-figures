# =============================================================================
# 14_forest.R —— 森林图 / 效应量图（每行一个对象：效应 + 95% CI + 右侧表格）
# 输入 d：label（行名）、est、lo、hi（效应和 95% CI；未检验为 NA）、p（显示用的 FDR 或 P）、n_text（如 "9 / 10"）、
#         lineage（上色用的分组，如谱系）、status（未检验的原因，显示在该行；已检验写 ""）
# 用户选择（2026-10-12，第一批盲选）：Tabula Muris Senis 式的谱系色半透明横条 + 同色加深细边 + 带表格森林图的右侧小鼠数 / FDR 表格。
# 开关（缝合体）：style = "bar"（默认）/ "point"（点 + 区间）；table；stripes（GTEx 的灰白条纹）；order = "effect" / "given"
# 科研护栏：未检验的行保留并写原因（不静默删行）；区间来自按样本的模型。
# 来源（读原代码后重写）：czbiohub/tabula-muris-senis（Nature 2020）横条 + CI；miargentieri/proteomic-age-ukb（Nature Medicine 2024）右侧表格；
#   broadinstitute/gtex-v8（Science 2020，BSD-3）条纹
# =============================================================================

.forest_y <- function(d, order = c("effect", "given")) {
  order <- match.arg(order)
  lv <- if (order == "given") rev(as.character(d$label)) else {
    tested <- d[!is.na(d$est), ]; untested <- d[is.na(d$est), ]
    c(rev(as.character(untested$label)), as.character(tested$label[order(tested$est)]))   # 由下到上：未检验在最下，效应从小到大
  }
  d$label <- factor(d$label, levels = lv)
  d
}
.forest_na_layer <- function(d, base_size) {
  na <- d[is.na(d$est), , drop = FALSE]
  if (!nrow(na)) return(NULL)
  ggplot2::geom_text(data = na, ggplot2::aes(x = 0, y = label, label = status), inherit.aes = FALSE, hjust = 0,
                     size = (base_size - 1.5) / ggplot2::.pt, colour = "grey55", fontface = "italic")
}

#' 森林图（默认：谱系色横条 + 细边 + 黑色 CI + 右侧表格）。pal：lineage → 颜色；xlab 可用 lab_gene_lfc("Gapdh")
forest_pub <- function(d, pal, style = c("bar", "point"), table = TRUE, stripes = FALSE, order = c("effect", "given"),
                       xlab = "Effect (95% CI)", col_titles = c("n", "FDR"), sig_cut = 0.05, base_size = FIG_PT,
                       alpha = 0.55, outline_dark = 0.35, widths = c(3, 1.25), legend_pos = c(0.98, 0.02)) {
  style <- match.arg(style); order <- match.arg(order)
  .need_cols(d, c("label", "est", "lo", "hi", "p", "n_text", "lineage", "status"), "forest_pub")
  d <- .forest_y(d, order)
  n <- nlevels(d$label); ylim <- c(0.5, n + 0.5)
  xr <- range(c(d$lo, d$hi, 0), na.rm = TRUE); xlim <- xr + c(-1, 1) * 0.05 * diff(xr)
  d$lineage <- factor(ifelse(is.na(d$est), NA, d$lineage), levels = intersect(names(pal), unique(d$lineage[!is.na(d$est)])))
  d$sig <- !is.na(d$p) & d$p < sig_cut
  edge <- stats::setNames(mix_col(pal, "black", outline_dark), names(pal))
  p <- ggplot2::ggplot(d, ggplot2::aes(est, label))
  if (stripes) p <- p + ggplot2::geom_rect(data = data.frame(y = seq(2, n, by = 2)), ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
                                           inherit.aes = FALSE, fill = "grey94", colour = NA)
  p <- p + ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey35")
  if (style == "bar") {
    p <- p + ggplot2::geom_col(ggplot2::aes(fill = lineage, colour = lineage), alpha = alpha, width = 0.68, linewidth = 0.3, na.rm = TRUE) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi), width = 0.28, linewidth = 0.3, colour = "black", orientation = "y", na.rm = TRUE)
  } else {
    p <- p + ggplot2::geom_linerange(ggplot2::aes(xmin = lo, xmax = hi, colour = lineage), linewidth = 0.45, na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(fill = lineage), shape = 21, colour = "black", stroke = 0.3, size = 1.6, na.rm = TRUE)
  }
  p <- p + ggplot2::scale_colour_manual(values = edge, guide = "none", na.translate = FALSE) + .forest_na_layer(d, base_size) +
    ggplot2::scale_fill_manual(values = pal, name = NULL, drop = TRUE, na.translate = FALSE) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(x = xlab, y = NULL) + ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey85", linetype = "dashed", linewidth = 0.25),
                   panel.border = ggplot2::element_rect(colour = "grey40", linewidth = 0.35), axis.text = ggplot2::element_text(colour = "black"),
                   axis.ticks = ggplot2::element_line(linewidth = 0.25), axis.ticks.length.y = ggplot2::unit(0, "mm"),
                   legend.position = "inside", legend.position.inside = legend_pos, legend.justification = c(1, 0),
                   legend.background = ggplot2::element_rect(fill = "white", colour = "grey80", linewidth = 0.25),
                   legend.key.size = ggplot2::unit(2.5, "mm"), legend.margin = ggplot2::margin(1, 1.5, 1, 1.5, "mm"),
                   plot.margin = ggplot2::margin(5, 0, 2, 2, "mm"))
  if (!table) return(p)
  tab <- rbind(data.frame(label = d$label, col = 1, txt = d$n_text, sig = FALSE),
               data.frame(label = d$label, col = 2.25, txt = fmt_p_txt(d$p), sig = d$sig))
  right <- ggplot2::ggplot(tab, ggplot2::aes(col, label)) +
    ggplot2::geom_text(ggplot2::aes(label = txt, fontface = ifelse(sig, "bold", "plain")), size = (base_size - 0.5) / ggplot2::.pt, hjust = 0) +
    ggplot2::geom_text(data = data.frame(col = c(1, 2.25), txt = col_titles), ggplot2::aes(col, n + 0.9, label = txt), inherit.aes = FALSE,
                       fontface = "bold", size = (base_size - 0.5) / ggplot2::.pt, hjust = 0) +
    ggplot2::scale_x_continuous(limits = c(0.9, 3.2)) + ggplot2::scale_y_discrete(limits = levels(d$label)) +
    ggplot2::coord_cartesian(ylim = ylim, expand = FALSE, clip = "off") +   # 两边 y 范围完全相同，行才对得齐；表头画在面板上方
    ggplot2::theme_void(base_size = base_size) + ggplot2::theme(plot.margin = ggplot2::margin(5, 2, 2, 1, "mm"))
  patchwork::wrap_plots(p, right, widths = widths)
}
