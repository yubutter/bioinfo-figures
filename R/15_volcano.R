# =============================================================================
# 15_volcano.R —— 火山图
# 输入 d：gene、logFC、p（原始 P）、fdr。应是按样本合并后（pseudobulk）或混合模型的差异结果，不能是逐细胞检验。
# 统一判定：显著 = fdr < fdr_cut 且 |logFC| ≥ lfc_cut；y = −log10 P；水平阈值线画在"最大的显著 P"处；
#   标签 = 每个方向按 P 排序的前 n_lab 个显著基因（按规则选，不挑）；target 基因加黑圈（事先定好的目标基因，不能事后挑）。
# 用户选择（2026-10-12，第一批盲选）：style = "highlight"（浅灰底、柔和红蓝、标签与点同色加深且斜体、目标基因黑圈、封顶三角）；
#   引线 connect = "auto"：标签离点较远（放不下）时才画，学 cboix（用户："只有点密集、画不下的时候考虑连线"）。
# 开关（缝合体）：style = "outline"（silverbushlab：每点下垫半透明黑点作描边、NEJM 红蓝、灰虚线、斜体黑标签、四边框）、
#   "light"（cboix：grey90 底、Paired 浅红浅蓝、标签同色、x = 0 实线）；connect = "all" / "none"；cap（−log10 P 封顶）
# 科研护栏：min(p) < 1e-30 或行数像逐细胞检验时给出警告。
# 来源（读原代码后重写）：cboix/admultiregion_analysis（Nature 2024，无许可证）；silverbushlab do_VolcanoPlot（Cell 2026，MIT；由 SCpubr（GPL-3）改写）；
#   突出式 = 初筛 V1–V9 的组合（aifimmunology、sukeniklab、zhumy09 等）
# =============================================================================

.volcano_prep <- function(d, fdr_cut, lfc_cut, n_lab, cap = Inf) {
  .need_cols(d, c("gene", "logFC", "p", "fdr"), "volcano_pub")
  if (min(d$p, na.rm = TRUE) < 1e-30)
    warning("volcano_pub：最小 P < 1e-30。按样本（pseudobulk）检验很少这么小，检查是否按细胞检验（伪重复）", call. = FALSE)
  d <- d[!is.na(d$p) & !is.na(d$logFC), ]
  d$y <- -log10(pmax(d$p, 1e-300))
  d$dir <- factor(ifelse(d$fdr < fdr_cut & abs(d$logFC) >= lfc_cut, ifelse(d$logFC > 0, "up", "down"), "ns"), levels = c("ns", "down", "up"))
  d$capped <- d$y > cap; d$y_plot <- pmin(d$y, cap)
  d <- d[order(d$dir), ]                                            # 灰点先画，显著点在上层
  lab <- do.call(rbind, lapply(c("up", "down"), function(k) { z <- d[d$dir == k, ]; utils::head(z[order(z$p), ], n_lab) }))
  p_thr <- if (any(d$fdr < fdr_cut, na.rm = TRUE)) max(d$p[d$fdr < fdr_cut], na.rm = TRUE) else NA
  list(d = d, lab = lab, y_thr = -log10(p_thr))
}
.vx <- quote(log[2]~"fold change")
.vy <- quote(-log[10]~italic(P))

#' 火山图。target：要加黑圈的目标基因（可多个或 NULL）；xlab 可改成 quote(log[2]~"fold change (old vs young)")
volcano_pub <- function(d, target = NULL, style = c("highlight", "outline", "light"), connect = c("auto", "all", "none"),
                        fdr_cut = 0.05, lfc_cut = log2(1.5), n_lab = 10, cap = Inf, xlab = .vx, base_size = FIG_PT,
                        point_size = NULL, label_pt = 5.5, raster_dpi = 600, seed = 42L, cols = DIR_COLS) {
  style <- match.arg(style); connect <- match.arg(connect)
  min_seg <- switch(connect, auto = 0.5, all = 0, none = Inf)
  v <- .volcano_prep(d, fdr_cut, lfc_cut, n_lab, cap)
  if (style == "highlight") {
    ps <- point_size %||% 0.45
    lab <- v$lab[!v$lab$gene %in% target, , drop = FALSE]; tg <- v$d[v$d$gene %in% target, , drop = FALSE]
    lab_cols <- c(ns = "grey35", down = unname(mix_col(cols[["down"]], "black", 0.35)), up = unname(mix_col(cols[["up"]], "black", 0.35)))
    lab$lc <- lab_cols[as.character(lab$dir)]                     # 标签用点色的加深版，小字才看得清
    ylab <- if (any(v$d$capped)) bquote(-log[10]~italic(P)~"(capped at"~.(cap)*")") else .vy
    return(ggplot2::ggplot(v$d, ggplot2::aes(logFC, y_plot)) +
      ggplot2::geom_hline(yintercept = v$y_thr, colour = "grey55", linewidth = 0.25, linetype = "dashed", na.rm = TRUE) +
      ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), colour = "grey55", linewidth = 0.25, linetype = "dashed") +
      .rast(ggplot2::geom_point(ggplot2::aes(colour = dir, shape = capped, size = dir), stroke = 0), raster_dpi) +
      ggplot2::geom_point(data = tg, shape = 21, fill = NA, colour = "black", size = ps * 4, stroke = 0.45) +
      ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene, colour = I(lc)), size = label_pt / ggplot2::.pt, fontface = "italic",
                               segment.size = 0.2, min.segment.length = min_seg, box.padding = 0.2, max.overlaps = Inf, seed = seed, show.legend = FALSE) +
      ggrepel::geom_text_repel(data = tg, ggplot2::aes(label = gene), size = (label_pt + 0.5) / ggplot2::.pt, fontface = "bold.italic",
                               colour = "black", segment.size = 0.25, nudge_y = 0.6, box.padding = 0.4, seed = seed) +
      ggplot2::scale_colour_manual(values = cols, guide = "none") +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
      ggplot2::scale_size_manual(values = c(ns = ps, down = ps * 1.5, up = ps * 1.5), guide = "none") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
      ggplot2::labs(x = xlab, y = ylab) + theme_pub(base_size))
  }
  if (style == "outline") {
    ps <- point_size %||% 0.53
    set.seed(seed); dd <- v$d[sample(nrow(v$d)), ]               # 原函数把点打乱后再画
    mx <- max(abs(dd$logFC)) + 0.6
    return(ggplot2::ggplot(dd, ggplot2::aes(logFC, y_plot)) +
      .rast(ggplot2::geom_point(size = ps * 1.5, alpha = 0.5, colour = "black", shape = 16, stroke = 0), raster_dpi) +
      .rast(ggplot2::geom_point(ggplot2::aes(colour = dir), size = ps, alpha = 0.9, shape = 16, stroke = 0), raster_dpi) +
      ggplot2::geom_hline(yintercept = v$y_thr, colour = "grey75", linewidth = 0.25, linetype = "dashed", na.rm = TRUE) +
      ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), colour = "grey75", linewidth = 0.25, linetype = "dashed") +
      ggrepel::geom_text_repel(data = v$lab, ggplot2::aes(label = gene), colour = "black", fontface = "italic", size = label_pt / ggplot2::.pt,
                               max.overlaps = 1000, min.segment.length = min_seg, force = 20, segment.size = 0.2, seed = seed) +
      ggplot2::scale_colour_manual(values = c(up = "#BC3C29", down = "#0072B5", ns = "grey"), guide = "none") +
      ggplot2::coord_cartesian(xlim = c(-mx, mx)) + ggplot2::labs(x = xlab, y = .vy) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.5), axis.line = ggplot2::element_line(linewidth = 0),
                     axis.text = ggplot2::element_text(colour = "black"), axis.ticks = ggplot2::element_line(linewidth = 0.3)))
  }
  ps <- point_size %||% 0.38                                       # light（cboix）
  pal <- c(ns = "grey90", down = "#A6CEE3", up = "#FB9A99")         # RColorBrewer Paired 1、5 号
  ggplot2::ggplot(v$d, ggplot2::aes(logFC, y_plot, colour = dir)) +
    .rast(ggplot2::geom_point(size = ps, shape = 16, stroke = 0), raster_dpi) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggrepel::geom_text_repel(data = v$lab, ggplot2::aes(label = gene), size = 5 / ggplot2::.pt, max.overlaps = 100, segment.size = 0.2,
                             min.segment.length = min_seg, seed = seed, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") + ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = xlab, y = .vy) + theme_pub(base_size)
}
