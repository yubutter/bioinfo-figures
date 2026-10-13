# =============================================================================
# 17_composition.R —— 细胞组成（每个样本一根 100% 堆叠条）
# 输入 d：sample、group（组，如年龄；因子顺序 = 排列顺序）、comp（成分，如细胞类型；因子顺序 = 堆叠顺序）、count；
#         facet_extra 用到的列（如 dataset）可选
# 用户选择（2026-10-12，第二批盲选）：mode = "bar_n"——每个样本一根条、条顶写细胞数、下方组色带；颜色调淡 25%（用户："有些过于鲜艳"）。
# 开关（缝合体）：mode = "facet"（按组分面、浅灰标题条）/ "highlight"（只给一个成分上色；占比小时自动放大 y 轴，轴上是真实百分比）
# 科研护栏：每个样本一根条，不把一组的细胞合在一起画（否则看不出样本间差异）；条顶的 n 是细胞数，组间比较的单位仍是样本。
# 来源（读原代码后重写）：heiniglab 心脏图谱 stacked_barplot（条顶 n）、ding-lab 泛癌（分面）、ventolab（单成分）
# =============================================================================

.comp_prep <- function(d) {
  .need_cols(d, c("sample", "group", "comp", "count"), "comp_pub")
  tot <- stats::aggregate(count ~ sample, data = d, FUN = sum); d$total <- tot$count[match(d$sample, tot$sample)]
  d$frac <- d$count / d$total
  ord <- unique(d[order(d$group, d$sample), "sample"]); d$sample <- factor(d$sample, levels = ord)
  d
}

#' 组成图。pal：成分 → 颜色；gcols：组 → 颜色；soften：成分色和白色按此比例混合（面积大时调淡）；target：mode = "highlight" 时看的成分
comp_pub <- function(d, pal, gcols, mode = c("bar_n", "facet", "highlight"), soften = 0.25, target = NULL, facet_extra = NULL,
                     ylab = "Fraction of cells", zoom_below = 0.3, base_size = FIG_PT) {
  mode <- match.arg(mode)
  if (soften > 0) pal <- stats::setNames(mix_col(pal, "white", soften), names(pal))
  d <- .comp_prep(d)
  th <- theme_pub(base_size) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = base_size - 1.5),
                                              legend.key.size = ggplot2::unit(2.5, "mm"), legend.text = ggplot2::element_text(size = base_size - 1))
  strip_grey <- ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey90", colour = NA), strip.text = ggplot2::element_text(face = "plain"))
  if (mode == "bar_n") {
    tot <- unique(d[, c("sample", "group", "total")])
    return(ggplot2::ggplot(d, ggplot2::aes(sample, frac)) +
      ggplot2::geom_col(ggplot2::aes(fill = comp), width = 0.85, colour = "white", linewidth = 0.1) +
      ggplot2::geom_text(data = tot, ggplot2::aes(sample, 1.015, label = format(total, big.mark = ",", trim = TRUE)), angle = 60, hjust = 0,
                         size = (base_size - 2) / ggplot2::.pt, colour = "grey25") +
      ggplot2::geom_tile(data = tot, ggplot2::aes(sample, -0.035, colour = group), fill = NA, height = 0.03, width = 0.85, linewidth = 1.2) +
      ggplot2::scale_fill_manual(values = pal, name = NULL) + ggplot2::scale_colour_manual(values = gcols, name = NULL) +
      ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25), expand = ggplot2::expansion(mult = c(0, 0.12))) +
      ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(x = NULL, y = ylab) + th)
  }
  if (mode == "facet") {
    fml <- if (is.null(facet_extra)) ggplot2::vars(group) else ggplot2::vars(.data[[facet_extra]], group)
    return(ggplot2::ggplot(d, ggplot2::aes(sample, frac)) + ggplot2::geom_col(ggplot2::aes(fill = comp), width = 0.9) +
      ggplot2::facet_grid(cols = fml, scales = "free_x", space = "free_x") + ggplot2::scale_fill_manual(values = pal, name = NULL) +
      ggplot2::scale_y_continuous(expand = c(0, 0), labels = scales::percent) + ggplot2::labs(x = NULL, y = ylab) + th + strip_grey +
      ggplot2::theme(panel.spacing.x = ggplot2::unit(0.8, "mm"), axis.line.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()))
  }
  if (is.null(target) || !target %in% levels(d$comp)) stop("comp_pub：mode = 'highlight' 需要 target = 一个成分名", call. = FALSE)
  s <- stats::aggregate(frac ~ sample + group, data = transform(d, frac = ifelse(comp == target, frac, 0)), FUN = sum)
  zoom <- max(s$frac) < zoom_below      # 成分占比小时，100% 的条里几乎看不见：只画该成分、y 轴放大到实际范围
  long <- rbind(if (!zoom) data.frame(sample = s$sample, group = s$group, part = "Other", frac = 1 - s$frac),
                data.frame(sample = s$sample, group = s$group, part = target, frac = s$frac))
  long$part <- factor(long$part, levels = c("Other", target))
  ggplot2::ggplot(long, ggplot2::aes(sample, frac)) +
    ggplot2::geom_col(ggplot2::aes(fill = part, colour = part), width = 0.75, linewidth = 0.25) +
    ggplot2::geom_text(data = s, ggplot2::aes(sample, frac, label = sprintf("%.1f", 100 * frac)), vjust = -0.4,
                       size = (base_size - 2) / ggplot2::.pt, colour = unname(mix_col(pal[[target]], "black", 0.4))) +
    ggplot2::facet_grid(cols = ggplot2::vars(group), scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = c(Other = "grey92", stats::setNames(pal[[target]], target)), name = NULL, breaks = target) +
    ggplot2::scale_colour_manual(values = c(Other = NA, stats::setNames(unname(mix_col(pal[[target]], "black", 0.35)), target)), guide = "none") +
    ggplot2::scale_y_continuous(expand = c(0, 0), labels = scales::percent, limits = if (zoom) c(0, max(s$frac) * 1.18) else c(0, 1)) +
    ggplot2::labs(x = NULL, y = paste0(target, " (% of cells)")) + th + strip_grey
}
