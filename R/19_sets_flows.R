# =============================================================================
# 19_sets_flows.R —— 集合交集（UpSet）、类别对应（桑基图）
# upset_pub 输入 sets：命名列表（集合名 → 元素向量），集合自上而下按列表顺序
# sankey_pub 输入 d：left、right 两列类别 + weight 列（如细胞数）
# 用户选择（2026-10-12）：UpSet 用 ComplexUpset（CRAN）；桑基图两侧文字不能挤。
# 科研护栏：UpSet 画出所有 ≥ min_size 个元素的交集，不用"前 n 个"——前 n 个会把大小相同的交集从中间切开，挑哪几个是任意的
#   （2026-10-12 手画版和 ComplexUpset 版因此显示了不同的 4 个交集）。图注写明 min_size、交集是"互斥交集"。
# 来源（读原代码后重写）：alleninstitute/biccn_m1_evo（UpSetR + queries 着色）；drieslab/giotto_suite_manuscript（ggalluvial）
# =============================================================================

#' UpSet 图。engine = "complexupset"（默认）/ "ggplot"（不依赖包的手画版）；highlight：要着色的交集（集合名向量）
upset_pub <- function(sets, min_size = 4, highlight = NULL, hl_col = "#D6604D", engine = c("complexupset", "ggplot"), base_size = FIG_PT) {
  engine <- match.arg(engine)
  if (!is.list(sets) || is.null(names(sets))) stop("upset_pub：sets 要是命名列表", call. = FALSE)
  genes <- unique(unlist(sets)); memb <- sapply(sets, function(s) genes %in% s); rownames(memb) <- genes
  if (engine == "complexupset") {
    .need_pkg("ComplexUpset")
    m <- as.data.frame(memb)
    th <- ComplexUpset::upset_default_themes(text = ggplot2::element_text(size = base_size),
                                             axis.text = ggplot2::element_text(size = base_size, colour = "black"), axis.title = ggplot2::element_text(size = base_size))
    q <- if (is.null(highlight)) list() else list(ComplexUpset::upset_query(intersect = highlight, color = hl_col, fill = hl_col,
                                                                              only_components = c("intersections_matrix", "Intersection size")))
    return(ComplexUpset::upset(
      m, rev(names(sets)), name = NULL, min_size = min_size, sort_sets = FALSE, width_ratio = 0.22, height_ratio = 0.7,   # 倒序传入：第一个集合在最上面
      base_annotations = list(`Intersection size` = ComplexUpset::intersection_size(counts = TRUE, bar_number_threshold = 1,
                                                                                    text = list(size = (base_size - 1.5) / ggplot2::.pt), width = 0.65)),
      set_sizes = ComplexUpset::upset_set_size(geom = ggplot2::geom_bar(fill = "grey55", width = 0.6)) + ggplot2::ylab("Set size") +
        ggplot2::scale_y_reverse(breaks = scales::breaks_pretty(3)),
      matrix = ComplexUpset::intersection_matrix(geom = ggplot2::geom_point(size = 1.6), segment = ggplot2::geom_segment(linewidth = 0.5),
                                                 outline_color = list(active = "grey20", inactive = "grey88")),
      queries = q, themes = th, stripes = c("white", "grey96")))
  }
  key <- apply(memb, 1, function(v) paste(names(sets)[v], collapse = "&"))
  tab <- sort(table(key), decreasing = TRUE); tab <- tab[tab >= min_size]
  hl_key <- if (is.null(highlight)) "" else paste(names(sets)[names(sets) %in% highlight], collapse = "&")
  ix <- data.frame(id = factor(names(tab), levels = names(tab)), size = as.integer(tab)); ix$hl <- as.character(ix$id) == hl_key
  mat <- do.call(rbind, lapply(names(tab), function(k) data.frame(id = k, set = names(sets), on = names(sets) %in% strsplit(k, "&")[[1]])))
  mat$id <- factor(mat$id, levels = names(tab)); mat$set <- factor(mat$set, levels = rev(names(sets))); mat$hl <- mat$id %in% ix$id[ix$hl] & mat$on
  seg <- do.call(rbind, lapply(split(mat[mat$on, ], mat$id[mat$on]), function(z) if (nrow(z) > 1)
    data.frame(id = z$id[1], lo = min(as.integer(z$set)), hi = max(as.integer(z$set)), hl = any(z$hl))))
  bars <- ggplot2::ggplot(ix, ggplot2::aes(id, size, fill = hl)) + ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = size), vjust = -0.3, size = (base_size - 1.5) / ggplot2::.pt) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = hl_col, `FALSE` = "grey30"), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) + ggplot2::labs(x = NULL, y = "Intersection size") +
    theme_pub(base_size) + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  dots <- ggplot2::ggplot(mat, ggplot2::aes(id, set)) +
    ggplot2::geom_point(ggplot2::aes(colour = ifelse(hl, "hl", ifelse(on, "on", "off"))), size = 1.6) +
    (if (!is.null(seg)) ggplot2::geom_segment(data = seg, ggplot2::aes(x = id, xend = id, y = lo, yend = hi, colour = ifelse(hl, "hl", "on")), linewidth = 0.5)) +
    ggplot2::scale_colour_manual(values = c(hl = hl_col, on = "grey20", off = "grey88"), guide = "none") + ggplot2::labs(x = NULL, y = NULL) +
    theme_pub(base_size) + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank())
  side <- ggplot2::ggplot(data.frame(set = factor(names(sets), levels = rev(names(sets))), n = lengths(sets)), ggplot2::aes(n, set)) +
    ggplot2::geom_col(fill = "grey55", width = 0.6) + ggplot2::scale_x_reverse(breaks = scales::breaks_pretty(3)) +
    ggplot2::labs(x = "Set size", y = NULL) + theme_pub(base_size) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank())
  patchwork::wrap_plots(patchwork::plot_spacer(), bars, side, dots, ncol = 2, widths = c(1, 4), heights = c(2.2, 1))
}

#' 桑基图（ggalluvial）：左 = 注释 A，右 = 注释 B，流按左侧类别着色。
#' 标签放在层块外侧，先对准各自层块的中点；相邻标签太近时上下推开（一维推挤），只有被推开的标签才画一段短折线（引线只在必要时出现）。
#' panel_mm：面板大约多高（mm），用来把"一行字多高"换算成数据单位（图高减去上下边距）
sankey_pub <- function(d, left, right, weight, pal, left_title, right_title, panel_mm = NULL, w = 1 / 8, base_size = FIG_PT) {
  .need_pkg("ggalluvial")
  .need_cols(d, c(left, right, weight), "sankey_pub")
  d <- d[d[[weight]] > 0, ]
  d$L <- factor(d[[left]], levels = intersect(names(pal), unique(d[[left]]))); d$R <- d[[right]]; d$W <- d[[weight]]
  if (anyNA(d$L)) stop("sankey_pub：左侧有不在 pal 名字中的类别", call. = FALSE)
  p <- ggplot2::ggplot(d, ggplot2::aes(axis1 = L, axis2 = R, y = W)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = L), alpha = 0.55, width = w, knot.pos = 0.4, linewidth = 0) +
    ggalluvial::geom_stratum(width = w, fill = "grey97", colour = "grey40", linewidth = 0.2)
  st <- ggplot2::ggplot_build(p)$data[[2]]; st$yc <- (st$ymin + st$ymax) / 2
  lab_pt <- base_size - 1.5
  g <- sum(d$W) * (lab_pt * 0.3528 * 1.25) / (panel_mm %||% 135)   # 一行字的高度（数据单位）
  spread <- function(y, lo, hi) {
    o <- order(-y); q <- y[o]
    for (it in seq_len(500)) {
      moved <- FALSE
      for (i in seq_along(q)[-1]) { dd <- q[i - 1] - q[i]; if (dd < g) { s <- (g - dd) / 2; q[i - 1] <- q[i - 1] + s; q[i] <- q[i] - s; moved <- TRUE } }
      q <- pmin(pmax(q, lo), hi); if (!moved) break
    }
    out <- y; out[o] <- q; out
  }
  labs_ <- do.call(rbind, lapply(1:2, function(side) {
    z <- st[st$x == side, c("x", "yc", "stratum")]
    z$yl <- spread(z$yc, min(st$ymin) - 2 * g, max(st$ymax) + 2 * g)
    z$edge <- if (side == 1) 1 - w / 2 else 2 + w / 2; z$xt <- if (side == 1) z$edge - 0.07 else z$edge + 0.07
    z$hj <- if (side == 1) 1 else 0; z$moved <- abs(z$yl - z$yc) > 0.3 * g
    z
  }))
  mv <- labs_[labs_$moved, ]
  p + ggplot2::geom_segment(data = mv, ggplot2::aes(x = edge, xend = edge + ifelse(hj == 1, -0.025, 0.025), y = yc, yend = yc), inherit.aes = FALSE,
                            linewidth = 0.2, colour = "grey50") +
    ggplot2::geom_segment(data = mv, ggplot2::aes(x = edge + ifelse(hj == 1, -0.025, 0.025), xend = xt + ifelse(hj == 1, 0.008, -0.008), y = yc, yend = yl),
                          inherit.aes = FALSE, linewidth = 0.2, colour = "grey50") +
    ggplot2::geom_text(data = labs_, ggplot2::aes(x = xt, y = yl, label = stratum, hjust = hj), inherit.aes = FALSE, size = lab_pt / ggplot2::.pt) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_x_discrete(limits = c(left_title, right_title), expand = c(0.3, 0.3)) + ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) + theme_pub(base_size) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank())
}
