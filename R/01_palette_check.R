# =============================================================================
# 01_palette_check.R —— 配色自检：UMAP 上相邻的类别对、色觉缺陷模拟与灰度下的色差（ΔE）、试色并排图
# 来源：single-cell-best-practice v0.5 S10（Palo 的思路，不装包；上色十条规则第 8、10 条）
# =============================================================================

#' UMAP 上相邻的类别对：把 UMAP 切成 n_bins × n_bins 的格子，两个类别在足够多的格子里同时出现（各自 ≥ min_cells 个细胞）就算相邻
umap_adjacent_pairs <- function(x, y, group, n_bins = 60L, min_cells = 3L, min_shared = 2L) {
  bx <- cut(x, n_bins, labels = FALSE); by <- cut(y, n_bins, labels = FALSE)
  tab <- table(paste(bx, by), as.character(group))
  pres <- tab >= min_cells
  co <- crossprod(pres * 1)
  idx <- which(upper.tri(co) & co >= min_shared, arr.ind = TRUE)
  if (!nrow(idx)) return(data.frame(a = character(), b = character(), shared_bins = integer()))
  out <- data.frame(a = rownames(co)[idx[, 1]], b = colnames(co)[idx[, 2]], shared_bins = co[idx])
  out[order(-out$shared_bins), , drop = FALSE]
}

#' 配色自检：正常视觉、三种色觉缺陷模拟（需要 colorspace）和灰度下，颜色两两之间的 CIELAB 距离（ΔE76）；ΔE < threshold 记为"易混"。
#' 阈值 12 是经验值【待验证】：ΔE 约 2 才刚能察觉，相邻小点要明显更大的差别。
palette_check <- function(pal, pairs = NULL, threshold = 12) {
  if (is.null(names(pal))) names(pal) <- pal
  views <- list(normal = pal)
  if (requireNamespace("colorspace", quietly = TRUE)) {
    views$deutan <- colorspace::deutan(pal); views$protan <- colorspace::protan(pal); views$tritan <- colorspace::tritan(pal)
  } else {
    message("colorspace 不可用，只检查正常视觉和灰度")
  }
  views$grey <- colorspace_free_grey(pal)
  if (is.null(pairs)) { cmb <- utils::combn(names(pal), 2); pairs <- data.frame(a = cmb[1, ], b = cmb[2, ]) }
  res <- do.call(rbind, lapply(names(views), function(v) {
    lab <- grDevices::convertColor(t(grDevices::col2rgb(views[[v]]) / 255), from = "sRGB", to = "Lab")
    rownames(lab) <- names(pal)
    d <- sqrt(rowSums((lab[pairs$a, , drop = FALSE] - lab[pairs$b, , drop = FALSE])^2))
    data.frame(view = v, a = pairs$a, b = pairs$b, delta_e = round(d, 1))
  }))
  res$confusable <- res$delta_e < threshold
  res[order(res$view, res$delta_e), ]
}

#' 灰度版本（按 CIELAB 亮度 L；不依赖 colorspace）
colorspace_free_grey <- function(pal) {
  lab <- grDevices::convertColor(t(grDevices::col2rgb(pal) / 255), from = "sRGB", to = "Lab")
  g <- grDevices::convertColor(cbind(lab[, 1], 0, 0), from = "Lab", to = "sRGB")
  stats::setNames(grDevices::rgb(pmin(pmax(g, 0), 1)), names(pal))
}

#' 试色：同一张 UMAP 用几套配色并排画出来，每格标出相邻类别对在绿色盲模拟下的最小 ΔE，给用户选
palette_audition <- function(df, group, palettes, x = "umap_1", y = "umap_2", n_max = 40000L, seed = 42L, base_size = FIG_PT) {
  set.seed(seed)
  sub <- if (nrow(df) > n_max) df[sample(nrow(df), n_max), , drop = FALSE] else df
  pairs <- umap_adjacent_pairs(df[[x]], df[[y]], df[[group]])
  plots <- lapply(names(palettes), function(nm) {
    chk <- palette_check(palettes[[nm]], pairs)
    worst <- chk[chk$view %in% c("normal", "deutan"), ]; worst <- worst[order(worst$delta_e), ][1, ]
    plot_umap_pub(sub, group, palettes[[nm]], x = x, y = y, outline = FALSE, labels = TRUE, arrows = FALSE, raster_dpi = 300,
                  base_size = base_size - 1,
                  title = sprintf("%s  (min adjacent dE %s: %.0f, %s-%s)", nm, worst$view, worst$delta_e, worst$a, worst$b))
  })
  patchwork::wrap_plots(plots, ncol = min(3L, length(plots)))
}
