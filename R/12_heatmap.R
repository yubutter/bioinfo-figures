# =============================================================================
# 12_heatmap.R —— 注释热图（ComplexHeatmap；卫冕方 v0.5 heatmap_pub；2026-10-12 盲选后加：行分组、显著格黑框开关）
# 输入 mat：行 = 基因、列 = 细胞类型/样本的数值矩阵
# 用户选择：第二批热图选"行分组"（不带黑框）；组名放右侧图例（"大块只写几个字、小块写不下"），左侧只留细色带。
# 来源：sscVis（z 截断 + colorRamp2）；行分组 mjgirgenti PTSD（Nature 2025）；黑框 allonkleinlab（读代码后重写）
# 依赖：ComplexHeatmap、circlize（不在默认库时把所在的库加进 .libPaths()）
# =============================================================================

#' 热图。z = TRUE 时按行 z 分数并截断到 ±clip；发散色以 0 为中心；格子之间细灰线；行名（基因）斜体。
#'   col_anno / anno_cols：列注释（列名 → 类别）及颜色，画成顶部色条；
#'   row_group / group_cols：行分组（基因 → 组）及颜色；group_labels = "legend"（默认：左侧细色带 + 右侧图例）或 "block"（组名写在色块里）；
#'   fdr + boxes = TRUE：fdr < 0.05 的格子加黑框（图例标题注明）；同形状的 fdr 矩阵（NA = 未检验）
heatmap_pub <- function(mat, col_anno = NULL, anno_cols = NULL, z = TRUE, clip = 2, base_size = FIG_PT, name = if (z) "z-score" else "value",
                        row_group = NULL, group_cols = NULL, group_labels = c("legend", "block"), group_title = NULL,
                        fdr = NULL, boxes = FALSE, fdr_cut = 0.05, italic_rows = TRUE, grid_col = "grey92", ...) {
  .need_pkg("ComplexHeatmap", bioc = TRUE); .need_pkg("circlize")
  group_labels <- match.arg(group_labels)
  m <- if (z) t(scale(t(mat))) else mat
  m[is.na(m)] <- 0
  col_fun <- circlize::colorRamp2(c(-clip, 0, clip), c("#2166AC", "white", "#B2182B"))
  gp <- grid::gpar(fontsize = base_size)
  lg <- function(title) list(title = title, title_gp = gp, labels_gp = gp, grid_width = grid::unit(2.5, "mm"), grid_height = grid::unit(2.5, "mm"))
  top <- NULL
  if (!is.null(col_anno)) {
    top <- ComplexHeatmap::HeatmapAnnotation(group = col_anno[colnames(m)], col = list(group = anno_cols[unique(col_anno[colnames(m)])]),
                                             show_annotation_name = FALSE, simple_anno_size = grid::unit(2, "mm"),
                                             annotation_legend_param = list(group = lg(NULL)))
  }
  la <- NULL; split <- NULL
  if (!is.null(row_group)) {
    split <- factor(row_group[rownames(m)], levels = intersect(names(group_cols) %||% unique(row_group), unique(row_group[rownames(m)])))
    gc <- group_cols %||% stats::setNames(cat_palette(nlevels(split)), levels(split))
    light <- stats::setNames(mix_col(gc[levels(split)], "white", 0.6), levels(split))
    la <- if (group_labels == "block") {
      ComplexHeatmap::rowAnnotation(grp = ComplexHeatmap::anno_block(gp = grid::gpar(fill = light, col = NA), labels = levels(split),
                                                                     labels_gp = grid::gpar(fontsize = base_size - 0.5, col = "grey15"),
                                                                     labels_rot = 0, width = grid::unit(16, "mm")))
    } else {
      ComplexHeatmap::rowAnnotation(grp = split, col = list(grp = light), show_annotation_name = FALSE, simple_anno_size = grid::unit(2, "mm"),
                                    annotation_legend_param = list(grp = lg(group_title)))
    }
  }
  cf <- NULL
  if (boxes) {
    if (is.null(fdr)) stop("heatmap_pub：boxes = TRUE 需要同形状的 fdr 矩阵", call. = FALSE)
    fdr <- fdr[rownames(m), colnames(m), drop = FALSE]
    cf <- function(j, i, x, y, w, h, fill) if (!is.na(fdr[i, j]) && fdr[i, j] < fdr_cut)
      grid::grid.rect(x, y, w * 0.88, h * 0.88, gp = grid::gpar(col = "black", fill = NA, lwd = 0.6))
  }
  ComplexHeatmap::Heatmap(m, name = name, col = col_fun, top_annotation = top, left_annotation = la,
                          row_split = split, cluster_row_slices = FALSE, row_title = if (is.null(split)) character(0) else NULL,
                          rect_gp = grid::gpar(col = grid_col, lwd = 0.3), cell_fun = cf,
                          row_names_gp = grid::gpar(fontsize = base_size - 0.5, fontface = if (italic_rows) "italic" else "plain"),
                          column_names_gp = gp, column_title_gp = gp,
                          heatmap_legend_param = list(title = if (boxes) sprintf("%s\n(box: FDR < %s)", name, fdr_cut) else name,
                                                      title_gp = gp, labels_gp = gp, at = c(-clip, 0, clip), grid_width = grid::unit(2, "mm")),
                          ...)
}
