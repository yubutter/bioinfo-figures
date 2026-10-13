# =============================================================================
# test_templates.R —— bioinfo-figures 回归测试：每个模板在模拟的"常规"和"难"数据上都能画，且科研护栏生效
# 运行（任意位置）：Rscript --vanilla <仓库目录>/tests/test_templates.R
# 包库：设了环境变量 BIOFIG_LIBS（分号分隔）就用它；否则仓库根目录有 local_libpaths.R（本机专用，不进仓库）就执行它；都没有就用默认 .libPaths()。
# 改了任何代码后都要跑，必须 ALL PASS。
# =============================================================================
skill <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  d <- if (length(f)) dirname(dirname(normalizePath(f[1], mustWork = FALSE))) else Sys.getenv("BIOFIG_DIR", getwd())
  if (!file.exists(file.path(d, "load.R"))) d <- Sys.getenv("BIOFIG_DIR", getwd())
  d
})
libs <- Sys.getenv("BIOFIG_LIBS", "")
if (nzchar(libs)) { libs <- strsplit(libs, ";")[[1]]; .libPaths(libs[dir.exists(libs)])
} else if (file.exists(file.path(skill, "local_libpaths.R"))) source(file.path(skill, "local_libpaths.R"), encoding = "UTF-8")
source(file.path(skill, "load.R"), encoding = "UTF-8")
source(file.path(skill, "examples", "demo_data.R"), encoding = "UTF-8")
cat("bioinfo-figures", FIG_SKILL_VERSION, "\n")
tmp <- file.path(tempdir(), "biofig_test"); dir.create(tmp, showWarnings = FALSE)
checks <- data.frame(check = character(), passed = logical(), detail = character())
add <- function(name, ok, detail = "") checks[nrow(checks) + 1, ] <<- list(name, isTRUE(ok), detail)
try_ok <- function(expr) !inherits(tryCatch(expr, error = function(e) e), "error")
layer_data_of <- function(p, geom) {           # 某类 geom 的作图数据（合并所有同类图层）
  b <- ggplot2::ggplot_build(p); i <- which(vapply(p$layers, function(l) inherits(l$geom, geom), TRUE))
  if (!length(i)) return(NULL); do.call(rbind, lapply(i, function(k) b$data[[k]][, intersect(c("x", "y", "label", "group", "PANEL"), names(b$data[[k]])), drop = FALSE]))
}
save_ok <- function(p, stem, w = 89, h = 70) try_ok(save_plot_pub(p, tmp, stem, w, h, dpi = 150)) && file.exists(file.path(tmp, paste0(stem, ".pdf")))

# ---- 核心 ----------------------------------------------------------------------------------------------------
e <- lab_italic("Gapdh"); add("lab_italic_quoted_string", is.character(e[[2]]), "italic() 的参数必须是字符串，否则字母会被排开")
add("fmt_p_math", identical(fmt_p_math(0.0085, "FDR"), "FDR == '0.0085'") && grepl("%\\*% 10\\^-4", fmt_p_math(2.4e-4)))
add("fmt_p_txt", identical(fmt_p_txt(c(0.0085, 2.4e-4)), c("0.0085", "2.4e-4")))
add("clean_term", identical(clean_term(c("HALLMARK_MTORC1_SIGNALING", "REACTOME_RRNA_PROCESSING")), c("MTORC1 signaling", "rRNA processing")))
lp <- lineage_palette(DEMO_LINEAGE, DEMO_BASE); add("lineage_palette", length(lp) == length(DEMO_LINEAGE) && identical(names(lp), names(DEMO_LINEAGE)))
add("need_cols_stops", inherits(tryCatch(forest_pub(data.frame(a = 1), DEMO_BASE), error = function(e) e), "error"))
p0 <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
save_plot_pub(p0, tmp, "with_data", 40, 30, dpi = 100, data = data.frame(x = 1:3, note = c("a", "b", "\u5e74\u8f7b")))
add("save_plot_pub_data_csv", file.exists(file.path(tmp, "with_data_data.csv")) &&
      identical(utils::read.csv(file.path(tmp, "with_data_data.csv"), encoding = "UTF-8")$note[3], "\u5e74\u8f7b"))

pv <- fig_skill_provenance()
has_git <- dir.exists(file.path(pv$path, ".git"))   # 从 GitHub 下载 ZIP 时没有 .git，提交号只能是 NA
add("provenance_version_commit", !is.na(pv$version) && (if (has_git) grepl("^[0-9a-f]{40}$", pv$commit) else is.na(pv$commit)),
    paste(pv$version, pv$commit, if (!has_git) "（无 .git）"))
fz <- file.path(tmp, "frozen"); unlink(fz, recursive = TRUE); fig_skill_freeze(fz)
add("freeze_copies_and_refuses_overwrite", file.exists(file.path(fz, "load.R")) && file.exists(file.path(fz, "R", "14_forest.R")) &&
      inherits(tryCatch(fig_skill_freeze(fz), error = function(e) e), "error"))

# ---- 动态参数 -------------------------------------------------------------------------------------------------
s1 <- umap_point_style(2.75e5); s2 <- umap_point_style(5e3); s3 <- umap_point_style(2.75e5, panel_mm = 155); s4 <- umap_point_style(1.5e5)
add("point_style_calibration_points", s1$size == 0.25 && s1$alpha == 0.4 && s2$size == 0.8 && s2$alpha == 1, "89 mm 面板上与 v0.5 标定一致")
add("point_style_continuous_and_panel_aware", s4$size < 0.35 && s4$size > 0.25 && s3$size > s1$size, "两档之间连续；面板变大点变大")
tg <- tune_grid(function(pr) ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) + ggplot2::geom_point(size = pr$size, alpha = pr$alpha),
                expand.grid(size = c(0.2, 0.4), alpha = c(0.4, 0.8)))
add("tune_grid_builds", inherits(tg, "patchwork") && length(tg$patches$plots) + 1 == 4)
add("fig_height_mm", fig_height_mm(10) == 64 && fig_height_mm(100) == 230)

# ---- UMAP ----------------------------------------------------------------------------------------------------
um <- demo_umap(); pal_u <- lineage_palette(DEMO_LINEAGE, DEMO_BASE)
add("plot_umap_pub", save_ok(plot_umap_pub(um, "celltype", pal_u), "umap"))
add("plot_feature_pub_all_zero", save_ok(plot_feature_pub(transform(um, GeneZ = 0), "GeneZ"), "feat0"))
um$is_d3 <- um$dataset == "D3"
add("umap_highlight_pub", save_ok(umap_highlight_pub(um, "is_d3", title = "Dataset D3"), "u1"))
add("umap_genes_pub", save_ok(umap_genes_pub(um, c("GeneA", "GeneB", "GeneC", "GeneD")), "u4", 120, 110))

# ---- 点图 ----------------------------------------------------------------------------------------------------
dp <- data.frame(features.plot = rep(c("g1", "g2"), each = 3), id = rep(c("a", "b", "c"), 2), pct.exp = c(80, 0.5, 20, 5, 90, 40),
                 avg.exp.scaled = c(4, -1, 0, -1, 1, 0))
b_dp <- ggplot2::ggplot_build(dotplot_pub(dp))$data[[1]]
add("dotplot_pub_compat_min_pct", sum(!is.na(b_dp$size)) == 5, "v0.5 旧调用：表达比例 < 1% 的格子不画点")
for (hard in c(FALSE, TRUE)) {
  z <- demo_dot(hard)
  p <- dotplot_pub(z$d, gene = "gene", group = "celltype", pct = "pct", value = "z", gene_order = z$genes, group_order = z$types,
                   gene_group = z$gene_group, ncell = z$ncell)
  add(paste0("dotplot_brackets_n_", if (hard) "hard" else "normal"), inherits(p, "patchwork") && save_ok(p, paste0("dot_", hard), 140, 80))
}

# ---- 热图 ----------------------------------------------------------------------------------------------------
if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  hz <- demo_heat()
  ht <- heatmap_pub(hz$m, col_anno = hz$col_anno, anno_cols = DEMO_BASE, z = FALSE, clip = 3, name = "log2FC",
                    row_group = hz$row_group, group_cols = DEMO_BASE, fdr = hz$fdr, boxes = TRUE, cluster_columns = FALSE)
  add("heatmap_rowgroup_boxes", inherits(ht, "Heatmap") && save_ok(ht, "heat", 110, 110))
  add("heatmap_boxes_need_fdr", inherits(tryCatch(heatmap_pub(hz$m, boxes = TRUE), error = function(e) e), "error"))
  add("heatmap_compat_v05", save_ok(heatmap_pub(hz$m[1:5, 1:5], col_anno = hz$col_anno, anno_cols = DEMO_BASE), "heat05", 60, 40))
} else add("heatmap", TRUE, "ComplexHeatmap 不在库中，跳过")

# ---- 箱线 / 小提琴（护栏）-------------------------------------------------------------------------------------
bn <- demo_box(); bh <- demo_box(hard = TRUE)
pb <- box_pub(bn$d, bn$s, bn$cols, "value"); txt <- layer_data_of(pb, "GeomText")
add("box_normal_p_and_n", sum(grepl("FDR", txt$label)) == 3 && sum(grepl("^n = ", txt$label)) == 6)
pbh <- box_pub(bh$d, bh$s, bh$cols, "value"); txth <- layer_data_of(pbh, "GeomText")
add("box_hard_no_p_when_n_lt_3", sum(grepl("FDR", txth$label)) == 1, "Type A1 只有 3 个年轻样本（刚够），Type C1 年轻组 n = 2 不显示 p")
bx <- ggplot2::ggplot_build(pbh); box_layers <- which(vapply(pbh$layers, function(l) inherits(l$geom, "GeomBoxplot"), TRUE))
n_boxes <- nrow(bx$data[[box_layers[1]]])
add("box_hard_no_box_when_n_lt_3", n_boxes == 5, sprintf("应画 5 个箱（A1×3、C1 的 Old/Very old），实际 %d", n_boxes))
pv <- violin_pub(bh$d, bh$s, bh$cols, "value"); vl <- which(vapply(pv$layers, function(l) inherits(l$geom, "GeomViolin"), TRUE))
n_vio <- length(unique(ggplot2::ggplot_build(pv)$data[[vl]]$group * 10 + as.integer(ggplot2::ggplot_build(pv)$data[[vl]]$PANEL)))
add("violin_hard_no_violin_when_n_lt_5", n_vio == 1, sprintf("只有 C1 的 Old（n = 5）应画小提琴，实际 %d 个", n_vio))
add("violin_framed_style", save_ok(violin_pub(bn$d, bn$s, bn$cols, "value", style = "framed"), "vio_framed", 120, 58))
add("box_points_pub_compat", try_ok(ggplot2::ggplot_build(box_points_pub(data.frame(age = rep(c("young", "old"), each = 3), rate = 1:6),
                                                                      "age", "rate", c(young = "#56B4E9", old = "#D55E00")))))

# ---- 森林图 ---------------------------------------------------------------------------------------------------
fh <- demo_forest(hard = TRUE)
pf <- forest_pub(fh, DEMO_BASE, xlab = lab_gene_lfc("Target1"))
left <- pf[[1]]; tx <- layer_data_of(left, "GeomText")
add("forest_keeps_untested_rows", nlevels(left$data$label) == nrow(fh) && all(c("too few samples", "gene too low to test") %in% tx$label))
add("forest_table_and_point_style", save_ok(pf, "forest", 120, 60) && save_ok(forest_pub(fh, DEMO_BASE, style = "point", stripes = TRUE, table = FALSE), "forest_pt"))

# ---- 火山图 ---------------------------------------------------------------------------------------------------
vd <- demo_volcano()
for (st in c("highlight", "outline", "light")) add(paste0("volcano_", st), save_ok(volcano_pub(vd, target = "Target1", style = st), paste0("vol_", st)))
vn <- demo_volcano(none_sig = TRUE)
add("volcano_no_significant_genes", save_ok(volcano_pub(vn, target = "Target1"), "vol_none"))
vt <- vd; vt$p[1] <- 1e-60
add("volcano_warns_tiny_p", inherits(tryCatch(volcano_pub(vt), warning = function(w) w), "warning"), "最小 P < 1e-30 时提醒检查逐细胞检验")

# ---- 富集 ----------------------------------------------------------------------------------------------------
for (v in c("fdr_bar_dot", "nes_bar_dot", "lollipop", "dot")) add(paste0("enrich_", v), save_ok(enrich_pub(demo_enrich(), v), paste0("en_", v), 110, 60))
add("enrich_one_direction_long_names", save_ok(enrich_pub(demo_enrich(25, one_direction = TRUE, long_names = TRUE), "nes_bar_dot"), "en_hard", 120, 130))
ora <- enrich_pub(demo_enrich(), effect_label = "log2 FE", effect_title = "log2 fold enrichment", size_title = "Overlap genes")
add("enrich_ora_labels", grepl("log2 FE", paste(deparse(ora$labels$x), collapse = "")) && identical(ora$scales$get_scales("size")$name, "Overlap genes"),
    "ORA 没有 NES：轴标题和点大小图例可改")
sz <- ggplot2::ggplot_build(enrich_pub(demo_enrich()))$data[[3]]$size
add("enrich_point_size_floor_cap", abs(min(sz) - 1.3) < 1e-6 && abs(max(sz) - 3.2) < 1e-6, "最小点 = 1.3、最大点 = 3.2（用户：最小的不能过小）")
for (st in c("dot", "tile", "triangle")) add(paste0("grid_", st), save_ok(enrich_grid_pub(demo_grid(30), st), paste0("grid_", st), 110, 120))
rl <- enrich_grid_pub(demo_grid())$layers[[2]]$aes_params$stroke
add("grid_ring_thin", isTRUE(all.equal(rl, 0.35)), "黑圈 0.35（用户：边框稍微细一点）")

# ---- 组成 ----------------------------------------------------------------------------------------------------
cn <- demo_comp(); ch <- demo_comp(hard = TRUE)
for (m in c("bar_n", "facet")) add(paste0("comp_", m), save_ok(comp_pub(ch$d, ch$pal, ch$gcols, mode = m), paste0("comp_", m), 183, 85))
chd <- ch$d; chd$dataset <- ifelse(as.integer(sub("S", "", chd$sample)) <= 15, "D1", "D2")
add("comp_facet_extra_dataset", save_ok(comp_pub(chd, ch$pal, ch$gcols, mode = "facet", facet_extra = "dataset"), "comp_facet2", 183, 85))
ph <- comp_pub(cn$d, cn$pal, cn$gcols, mode = "highlight", target = "Lineage E")
add("comp_highlight_zoom_small_fraction", max(ggplot2::ggplot_build(ph)$layout$panel_params[[1]]$y.range) < 0.5, "占比小时 y 轴放大")
add("comp_soften_default", !identical(ggplot2::ggplot_build(comp_pub(cn$d, cn$pal, cn$gcols))$data[[1]]$fill[1], unname(cn$pal[1])))

# ---- 相关散点 -------------------------------------------------------------------------------------------------
sd_ <- demo_scatter()
for (st in c("black", "group", "red")) add(paste0("scatter_", st), save_ok(scatter_pub(sd_, "x", "y", DEMO_GROUP_COLS, style = st), paste0("sc_", st)))
add("scatter_n_in_label", grepl("italic(n)", paste(deparse(scatter_pub(sd_, "x", "y", DEMO_GROUP_COLS)$labels$subtitle), collapse = ""), fixed = TRUE))

# ---- UpSet / 桑基 ---------------------------------------------------------------------------------------------
sets <- demo_sets()
pg <- upset_pub(sets, min_size = 4, engine = "ggplot", highlight = c("A", "B"))
n_bars <- nrow(pg[[2]]$data)
add("upset_ggplot_keeps_ties", n_bars == 9, sprintf("≥ 4 的互斥交集应为 9 个（含 3 个并列为 4 的），实际 %d", n_bars))
if (requireNamespace("ComplexUpset", quietly = TRUE)) {
  add("upset_complexupset_draws", save_ok(upset_pub(sets, min_size = 4, highlight = c("A", "B")), "upset_cu", 120, 80))
} else add("upset_complexupset_draws", TRUE, "ComplexUpset 不在库中，跳过")
if (requireNamespace("ggalluvial", quietly = TRUE)) {
  ps <- sankey_pub(demo_sankey(), "left", "right", "cells", lineage_palette(DEMO_LINEAGE, DEMO_BASE), "Ours", "Author")
  tl <- layer_data_of(ps, "GeomText")
  gap_ok <- all(vapply(split(tl, round(tl$x)), function(z) { y <- sort(z$y); length(y) < 2 || min(diff(y)) > 0 }, TRUE))
  add("sankey_labels_spread", gap_ok && save_ok(ps, "sankey", 150, 150), "同一侧的标签上下不重合")
} else add("sankey_labels_spread", TRUE, "ggalluvial 不在库中，跳过")

# ---- GSEA / PCA -----------------------------------------------------------------------------------------------
gz <- demo_gsea()
add("gsea_curve_pub", save_ok(gsea_curve_pub(gz$ranks, gz$set, "Demo set", 1.9, 3e-4), "gsea"))
pz <- demo_pca()
add("pca_pub", save_ok(pca_pub(pz$logcpm, pz$info, "Group", "Dataset", c(Young = "#56B4E9", Old = "#D55E00")), "pca"))

# ---- 结果 ----------------------------------------------------------------------------------------------------
print(checks, row.names = FALSE, right = FALSE)
if (all(checks$passed)) cat("\nALL PASS (", nrow(checks), " checks)\n", sep = "") else {
  cat("\nFAILED:", paste(checks$check[!checks$passed], collapse = ", "), "\n"); quit(status = 1)
}
