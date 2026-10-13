# =============================================================================
# render_examples.R —— 用模拟数据把每个模板的默认样子（和主要开关）画出来，看一眼是否正常
# 运行（任意位置）：Rscript --vanilla <仓库目录>/examples/render_examples.R [输出目录]
# 输出：默认写到 examples/out/（.gitignore 已排除，不进仓库）
# 包库：同 tests/test_templates.R（BIOFIG_LIBS → 本机的 local_libpaths.R → 默认）
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
source(file.path(skill, "load.R"), encoding = "UTF-8"); source(file.path(skill, "examples", "demo_data.R"), encoding = "UTF-8")
out <- commandArgs(TRUE)[1]; if (is.na(out)) out <- file.path(skill, "examples", "out")
sv <- function(p, stem, w, h) save_plot_pub(p, out, stem, w, h, dpi = 300)
pal_u <- lineage_palette(DEMO_LINEAGE, DEMO_BASE)

um <- demo_umap(); um$is_d3 <- um$dataset == "D3"
sv(plot_umap_pub(um, "celltype", pal_u), "01_umap", 89, 85)
sv(umap_highlight_pub(um, "is_d3", title = "Dataset D3"), "02_umap_highlight", 89, 85)
sv(umap_genes_pub(um, c("GeneA", "GeneB", "GeneC", "GeneD")), "03_umap_genes", 120, 110)
z <- demo_dot(TRUE)
sv(dotplot_pub(z$d, gene = "gene", group = "celltype", pct = "pct", value = "z", gene_order = z$genes, group_order = z$types,
               gene_group = z$gene_group, ncell = z$ncell), "04_dotplot", 150, 80)
hz <- demo_heat()
sv(heatmap_pub(hz$m, col_anno = hz$col_anno, anno_cols = DEMO_BASE, z = FALSE, clip = 3, name = "log2FC", row_group = hz$row_group,
               group_cols = DEMO_BASE, group_title = "Largest change in", cluster_columns = FALSE), "05_heatmap", 110, 110)
bn <- demo_box(); bh <- demo_box(TRUE)
sv(box_pub(bn$d, bn$s, bn$cols, "Value per sample"), "06_box", 120, 58)
sv(violin_pub(bn$d, bn$s, bn$cols, "Value per sample"), "07_violin", 120, 58)
sv(box_pub(bh$d, bh$s, bh$cols, "Value per sample"), "08_box_hard", 120, 58)
sv(forest_pub(demo_forest(TRUE), DEMO_BASE, xlab = lab_gene_lfc("Target1")), "09_forest", 120, 60)
sv(volcano_pub(demo_volcano(), target = "Target1"), "10_volcano", 89, 80)
sv(enrich_pub(demo_enrich(), dir_labels = c("Higher in old", "Lower in old")), "11_enrich", 120, 66)
sv(enrich_grid_pub(demo_grid()), "12_enrich_grid", 110, 80)
cn <- demo_comp()
sv(comp_pub(cn$d, cn$pal, cn$gcols), "13_composition", 120, 70)
sv(scatter_pub(demo_scatter(), "Gene X (log2 CPM + 1)", "Gene Y (log2 CPM + 1)", DEMO_GROUP_COLS), "14_scatter", 89, 80)
sv(upset_pub(demo_sets(), min_size = 4, highlight = c("A", "B")), "15_upset", 120, 80)
sv(sankey_pub(demo_sankey(), "left", "right", "cells", pal_u, "Ours", "Author"), "16_sankey", 120, 120)
gz <- demo_gsea(); sv(gsea_curve_pub(gz$ranks, gz$set, "Demo gene set", 1.9, 3e-4), "17_gsea", 89, 80)
pz <- demo_pca(); sv(pca_pub(pz$logcpm, pz$info, "Group", "Dataset", c(Young = "#56B4E9", Old = "#D55E00")), "18_pca", 89, 80)
message("示例图写到 ", out)
