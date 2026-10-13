# 模板目录：数据契约、默认值、开关、来源

每个函数输入一张整理好的表（tidy data frame），输出 ggplot / patchwork / ComplexHeatmap 对象，再用 `save_plot_pub()` 存。
缺列时函数会停下并列出缺哪些列。来源一栏的代码都是**读后重写**的，没有运行或复制原代码；GPL 或无许可证的只借样式。

## forest_pub —— 森林图 / 效应量图
- **输入**：label、est、lo、hi、p、n_text、lineage、status（未检验的行 est 等为 NA，status 写原因）
- **默认**：style = "bar"（谱系色半透明横条 alpha 0.55 + 同色加深 35% 的细边 + 黑色 CI），table = TRUE（n、FDR，显著的加粗），order = "effect"
- **开关**：style = "point"（点 + 区间）；stripes（灰白条纹）；table = FALSE；order = "given"；col_titles
- **来源**：czbiohub/tabula-muris-senis（Nature 2020）；miargentieri/proteomic-age-ukb（Nature Medicine 2024）；broadinstitute/gtex-v8（Science 2020，BSD-3）
- **盲选**：第一批（2026-10-12）用户选 TMS 式横条 + 带表格森林图的表格；要求"柱状图精致一些，加个边"

## volcano_pub —— 火山图
- **输入**：gene、logFC、p、fdr（按样本或混合模型的差异结果）
- **默认**：style = "highlight"；显著 = fdr < 0.05 且 |logFC| ≥ log2(1.5)；y = −log10 P；每个方向标前 10 个；connect = "auto"
- **开关**：style = "outline"（silverbushlab）、"light"（cboix）；connect = "all" / "none"；cap；target；n_lab
- **来源**：cboix/admultiregion_analysis（Nature 2024，无许可证）；silverbushlab（Cell 2026，MIT；改自 SCpubr，GPL-3）；突出式综合了 aifimmunology、sukeniklab、zhumy09
- **盲选**：第一批用户选突出式（和我的推荐一致），引线学 cboix

## enrich_pub / enrich_grid_pub —— 富集图
- **enrich_pub 输入**：term、nes、fdr、n_le；dir_labels 定方向说法
- **默认**：v = "fdr_bar_dot"；点的大小 size_range = c(1.3, 3.2)（数据最小 → 1.3，最大 → 3.2）
- **开关**：v = "nes_bar_dot"（x = NES，末端点颜色 = FDR）、"lollipop"、"dot"
- **ORA（超几何，fgsea::fora 等）没有 NES**：传 nes = 带方向的 log2 富集倍数（下调基因集的结果取负）、n_le = 重叠基因数，
  并设 `effect_label = "log2 FE"`、`effect_title = "log2 fold enrichment"`、`size_title = "Overlap genes"`；上调、下调分开做 ORA 后合并成一张表
- **enrich_grid_pub 输入**：row、col、nes、fdr；默认 style = "dot"（ring_stroke 0.35）；开关 "tile"、"triangle"
- **来源**：linnarsson-lab go_plot、alinaisakova/totalx、heiniglab/dcm_heart_cell_atlas、allonkleinlab
- **盲选**：单个比较——用户对候选都不太满意，要求把点图和横条合并，选 ② FDR 横条；网格——选圆点 + 黑圈，要求黑圈更细

## comp_pub —— 细胞组成
- **输入**：sample、group、comp、count（可选 facet_extra 列，如 dataset）
- **默认**：mode = "bar_n"（每个样本一根条 + 条顶细胞数 + 组色带），soften = 0.25
- **开关**：mode = "facet"、"highlight"（target；占比 < 30% 时自动放大 y 轴）
- **来源**：heiniglab、ding-lab 泛癌、ventolab
- **盲选**：第二批用户选 bar_n，嫌颜色过于鲜艳

## box_pub / violin_pub —— 每个样本一个点的组间比较
- **输入 d**：panel、group、value；**s**：panel、p、g1、g2（按样本模型的结果）；cols = 组 → 颜色（名字定组的顺序）
- **默认**：show_p、show_n；min_n_box = 3、min_n_violin = 5、min_n_p = 3
- **开关**：stripes；violin_pub(style = "framed")（silverbushlab 式）
- **兼容**：box_points_pub(d, x, y, cols) 保留 v0.5 的旧接口
- **来源**：v0.5（Zhuang 组、人脂肪图谱）；gandallab/devbrain_xqtl（Nature 2024）；silverbushlab
- **盲选**：第一批用户要两个版本——箱线用现状，小提琴用 gandal 式

## scatter_pub —— 样本层面相关散点
- **输入**：x、y、group（style = "red" 时可不要），dataset（可选）
- **默认**：style = "black"，统计写在面板上方；method = "pearson"
- **开关**：style = "group"（形状 = 数据集）、"red"（silverbushlab，统计写在点最少的角）；method = "spearman"
- **盲选**：第二批用户选黑线 + 灰带

## heatmap_pub —— 热图
- **输入**：数值矩阵（行 = 基因）；可选 col_anno、row_group、fdr
- **默认**：z = TRUE、clip = 2；行分组用左侧细色带 + 右侧图例（group_title）；行名斜体；格子细灰线
- **开关**：boxes = TRUE（显著黑框，需要 fdr）；group_labels = "block"；z = FALSE + clip（log2FC 矩阵）
- **来源**：sscVis；mjgirgenti PTSD（Nature 2025）；allonkleinlab
- **盲选**：第二批用户选行分组、不要黑框；组名放右侧图例

## dotplot_pub —— 标志基因点图
- **输入**：默认按 Seurat::DotPlot()$data 的列名，可用 gene/group/pct/value 指定；gene_group、ncell 可选
- **默认**：scaled = TRUE（发散色），min_pct = 1；给了 gene_group 就按组分面加括号，给了 ncell 就加右侧细胞数条
- **来源**：v0.5（Chu 2023、KPMP）；HLCA、heiniglab（括号）；aifimmunology（细胞数条）
- **盲选**：第二批用户选括号 + 细胞数条（和我的推荐一致）

## UMAP —— plot_umap_pub、plot_feature_pub、umap_highlight_pub、umap_genes_pub、umap_caption
- v0.5 的样式是卫冕方（大多数图库里的 UMAP 不如它）；新增两个手法：突出子集（theislab、HLCA），斜体基因小图阵列（tgoliver）

## upset_pub、sankey_pub、gsea_curve_pub、pca_pub —— 规范版
- **upset_pub**：输入命名列表；engine = "complexupset"（默认，需要 CRAN 包 ComplexUpset）或 "ggplot"；min_size；highlight
- **sankey_pub**：输入 left、right、weight 三列；需要 ggalluvial；panel_mm 约等于图高减 15 mm
- **gsea_curve_pub**：ranks（命名向量）、set、NES、FDR
- **pca_pub**：logcpm 矩阵、info、colour、shape
