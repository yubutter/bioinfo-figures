---
name: bioinfo-figures
description: 用户本人的生信作图模板库（R / ggplot2，投稿尺寸、7 pt、Arial），来自 438 个高分文章代码库的筛选和用户两轮盲选，每类图一个"缝合体"函数：默认样子 = 用户选中的，其他好手法做成开关。凡是要画或改：森林图/效应量图、火山图、富集图（GSEA/ORA 单个比较、多组网格）、细胞组成堆叠条、每个样本一个点的箱线/小提琴、样本层面相关散点、热图（ComplexHeatmap）、标志基因点图、UMAP（细胞类型、单基因、突出子集、多基因小图阵列）、UpSet、桑基图、GSEA 曲线、PCA，或要定配色、存投稿图、做"作图盲选"，都先用这个 skill——即使用户只说"画个图""出图""改一下图""配色""图太丑"。单细胞分析流程本身见 single-cell-best-practice。
---

# 生信作图模板（bioinfo-figures）

目标：项目里出图时直接调模板、微调参数，不再每次从零写 ggplot；图既**好看**（用户的审美档案）又**经得起科研追问**（护栏写在代码里）。对话和文档用中文，**图里的文字一律英文**。

## 一、怎么用

```r
# 项目脚本里，先固定 .libPaths()（由项目的 _common.R 负责），再载入（路径换成本 skill 所在目录）：
source("<本 skill 目录>/load.R", encoding = "UTF-8")
p <- forest_pub(d, pal = lineage_cols, xlab = lab_gene_lfc("Gapdh"))
save_plot_pub(p, out_dir, "fig_forest_gapdh", width_mm = FIG_MM[["onehalf"]], height_mm = 70, data = d)   # 同时写 _data.csv
```

- 包：ggplot2、scales、ggrepel、patchwork 必需；ggrastr、cowplot、colorspace、circlize、ComplexUpset、ggalluvial、ComplexHeatmap 可选。`Rscript install_deps.R` 列出缺哪些、各影响哪个函数；**不自动安装**，用户确认后才加 `--install`。缺包时函数会停下并给出安装命令。
- **作者本机**（他人不需要）：可选包分在几个库里，项目要把它们加进 `.libPaths()`：ComplexHeatmap / circlize / colorspace 在 `C:/Rlibs/task2_extra`；ComplexUpset 在 `C:/Rlibs/fig_extra`；ggalluvial、ggpubr 在用户库 `AppData/Local/R/win-library/4.5`（只放在最后补缺）。本 skill 的测试和示例读仓库根目录的 `local_libpaths.R`（不进仓库）。
- 没有 Arial 的机器（常见于 Linux）：设环境变量 `BIOFIG_FONT`（如 `Liberation Sans`）。
- 每个函数的数据契约、默认值、开关、来源见 `references/catalog.md`；视觉约定和已知的坑见 `references/conventions.md`。
- 先看样子：`examples/render_examples.R` 用模拟数据把每类图画一遍。

## 二、必须遵守的规则

1. **按最终尺寸画**：宽度用 `FIG_MM`（单栏 89、1.5 栏 120、双栏 183、Word 155 mm），字号 `FIG_PT = 7`，不画大图再缩小。
2. **统计单位是样本，不是细胞**。箱线/小提琴、森林图、火山图、相关散点的输入都应是按样本（pseudobulk）或混合模型的结果；模板不现场做检验。`volcano_pub` 在 P < 1e-30 时会警告，提醒检查是否按细胞检验了。
3. **护栏不要关掉**：n < 3 不显示 p、不画箱；n < 5 不画小提琴；未检验的行保留并写原因；UpSet 不把大小相同的交集切开。确实要改，在图注里说明。
4. **斜体基因名用 `lab_italic()` / `lab_gene_lfc()`**（plotmath 里必须写成带引号的 `italic("Gapdh")`，否则字母会被排开）。
5. **颜色从项目颜色表来**（细胞类型用谱系色族 `lineage_palette()`），同一个类别在所有图里同一个颜色；不用纯红纯蓝。
6. **每张图都存作图数据**：`save_plot_pub(..., data = )` 写 `<stem>_data.csv`。
7. **不改模板去迁就一次性需求**：先用开关；开关不够，在项目脚本里叠加图层。要改模板本身，按第五节流程。

## 三、模板一览（默认 = 用户 2026-10-12 盲选的结果）

| 函数 | 图 | 默认样子 | 主要开关 |
|---|---|---|---|
| `forest_pub` | 森林图 / 效应量图 | 谱系色半透明横条 + 同色细边 + 黑色 95% CI + 右侧 n / FDR 表格 | `style = "point"`、`stripes`、`table`、`order` |
| `volcano_pub` | 火山图 | 突出式：浅灰底、柔和红蓝、斜体同色标签、目标基因黑圈；引线只在放不下时画 | `style = "outline" / "light"`、`connect`、`cap`、`n_lab` |
| `enrich_pub` | 富集（单个比较） | −log10 FDR 横条（按方向着色）+ 末端空心点（大小 = 基因数，1.3–3.2 mm）+ NES 文字 | `v = "nes_bar_dot" / "lollipop" / "dot"` |
| `enrich_grid_pub` | 富集（多组网格） | 圆点：大小 = −log10 FDR，填充 = NES，显著格细黑圈 | `style = "tile" / "triangle"` |
| `comp_pub` | 细胞组成 | 每个样本一根 100% 堆叠条 + 条顶细胞数 + 组色带，颜色调淡 25% | `mode = "facet" / "highlight"`、`soften` |
| `box_pub` / `violin_pub` | 每个样本一个点的组间比较 | 箱线 = v0.5 现状样式；小提琴 = 半透明 + 同色描边；都带 FDR 括号和每格 n | `min_n_*`、`stripes`、`violin_pub(style = "framed")` |
| `scatter_pub` | 样本层面相关 | 黑线 + 灰带，点按组填色，R / P / n 写在面板上方 | `style = "group" / "red"`、`method` |
| `heatmap_pub` | 热图 | 行分组（左细色带 + 右侧图例）、格子细灰线、斜体基因名 | `boxes`（显著黑框）、`group_labels = "block"`、`z`、`clip` |
| `dotplot_pub` | 标志基因点图 | v0.5 现状 + 基因分组括号 + 右侧细胞数条 | 不传 `gene_group` / `ncell` 即关 |
| `plot_umap_pub` / `plot_feature_pub` | UMAP | v0.5（卫冕方） | `outline`、`labels` |
| `umap_highlight_pub` / `umap_genes_pub` | UMAP 新手法 | 子集着色其余浅灰；斜体基因小图阵列 | — |
| `upset_pub` | UpSet | ComplexUpset，所有 ≥ `min_size` 的交集，关注的交集着色 | `engine = "ggplot"` |
| `sankey_pub` | 桑基图 | 标签在两侧，太近时推开，只有被推开的画折线 | `panel_mm` |
| `gsea_curve_pub` / `pca_pub` | GSEA 曲线、PCA | 规范版 | — |
| `save_plot_pub` | 存图 | PDF（cairo、Arial）+ 600 dpi PNG（+ `_data.csv`） | `legend_separate` |

兼容旧接口：`box_points_pub`、`heatmap_pub` 和 `dotplot_pub` 的旧参数照常可用（single-cell-best-practice v0.5 的脚本不用改）。

## 四、用户的审美与原则（详见 `references/conventions.md` 第二节）

- **整合各家优点**，模板是"缝合体"：好手法都留着，用开关控制，默认值取用户选的。
- 喜欢有颜色层次但收得住；线要细（描边、黑圈偏细）；面积色不要太鲜艳。
- 信息多比少好（n、FDR、细胞数都写上，不要时关掉）。
- 引线只在必要时出现；文字放不下就进图例，不在色块里硬塞。
- 点的大小有上下限：最大的不过大、最小的不过小。
- 同一类图保持一种画法，不因为数据多就换样式。

## 五、新增或修改图类的流程

1. **筛选**：从代码语料（`D:/实战/code_corpus`，图库 `D:/实战/code_corpus_gallery`）找该类图的做法，先过科研检查（作者本机路径，他人换成自己的代码语料）。
2. **候选**：3 个配方 + 现状或默认；**候选太像时有意加 1–2 个明显不同的**（编码、布局或配色不同）。
3. **同数据盲选**：同一份数据（真实项目数据只在本机、不上传）、同一成图尺寸，常规 + 难两种数据；我的推荐先封存（哈希）再出页面。
4. **合并**：用户选中的做默认，其他好手法做开关；按用户的意见改。
5. **写进本 skill**：函数头写来源和许可证（GPL 或无许可证的只借样式、自己重写，不复制代码）；在 `examples/demo_data.R` 加模拟数据；在 `tests/test_templates.R` 加检查（**检查护栏，不只是"能画出来"**）。
6. **跑测试**：`Rscript --vanilla tests/test_templates.R` 必须 ALL PASS；更新 `CHANGELOG.md` 和 `references/catalog.md`；git 提交。

调研和盲选的全部记录在 `D:/实战/skill升级_202610/作图调研二期_202610/`（决策日志与审美档案.md、bench/；作者本机路径，不在仓库里）。
