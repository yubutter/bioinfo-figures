# bioinfo-figures

R 生信作图模板库：按投稿尺寸直接出图（宽 89 / 120 / 183 mm，字号 7 pt，Arial），配色统一，常见的统计错误由代码检测并拦下。
每类图一个函数。其他手法保留为参数开关，需要时打开。

Publication-ready figure templates for bioinformatics in R (ggplot2 / ComplexHeatmap), drawn at final journal size with statistical guardrails built in.

## 部署（约 5 分钟）

需要 R ≥ 4.1。测试环境是 Windows 11 + R 4.5.3，用的包版本写在 `install_deps.R` 开头。

**1. 下载**

```bash
git clone https://github.com/yubutter/bioinfo-figures.git
```

没装 git 也可以：仓库页面点 Code → Download ZIP，解压后把文件夹改名为 `bioinfo-figures`。ZIP 里没有 git 记录，`fig_skill_provenance()` 的提交号一栏会是 NA。

**2. 检查缺哪些包**（只检查，不安装）

```bash
Rscript bioinfo-figures/install_deps.R
```

会列出每个包有没有、影响哪个函数。确认后安装：

```bash
Rscript bioinfo-figures/install_deps.R --install
```

- 只装必需的 4 个包：加 `--required`。
- 换 CRAN 镜像：加 `--repo=https://mirrors.tuna.tsinghua.edu.cn/CRAN/`。
- 装到指定的库：加 `--lib=<目录>`。
- ComplexHeatmap 来自 Bioconductor，脚本会先装 BiocManager。

**3. 跑测试**，最后一行应是 `ALL PASS`：

```bash
Rscript --vanilla bioinfo-figures/tests/test_templates.R
```

**4. 在自己的脚本里使用**

```r
source("path/to/bioinfo-figures/load.R", encoding = "UTF-8")   # 嵌套 source 也能找到自己；或先设环境变量 BIOFIG_DIR
p <- volcano_pub(de)                                            # de 的列：gene、logFC、p、fdr
save_plot_pub(p, "out", "fig_volcano", width_mm = FIG_MM[["single"]], height_mm = 75, data = de)
# 得到 fig_volcano.pdf、fig_volcano.png（600 dpi）和作图数据 fig_volcano_data.csv
```

先看样子：`Rscript --vanilla bioinfo-figures/examples/render_examples.R <输出目录>` 用模拟数据把每类图画一遍。

**可能用到的设置**

- **包装在别的库里**：设环境变量 `BIOFIG_LIBS`（多个库用分号隔开）。也可以在仓库根目录建 `local_libpaths.R`，里面写 `.libPaths(c("库1", "库2"))`。这个文件只在本机用，已被 `.gitignore` 排除。测试、示例和 `install_deps.R` 都会读它。
- **没有 Arial**（常见于 Linux）：设环境变量 `BIOFIG_FONT`，如 `Liberation Sans`；或者装上 Arial。
- **更新**：`git pull`，改了什么看 `CHANGELOG.md`。

**作为 Claude Code 的 skill 使用（可选）**：把仓库克隆到 `~/.claude/skills/bioinfo-figures`。之后让 Claude 画这些图时，它会按 `SKILL.md` 调用模板。`SKILL.md` 是中文的，里面有作者的审美偏好；标着"作者本机"的路径别人不需要。

## 图类

| 函数 | 图 |
|---|---|
| `plot_umap_pub`、`plot_feature_pub`、`umap_highlight_pub`、`umap_genes_pub` | UMAP：细胞类型、单基因、突出子集、多基因小图阵列 |
| `dotplot_pub` | 标志基因点图（可加基因分组括号、细胞数条） |
| `heatmap_pub` | 热图（ComplexHeatmap；行分组、显著格黑框） |
| `box_pub`、`violin_pub` | 每个样本一个点的组间比较 |
| `forest_pub` | 森林图 / 效应量图（带 n、FDR 表格） |
| `volcano_pub` | 火山图 |
| `enrich_pub`、`enrich_grid_pub`、`gsea_curve_pub` | 富集：单个比较、多组网格、GSEA 曲线 |
| `comp_pub` | 细胞组成堆叠条 |
| `scatter_pub`、`pca_pub` | 样本层面相关散点、PCA |
| `upset_pub`、`sankey_pub` | UpSet、桑基图 |

辅助函数：
- `umap_point_style()`：按每平方毫米的细胞数自动定点的大小和透明度；
- `tune_grid()`：把几组参数各画一张拼成一页，看图挑值；
- `fig_height_mm()`：按行数定图高；
- `lineage_palette()`：谱系色族；`palette_check()`：配色和色盲检查；
- `lab_italic()`：斜体基因名；
- `fig_skill_provenance()`、`fig_skill_freeze()`：记录版本和提交号；项目交付时冻结一份模板副本。

每个函数的输入列、默认值和开关见 [`references/catalog.md`](references/catalog.md)；视觉约定和已知的坑见 [`references/conventions.md`](references/conventions.md)。

## 内置的护栏

- 统计单位是样本，不是细胞：模板只画传进来的检验结果，不现场检验。火山图 P < 1e-30 时会警告，提醒检查是否按细胞做了检验。
- 每组 n < 3 时不显示 p、不画箱；n < 5 时不画小提琴。
- 森林图里没做检验的行也保留，并写明原因。
- UpSet 不把大小相同的交集切开。
- 缺列时停下，并用中文说明缺哪一列；缺包时停下并给出安装命令，不自动安装。

## 说明

- 测试检查的是护栏有没有生效，不只是"能画出来"。
- 示例和测试数据全部是模拟的，仓库里不放任何项目数据。
- 参考过的公开代码库写在各函数头里。参考代码只读、不运行；所有代码都是重写的。GPL 或无许可证的代码库只借鉴样式。

## 许可证

MIT，见 [LICENSE](LICENSE)。
