# CHANGELOG

## 0.2.1（2026-10-12）

让别人也能部署：
- `load.R` 不再默认作者本机路径：从 `source()` 的路径找自己（嵌套 source 也行），再看 `BIOFIG_DIR`，再看当前目录；找不到或缺必需包时给出中文说明和下一步命令。
- 新增 `install_deps.R`：默认只检查（列出缺的包、影响的函数和将安装到哪个库），加 `--install` 才安装；默认装到个人库（R_LIBS_USER，刚装好 R 时目录还不存在也行）；支持 `--required`、`--lib=`、`--repo=`。
- 缺可选包时的提示从"在作者本机哪个库"改成安装命令（区分 CRAN / Bioconductor）；`save_plot_pub(legend_separate = TRUE)` 缺 cowplot 时也会提示。
- 测试和示例脚本：包库按 `BIOFIG_LIBS` → 本机的 `local_libpaths.R`（不进仓库）→ 默认 `.libPaths()`；从 ZIP 下载（没有 .git）时溯源测试不再误报。
- 字体可用环境变量 `BIOFIG_FONT` 改（没有 Arial 的机器）。
- README 加部署步骤；文档里的示例基因换成中性的 Gapdh。

## 0.2.0（2026-10-12）

- 动态参数：`umap_point_style(n, panel_mm)` 改为按"每平方毫米细胞数"连续插值（89 mm 时与 v0.5 的四个标定点一致），UMAP 三个函数加 `panel_mm`；新增 `tune_grid()`（几组参数并排画一页，看了再挑）和 `fig_height_mm()`（按行数给图高）。测试 61 项。

## 0.1.1（2026-10-12）

- 新增 `fig_skill_provenance()`（版本号 + git 提交号，写进项目的 software_versions.csv）和 `fig_skill_freeze(dest)`（交付前把作图代码复制进项目冻结）。用户选了方案 A：项目通过 single-cell-best-practice 的 sc_helpers.R 指向本 skill。
- `enrich_pub` 支持 ORA 的标签；测试增加到 57 项。

## 0.1.0（2026-10-12）

首个版本。来源：作图调研二期——438 个高分文章代码库的普查和图库初筛，以及用户两轮盲选（记录在作者本机 `D:/实战/skill升级_202610/作图调研二期_202610/`，不在仓库里）。

- 从 single-cell-best-practice v0.5 的 S10 搬来：尺寸与主题、谱系色族、配色自检、UMAP、点图、热图、箱线（保留旧接口）、存图。
- 新增模板，默认样子都是用户盲选的结果：
  - 森林图、火山图；
  - 富集图（单个比较、多组网格）、细胞组成；
  - 箱线 / 小提琴、相关散点；
  - UpSet（ComplexUpset）、桑基图、GSEA 曲线、PCA；
  - UMAP 的突出子集和斜体基因小图阵列。
- 点图加了基因分组括号和细胞数条；热图加了行分组（右侧图例）和显著黑框开关。
- `save_plot_pub()` 加了 `data =`，同时写作图数据 CSV。
- 科研护栏写进代码：
  - n 不够时不显示 p、不画箱或小提琴；
  - 森林图保留未检验的行；
  - 火山图遇到极小 P 时提醒检查是否按细胞检验；
  - UpSet 不把并列的交集切开；
  - 缺列时用中文报错。
- 模拟示例数据（`examples/demo_data.R`）和回归测试（`tests/test_templates.R`，55 项，ALL PASS）。
