# =============================================================================
# demo_data.R —— 模拟示例数据（固定种子，全部随机生成；不含任何项目数据）
# 每个 demo_*() 返回对应模板的输入；hard = TRUE 给"难"数据（类别多、n 小、未检验的行、并列的交集等），回归测试用它检查护栏。
# =============================================================================

DEMO_LINEAGE <- c("Type A1" = "Lineage A", "Type A2" = "Lineage A", "Type A3" = "Lineage A", "Type B1" = "Lineage B", "Type B2" = "Lineage B",
                  "Type C1" = "Lineage C", "Type D1" = "Lineage D", "Type E1" = "Lineage E")
DEMO_BASE <- c("Lineage A" = "#2166AC", "Lineage B" = "#1B9E77", "Lineage C" = "#E69F00", "Lineage D" = "#CC79A7", "Lineage E" = "#A6761D")
DEMO_GROUP_COLS <- c(Young = "#56B4E9", Middle = "#999999", Old = "#D55E00", `Very old` = "#7F2704")

demo_umap <- function(n = 20000, seed = 42L) {
  set.seed(seed); types <- names(DEMO_LINEAGE)
  ct <- sample(types, n, replace = TRUE, prob = seq(3, 1, length.out = length(types)))
  cx <- stats::setNames(c(-6, -4, -5, 3, 4, 0, 6, -1), types); cy <- stats::setNames(c(3, 5, 1, 4, 2, -4, -2, -6), types)
  d <- data.frame(umap_1 = cx[ct] + stats::rnorm(n, sd = 0.8), umap_2 = cy[ct] + stats::rnorm(n, sd = 0.8), celltype = ct,
                  dataset = sample(c("D1", "D2", "D3"), n, replace = TRUE, prob = c(0.6, 0.3, 0.1)))
  expr <- function(on) ifelse(stats::runif(n) < ifelse(on, 0.7, 0.05), stats::rexp(n, 0.7), 0)
  d$GeneA <- expr(ct %in% c("Type A1", "Type A2")); d$GeneB <- expr(ct == "Type C1")
  d$GeneC <- expr(ct %in% c("Type B1", "Type B2")); d$GeneD <- expr(ct == "Type E1")
  d
}

demo_dot <- function(hard = FALSE, seed = 42L) {
  set.seed(seed); types <- names(DEMO_LINEAGE); if (!hard) types <- types[1:5]
  genes <- unlist(lapply(types, function(t) paste0(sub("Type ", "g", t), c("a", "b", if (hard) "c"))))
  own <- stats::setNames(rep(types, each = if (hard) 3 else 2), genes)
  d <- expand.grid(gene = genes, celltype = types, stringsAsFactors = FALSE)
  on <- own[d$gene] == d$celltype
  d$pct <- round(ifelse(on, stats::runif(nrow(d), 60, 98), stats::runif(nrow(d), 0, 25)), 1)
  d$z <- ifelse(on, stats::runif(nrow(d), 1, 2.5), stats::rnorm(nrow(d), -0.3, 0.4))
  list(d = d, genes = genes, types = types, gene_group = stats::setNames(DEMO_LINEAGE[own], genes),
       ncell = stats::setNames(round(10^stats::runif(length(types), 1.5, 4.5)), types))
}

demo_heat <- function(n_genes = 30, seed = 42L) {
  set.seed(seed); types <- names(DEMO_LINEAGE)
  g <- sprintf("Gene%02d", seq_len(n_genes)); home <- sample(types, n_genes, replace = TRUE)
  m <- matrix(stats::rnorm(n_genes * length(types), 0, 0.5), n_genes, dimnames = list(g, types))
  for (i in seq_len(n_genes)) m[i, home[i]] <- m[i, home[i]] + sample(c(-1, 1), 1) * stats::runif(1, 1.5, 3)
  fdr <- matrix(pmin(1, 10^(-abs(m) * 1.5) * 5), n_genes, dimnames = dimnames(m)); fdr[sample(length(fdr), 5)] <- NA
  list(m = m, fdr = fdr, row_group = stats::setNames(DEMO_LINEAGE[home], g), col_anno = DEMO_LINEAGE)
}

demo_box <- function(hard = FALSE, seed = 42L) {
  set.seed(seed)
  if (!hard) {
    d <- expand.grid(rep = 1:9, group = c("Young", "Old"), panel = c("Type A1", "Type B1", "Type C1"), stringsAsFactors = FALSE)
    d$value <- stats::rnorm(nrow(d), 5, 1) + ifelse(d$group == "Old", 1.2, 0)
    s <- data.frame(panel = c("Type A1", "Type B1", "Type C1"), p = c(2.3e-4, 0.0085, 0.21), g1 = "Young", g2 = "Old")
    return(list(d = d, s = s, cols = DEMO_GROUP_COLS[c("Young", "Old")]))
  }
  # 难：有的组只有 1–2 个样本、有的组缺失、小提琴门槛（5）以下的组
  d <- rbind(data.frame(panel = "Type A1", group = rep(c("Young", "Old", "Very old"), c(3, 4, 4))),
             data.frame(panel = "Type B1", group = rep(c("Old", "Very old"), c(1, 1))),
             data.frame(panel = "Type C1", group = rep(c("Young", "Old", "Very old"), c(2, 5, 4))))
  d$value <- stats::rnorm(nrow(d), 8, 1.5)
  s <- data.frame(panel = c("Type A1", "Type C1"), p = c(0.17, 7.5e-4), g1 = "Young", g2 = "Old")
  list(d = d, s = s, cols = DEMO_GROUP_COLS[c("Young", "Old", "Very old")])
}

demo_forest <- function(hard = FALSE, seed = 42L) {
  set.seed(seed); types <- names(DEMO_LINEAGE)
  est <- stats::rnorm(length(types), 1, 1.2); se <- stats::runif(length(types), 0.2, if (hard) 1.5 else 0.6)
  d <- data.frame(label = types, est = est, lo = est - 1.96 * se, hi = est + 1.96 * se, p = pmin(1, 2 * stats::pnorm(-abs(est / se)) * 3),
                  n_text = "9 / 10", lineage = unname(DEMO_LINEAGE), status = "")
  if (hard) {
    na <- c("Type C1", "Type E1")
    d[d$label %in% na, c("est", "lo", "hi", "p")] <- NA; d$status[d$label %in% na] <- c("too few samples", "gene too low to test")
  }
  d
}

demo_volcano <- function(n = 8000, n_sig = 150, none_sig = FALSE, seed = 42L) {
  set.seed(seed)
  lfc <- stats::rnorm(n, 0, 0.6); p <- stats::runif(n)
  if (!none_sig) { i <- sample(n, n_sig); lfc[i] <- lfc[i] + sample(c(-1, 1), n_sig, TRUE) * stats::runif(n_sig, 0.8, 4); p[i] <- 10^-stats::runif(n_sig, 3, 12) }
  d <- data.frame(gene = sprintf("Gene%04d", seq_len(n)), logFC = lfc, p = p)
  d$fdr <- stats::p.adjust(d$p, "BH"); d$gene[which.min(d$p) + 1] <- "Target1"
  d
}

demo_enrich <- function(n = 10, one_direction = FALSE, long_names = FALSE, seed = 42L) {
  set.seed(seed)
  nm <- if (long_names) paste("Regulation of a fairly long biological process name number", seq_len(n), "in a tissue") else paste("Pathway", LETTERS[seq_len(n)])
  nes <- stats::runif(n, 1.4, 2.6) * if (one_direction) 1 else sample(c(-1, 1), n, TRUE)
  data.frame(term = nm, nes = nes, fdr = 10^-stats::runif(n, 1.4, 6), n_le = sample(8:180, n))
}

demo_grid <- function(n_rows = 12, seed = 42L) {
  set.seed(seed)
  d <- expand.grid(row = paste("Pathway", seq_len(n_rows)), col = names(DEMO_LINEAGE)[1:6], stringsAsFactors = FALSE)
  d$nes <- stats::rnorm(nrow(d), 0, 1.3); d$fdr <- pmin(1, 10^(-abs(d$nes) * 1.6) * 3)
  d$row <- factor(d$row, levels = rev(unique(d$row))); d$col <- factor(d$col, levels = unique(d$col))
  d
}

demo_comp <- function(hard = FALSE, seed = 42L) {
  set.seed(seed)
  comps <- if (hard) paste("Cell type", sprintf("%02d", 1:25)) else names(DEMO_BASE)
  samples <- data.frame(sample = sprintf("S%02d", 1:(if (hard) 30 else 12)))
  samples$group <- factor(rep(c("Young", "Middle", "Old"), length.out = nrow(samples)), levels = c("Young", "Middle", "Old"))
  d <- merge(samples, data.frame(comp = comps))
  d$comp <- factor(d$comp, levels = comps)
  w <- seq(2, 0.1, length.out = length(comps))[as.integer(d$comp)]   # 前面的成分多、后面的少（最后几个是稀有成分）
  d$count <- stats::rpois(nrow(d), stats::runif(nrow(d), 200, 1500) * w)
  pal <- if (hard) stats::setNames(grDevices::hcl.colors(25, "Dark 3"), comps) else DEMO_BASE
  list(d = d, pal = pal, gcols = DEMO_GROUP_COLS[c("Young", "Middle", "Old")])
}

demo_scatter <- function(n = 30, seed = 42L) {
  set.seed(seed)
  d <- data.frame(group = factor(sample(c("Young", "Middle", "Old"), n, TRUE), levels = c("Young", "Middle", "Old")),
                  dataset = sample(c("D1", "D2"), n, TRUE))
  d$x <- stats::rnorm(n, 2, 1.2); d$y <- 0.6 * d$x + stats::rnorm(n, 4, 1)
  d
}

demo_sets <- function(seed = 42L) {
  # 构造了 3 个大小同为 4 的互斥交集（检查 UpSet 不把并列的交集切开）
  set.seed(seed); u <- sprintf("g%03d", 1:400); k <- 0
  take <- function(n) { k <<- k + n; u[(k - n + 1):k] }
  only <- lapply(c(A = 60, B = 45, C = 30, D = 20), take)
  ab <- take(12); abc <- take(8); cd <- take(4); ad <- take(4); bd <- take(4); bc <- take(2)
  list(A = c(only$A, ab, abc, ad), B = c(only$B, ab, abc, bd, bc), C = c(only$C, abc, cd, bc), D = c(only$D, cd, ad, bd))
}

demo_sankey <- function(seed = 42L) {
  set.seed(seed); left <- names(DEMO_LINEAGE)
  right <- c("Author 1", "Author 2", "Author 3", "Author 4", "Author 5", "Author 6", "Author 7a", "Author 7b", "Rare x", "Rare y")
  d <- expand.grid(left = left, right = right, stringsAsFactors = FALSE)
  d$cells <- ifelse(match(d$left, left) == pmin(match(d$right, right), 8), round(stats::runif(nrow(d), 500, 3000)), stats::rpois(nrow(d), 3))
  d
}

demo_gsea <- function(seed = 42L) {
  set.seed(seed); g <- sprintf("Gene%04d", 1:5000)
  r <- stats::setNames(sort(stats::rnorm(5000, 0, 2), decreasing = TRUE), g)
  list(ranks = r, set = c(sample(g[1:600], 40), sample(g, 20)))
}

demo_pca <- function(seed = 42L) {
  set.seed(seed)
  info <- data.frame(sample = sprintf("S%02d", 1:24), Group = rep(c("Young", "Old"), 12), Dataset = rep(c("D1", "D2", "D3"), each = 8))
  m <- matrix(stats::rnorm(3000 * 24, 5, 1), 3000)
  m[1:300, info$Group == "Old"] <- m[1:300, info$Group == "Old"] + 1.5
  m[301:900, info$Dataset == "D3"] <- m[301:900, info$Dataset == "D3"] + 2
  list(logcpm = m, info = info)
}
