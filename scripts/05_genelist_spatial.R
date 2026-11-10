# ============================================================
# 05_genelist_spatial.R —— 重点基因空间表达与共定位
# 对应：12.genelist_plot 与 13.interaction_gene_plot 两批图
# ============================================================

source("functions/utils.R")
load_config("config/config.R")
suppressPackageStartupMessages({
  library(Seurat); library(viridis); library(ggplot2)
})

st <- readRDS(file.path(cache_dir, "01_st_clustered.rds"))

# ---- 1. 重点基因空间表达 -------------------------------------
genes <- intersect(gene_list, rownames(st))
ensure_dirs(file.path(figures_dir, "genelist"))

for (each in genes) {
  p1 <- SpatialFeaturePlot(st, features = each, pt.size.factor = 1.5,
                           ncol = 3, crop = TRUE) &
    scale_fill_viridis(option = "D", direction = 1)
  ggsave(file.path(figures_dir, "genelist",
                   paste0(each, "_spatial_expression.pdf")),
         p1, width = 25, height = 25)
  ggsave(file.path(figures_dir, "genelist",
                   paste0(each, "_spatial_expression.png")),
         p1, width = 25, height = 25, dpi = fig_dpi)
}

# ---- 2. 共定位基因（JUN-CCL3-CCL4 轴）------------------------
coloc_genes <- c("JUN", "CCL3", "CCL4", "CCL3L1", "CCL4L2",
                 "CCR5", "TGFB1", "TNF", "IL1B")
coloc_genes <- intersect(coloc_genes, rownames(st))
ensure_dirs(file.path(figures_dir, "colocalization"))

for (each in coloc_genes) {
  p1 <- SpatialFeaturePlot(st, features = each, pt.size.factor = 1.5,
                           ncol = 3, crop = TRUE) &
    scale_fill_viridis(option = "D", direction = 1)
  ggsave(file.path(figures_dir, "colocalization",
                   paste0(each, "_spatial_expression.pdf")),
         p1, width = 25, height = 25)
}

# ---- 3. 共表达比例表（spot 层面）-----------------------------
if (length(coloc_genes) >= 2) {
  expr <- GetAssayData(st, assay = "Spatial", slot = "data")[coloc_genes, ,
                                                             drop = FALSE]
  bin <- as.matrix(expr) > 0
  co <- as.data.frame(as.table(bin %*% t(bin)))
  co$pct <- round(co$Freq / ncol(st) * 100, 2)
  colnames(co)[1:2] <- c("gene1", "gene2")
  save_tab(co[, c("gene1", "gene2", "pct")],
           "05_gene_colocalization_pct.tsv")
}

message("基因空间图完成。")