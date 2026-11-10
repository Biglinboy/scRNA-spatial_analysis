# ============================================================
# 07_sc_merge_annotate.R —— 单细胞合并、Harmony 聚类与注释
# 对应：new8pair_sample_merge 的 majortype 部分
# 说明：cluster 到细胞类型的人工映射请按自己的结果替换
# ============================================================

source("functions/utils.R")
load_config("config/config.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

sc <- readRDS(sc_merged_rds)

# 清理过大的 scale.data，降低内存占用
if ("scale.data" %in% slotNames(sc@assays$RNA)) {
  sc@assays$RNA@scale.data <- matrix()
}

# ---- 1. 聚类 -------------------------------------------------
sc <- FindNeighbors(sc, reduction = "harmony", dims = dims_use)
sc <- FindClusters(sc, resolution = cluster_res)
sc <- RunUMAP(sc, reduction = "harmony", dims = dims_use)

p_umap <- DimPlot(sc, label = TRUE, raster = FALSE, cols = color_panel) +
  ggplot2::ggtitle("Merged clusters")
save_fig(p_umap, "07_sc_umap_clusters", subdir = "singlecell")

# ---- 2. 标志基因特征图 ---------------------------------------
marker_sets <- list(
  Microglia        = c("TMEM119", "P2RY12"),
  Macrophages      = c("CD68", "CD86"),
  `T cells`        = c("CD3G", "CD3E"),
  `B cells`        = c("CD79A", "BANK1"),
  Oligodendrocytes = c("PLP1", "MOG"),
  Endothelial      = c("PECAM1", "CDH5"),
  TAFs             = c("DCN", "FN1"),
  Astrocytes       = c("GFAP", "AQP4"),
  NSCs             = c("SOX4", "DCX"),
  Proliferating    = c("MKI67", "TOP2A")
)

for (nm in names(marker_sets)) {
  feats <- intersect(marker_sets[[nm]], rownames(sc))
  if (length(feats) == 0) next
  p <- FeaturePlot(sc, features = feats, ncol = length(feats),
                   raster = FALSE, label = TRUE) &
    scale_colour_gradientn(colours = c("lightblue", "orange", "red"))
  save_fig(p, paste0("07_marker_", safe_name(nm)),
           width = 10 * length(feats), height = 10, subdir = "singlecell")
}

# ---- 3. 人工注释 ---------------------------------------------
# 按实际 cluster 编号替换；未列出的标为 Others
cluster2celltype <- c(
  "0"  = "Microglia",   "1"  = "Malignant cells", "2"  = "Macrophages",
  "3"  = "Malignant cells", "4" = "T cells",       "5"  = "Microglia",
  "6"  = "TAFs",        "7"  = "Macrophages",     "8"  = "Oligodendrocytes",
  "9"  = "NSCs",        "10" = "Malignant cells", "11" = "B cells",
  "12" = "Endothelial", "13" = "Astrocytes"
)
sc$annot <- unname(cluster2celltype[as.character(sc$seurat_clusters)])
sc$annot[is.na(sc$annot)] <- "Others"
sc$annot <- factor(sc$annot)

celltype_levels <- c("Microglia", "Macrophages", "T cells", "B cells", "TAFs",
                     "Endothelial", "NSCs", "Malignant cells", "Astrocytes",
                     "Oligodendrocytes")
sc$annot <- factor(sc$annot, levels = intersect(celltype_levels,
                                                levels(sc$annot)))

p_anno <- DimPlot(sc, cols = color_panel, label = TRUE, raster = FALSE)
p_split <- DimPlot(sc, cols = color_panel, split.by = "group",
                   label = TRUE, raster = FALSE)
save_fig(p_anno,  "07_umap_annot", subdir = "singlecell")
save_fig(p_split, "07_umap_annot_group", width = 12, height = 7,
         subdir = "singlecell")

# ---- 4. 标志基因点图 -----------------------------------------
dot_genes <- unique(unlist(marker_sets))
dot_genes <- intersect(dot_genes, rownames(sc))
sc$annot <- factor(sc$annot, levels = rev(levels(sc$annot)))
Idents(sc) <- sc$annot

p_dot <- DotPlot(sc, features = dot_genes) +
  scale_color_gradientn(colors = c("white", "#F8766D", "#B2182B")) +
  RotatedAxis() + theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA,
                                    linewidth = 0.8),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)) +
  labs(x = "", y = "")
save_fig(p_dot, "07_annot_dotplot", width = 10, height = 6,
         subdir = "singlecell")

save_tab(table(sc$seurat_clusters, sc$annot), "07_cluster_to_celltype.tsv",
         row.names = TRUE)
saveRDS(sc, file.path(cache_dir, "07_sc_annotated.rds"))
message("单细胞注释完成。")