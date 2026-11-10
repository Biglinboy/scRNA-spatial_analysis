# ============================================================
# 01_spatial_merge_cluster.R —— 空间转录组 RPCA 整合与聚类
# 对应：空间 UMAP、空间域图
# ============================================================

source("functions/utils.R")
load_config("config/config.R")
suppressPackageStartupMessages({ library(Seurat) })

st <- readRDS(st_merged_rds)

# 统一降维与聚类参数，保证与论文一致
st <- FindNeighbors(st, reduction = "pca", dims = st_dims)
st <- FindClusters(st, resolution = st_cluster_res)
st <- RunUMAP(st, reduction = "pca", dims = st_dims)

p_umap <- DimPlot(st, group.by = "seurat_clusters", label = TRUE,
                  raster = FALSE, cols = color_panel) +
  ggplot2::ggtitle("Spatial clusters")
save_fig(p_umap, "01_spatial_umap", subdir = "spatial")

p_dom <- SpatialDimPlot(st, group.by = "seurat_clusters", label = TRUE,
                        pt.size.factor = 1.5)
save_fig(p_dom, "01_spatial_domains", width = 10, height = 8,
         subdir = "spatial")

save_tab(table(st$seurat_clusters, st$orig.ident),
         "01_spatial_cluster_by_sample.tsv")
saveRDS(st, file.path(cache_dir, "01_st_clustered.rds"))