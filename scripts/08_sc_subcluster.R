# ============================================================
# 08_sc_subcluster.R —— 肿瘤与髓系亚群细分
# 对应：new8pair_sample_merge 中 Tumor / Myeloid 两个子目录
# ============================================================

source("functions/utils.R")
load_config("config/config.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

sc <- readRDS(file.path(cache_dir, "07_sc_annotated.rds"))

#' 对指定细胞类型做二次聚类
#' @param obj Seurat 对象
#' @param celltype 要细分的细胞类型
#' @param res 聚类分辨率
#' @param tag 输出文件前缀
subcluster <- function(obj, celltype, res, tag) {
  sub <- subset(obj, subset = annot == celltype)
  if (ncol(sub) < 100) {
    message("细胞数过少，跳过：", celltype)
    return(NULL)
  }
  sub <- NormalizeData(sub, verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = n_hvg, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = min(n_pcs, 30), verbose = FALSE) %>%
    FindNeighbors(dims = 1:min(n_pcs, 30)) %>%
    FindClusters(resolution = res) %>%
    RunUMAP(dims = 1:min(n_pcs, 30))

  sub[[paste0(tag, "_cluster")]] <- sub$seurat_clusters

  p <- DimPlot(sub, label = TRUE, raster = FALSE, cols = color_panel) +
    ggplot2::ggtitle(paste0(celltype, " subclusters"))
  save_fig(p, paste0("08_", tag, "_umap"), subdir = "singlecell")

  save_tab(table(sub$seurat_clusters, sub$group),
           paste0("08_", tag, "_cluster_by_group.tsv"), row.names = TRUE)
  sub
}

tumor   <- subcluster(sc, "Malignant cells", subcluster_res, "tumor")
myeloid <- subcluster(sc, "Myeloid",         subcluster_res, "myeloid")

if (!is.null(tumor))   saveRDS(tumor,   file.path(cache_dir, "08_tumor.rds"))
if (!is.null(myeloid)) saveRDS(myeloid, file.path(cache_dir, "08_myeloid.rds"))
message("亚群细分完成。")