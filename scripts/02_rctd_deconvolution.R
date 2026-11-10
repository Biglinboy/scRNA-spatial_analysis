# ============================================================
# 02_rctd_deconvolution.R —— RCTD 反卷积
# 参考：单细胞 celltype；目标：每个空间样本
# 对应：major celltype / subtype 权重表与散点饼图
#
# 用法（与原始脚本一致的命令行形式）：
#   Rscript scripts/02_rctd_deconvolution.R \
#     --rds1 <sc_rds> --rds2 <st_rds> \
#     --cluster1 celltype --cluster2 seurat_clusters \
#     --color Y --outdir results/rctd
# ============================================================

source("functions/utils.R")
load_config("config/config.R")
suppressPackageStartupMessages({
  library(Seurat); library(spacexr); library(STdeconvolve)
  library(ggplot2); library(argparser)
})

argv <- arg_parser("RCTD deconvolution")
argv <- add_argument(argv, "--rds1", default = sc_merged_rds,
                     help = "single-cell Seurat rds")
argv <- add_argument(argv, "--rds2", default = st_merged_rds,
                     help = "spatial Seurat rds")
argv <- add_argument(argv, "--cluster1", default = "celltype",
                     help = "cell type column in sc object")
argv <- add_argument(argv, "--cluster2", default = "seurat_clusters",
                     help = "cluster column in spatial object")
argv <- add_argument(argv, "--color", default = "Y",
                     help = "use project palette, Y or N")
argv <- add_argument(argv, "--outdir", default = file.path(results_dir, "rctd"),
                     help = "output dir")
argv <- parse_args(argv)

ensure_dirs(argv$outdir)

# ---- 1. 读入参考与目标 ---------------------------------------
sc <- readRDS(argv$rds1)
sc <- subset(sc, cells = colnames(sc)[!is.na(sc@meta.data[, argv$cluster1])])
sc@meta.data[, argv$cluster1] <- factor(sc@meta.data[, argv$cluster1])
st <- readRDS(argv$rds2)

# 参考矩阵要求整数 counts
ref_counts <- ceiling(sc@assays$RNA@counts)
ref_meta   <- sc@meta.data[, c(argv$cluster1, "nCount_RNA")]
cell_types <- ref_meta[, argv$cluster1]
names(cell_types) <- rownames(ref_meta)
nUMI_ref <- ceiling(ref_meta$nCount_RNA)
names(nUMI_ref) <- rownames(ref_meta)

reference <- Reference(ref_counts, cell_types, nUMI_ref)

# ---- 2. 逐样本运行 -------------------------------------------
index <- 1
for (each_sample in unique(st$orig.ident)) {
  message("RCTD: ", each_sample)
  d1 <- subset(st, orig.ident == each_sample)

  counts <- d1@assays$Spatial@counts
  coords <- d1@images[[index]]@coordinates[, c("imagerow", "imagecol")]
  index <- index + 1

  nUMI <- colSums(counts)
  puck <- SpatialRNA(coords, counts, nUMI)

  myRCTD <- create.RCTD(puck, reference, max_cores = rctd_max_cores)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = rctd_mode)

  saveRDS(myRCTD, file.path(argv$outdir,
                            paste0(each_sample, "_RCTD_", rctd_mode, ".rds")))

  weights      <- myRCTD@results$weights
  norm_weights <- normalize_weights(weights)

  # 列名以数字开头时补 C 前缀，避免下游取列出错
  if (any(grepl("^[0-9]", colnames(norm_weights)))) {
    colnames(norm_weights) <- paste0("C", colnames(norm_weights))
  }

  write.table(norm_weights,
              file.path(argv$outdir,
                        paste0(each_sample, "_spot_norm_celltype_weights.xls")),
              quote = FALSE, sep = "\t", col.names = NA)

  # ---- 3. 散点饼图 -------------------------------------------
  m <- as.matrix(norm_weights)
  p <- as.data.frame(coords)
  colnames(p) <- c("x", "y")

  plt <- vizAllTopics(theta = m, pos = p,
                      topicOrder = seq(ncol(m)),
                      topicCols = if (identical(argv$color, "Y")) color_panel else NULL,
                      groups = NA, group_cols = NA,
                      r = 60, lwd = 0.2,
                      showLegend = TRUE,
                      plotTitle = "scatterpies") +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2))

  ggsave(file.path(argv$outdir, paste0(each_sample, "_spatial_scatterpies.pdf")),
         width = 20, height = 15, plot = plt, bg = "white")
  ggsave(file.path(argv$outdir, paste0(each_sample, "_spatial_scatterpies.png")),
         width = 20, height = 15, plot = plt, bg = "white", dpi = fig_dpi)
}

message("RCTD 完成，结果目录：", argv$outdir)