# ============================================================
# config.R —— 全项目唯一的路径与参数入口
#
# 项目根目录默认取当前工作目录（在仓库根目录运行 Rscript 即可），
# 也可以用环境变量覆盖：export BJXW_PROJECT_ROOT=/path/to/repo
# 数据目录默认为代码目录的上一级，同样可用 BJXW_DATA_ROOT 覆盖。
# 脚本内不应出现任何绝对路径。
# ============================================================

project_root <- Sys.getenv("BJXW_PROJECT_ROOT", unset = getwd())
data_root    <- Sys.getenv("BJXW_DATA_ROOT",
                           unset = normalizePath(file.path(project_root, ".."),
                                                 mustWork = FALSE))

# ---- 输出目录 -------------------------------------------------
results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")
cache_dir   <- file.path(project_root, "cache")

# ---- 输入数据（按实际布局调整相对位置）------------------------
# 单细胞合并对象
sc_merged_rds   <- file.path(data_root, "Seurat", "new16Harmony",
                             "merge_rds", "merge_data.rds")
# 空间合并对象
st_merged_rds   <- file.path(data_root, "Seurat", "RPCA",
                             "merge_spatial_cluster.rds")
# 上皮/髓系亚群对象（SCENIC 使用）
sc_cna_rds      <- file.path(data_root, "Seurat", "RPCA", "5.RCTD",
                             "CNAs_result.rds")
sc_myeloid_rds  <- file.path(data_root, "Seurat", "RPCA", "5.RCTD",
                             "Myeloid_result.rds")
# 通路基因集（KEGG，两列：Pathway / Gene）
pathway_file    <- file.path(data_root, "script", "getGoTerm",
                             "Human_kegg_gene_sets_symbols.xls")
# SCENIC 结果 loom
scenic_loom     <- file.path(data_root, "Seurat", "RPCA", "14.Scenic",
                             "CNAs", "result", "CNAs_Finished.loom")
# 样本信息表
meta_file       <- file.path(project_root, "metadata", "sample_sheet.csv")

# ---- 全局参数 -------------------------------------------------
seed       <- 42
species    <- "human"
min_cells  <- 3
min_genes  <- 200
max_genes  <- 6000
max_mt_pct <- 15

# ---- 单细胞：整合与聚类 ---------------------------------------
n_pcs        <- 25
n_hvg        <- 2000
dims_use     <- 1:25
batch_key    <- "orig.ident"
cluster_res  <- 0.3
subcluster_res <- 0.3

# ---- 空间：整合、聚类与去卷积 ---------------------------------
st_dims        <- 1:30
st_cluster_res <- 0.6
rctd_max_cores <- 2
rctd_mode      <- "full"          # doublet / full

# ---- 分组：对照与病变（与 sample_sheet.csv 一致）--------------
group_levels <- c("C", "P")
control_ids  <- c("C1", "C2", "C3", "C6", "C7", "C8", "C9", "C11")
case_ids     <- c("P1", "P2", "P3", "P6", "P7", "P8", "P9", "P11")
paired_compare <- c("C1vsP1", "C2vsP2", "C3vsP3", "C6vsP6",
                    "C7vsP7", "C8vsP8", "C9vsP9", "C11vsP11")

# ---- 通路与基因（论文展示用）----------------------------------
pathway_list <- c(
  "hsa04010" = "MAPK signaling pathway",
  "hsa04066" = "HIF-1 signaling pathway",
  "hsa04064" = "NF-kappa B signaling pathway",
  "hsa04714" = "Thermogenesis",
  "hsa04620" = "Toll-like receptor signaling pathway",
  "hsa00010" = "Glycolysis / Gluconeogenesis"
)

gene_list <- c("JUN", "FOS", "CCL3", "CCL4", "CCL3L1", "CCL4L2", "CCR1",
               "CCR5", "CCR7", "TGFB1", "TNF", "IL1B", "TREM2", "PDCD1",
               "CD274", "ARG1", "SPP1", "CD74", "APOE", "NFKB1", "EGFR",
               "SOX2", "PROM1", "CD44", "MKI67", "S100A9", "NOS2", "MRC1",
               "HIF1A", "SLC2A1")

# ---- 绘图 -----------------------------------------------------
fig_dpi    <- 300
fig_width  <- 8
fig_height <- 6
fig_format <- "pdf"               # pdf / png / both

# 与论文一致的离散配色
color_panel <- c("#53A85F", "#58A4C3", "#AB3282", "#8C549C", "#BD956A",
                 "#57C3F3", "#6778AE", "#F3B1A0", "#F1BB72", "#DCC1DD",
                 "#E95C59", "#625D9E", "#F7F398", "#E63863", "#5F3D69",
                 "#C5DEBA", "#CCE0F5", "#B53E2B", "#AA9A59", "#E39A35")

# ---- 依赖 -----------------------------------------------------
required_pkgs <- c(
  "Seurat", "SeuratObject", "dplyr", "ggplot2", "cowplot", "patchwork",
  "spacexr", "STdeconvolve", "reshape2", "argparser", "viridis",
  "RColorBrewer", "paletteer", "Matrix", "pheatmap"
)
