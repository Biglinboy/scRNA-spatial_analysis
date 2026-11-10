# ============================================================
# 03_rctd_spatial_plot.R —— 把 RCTD 权重贴回空间对象并出图
# 对应：各细胞类型空间比例分布图（major / subtype）
# ============================================================

source("functions/utils.R")
load_config("config/config.R")
suppressPackageStartupMessages({
  library(Seurat); library(viridis); library(paletteer); library(ggplot2)
})

st <- readRDS(st_merged_rds)
rctd_dir <- file.path(results_dir, "rctd")

# ---- 1. 汇总所有样本的权重 -----------------------------------
all_w <- data.frame()
for (each in levels(st$orig.ident)) {
  f <- file.path(rctd_dir, paste0(each, "_spot_norm_celltype_weights.xls"))
  if (!file.exists(f)) next
  d <- read.delim(f, header = TRUE, sep = "\t", row.names = 1)
  all_w <- rbind(all_w, d)
}

# ---- 2. 贴回 meta.data ---------------------------------------
raw_ncol <- ncol(st@meta.data)
st <- AddMetaData(st, metadata = all_w)
new_cols <- colnames(st@meta.data)[(raw_ncol + 1):ncol(st@meta.data)]
# 点号替换成下划线，避免 SpatialFeaturePlot 取列名出问题
colnames(st@meta.data)[(raw_ncol + 1):ncol(st@meta.data)] <-
  gsub("\\.", "_", new_cols)
prop_cols <- colnames(st@meta.data)[(raw_ncol + 1):ncol(st@meta.data)]

# ---- 3. 配色与出图 -------------------------------------------
inferno_pal <- c("black",
                 as.vector(paletteer::paletteer_c("viridis::inferno",
                                                  n = 250, direction = 1)))

ensure_dirs(file.path(figures_dir, "rctd"))
for (each in prop_cols) {
  vals <- st@meta.data[, each]
  max_value <- round(max(vals[!is.na(vals)]), 2)

  p1 <- SpatialFeaturePlot(st, features = each, pt.size.factor = 1.5,
                           ncol = 3, crop = TRUE) &
    scale_fill_gradientn(colours = inferno_pal, na.value = "black",
                         limits = c(0, max_value))
  ggsave(file.path(figures_dir, "rctd", paste0(each, "_spatial.pdf")),
         p1, width = 25, height = 20)
  ggsave(file.path(figures_dir, "rctd", paste0(each, "_spatial.png")),
         p1, width = 25, height = 20, dpi = fig_dpi)
}

# ---- 4. 各空间域的平均比例 -----------------------------------
if ("seurat_clusters" %in% colnames(st@meta.data)) {
  avg <- sapply(split(st@meta.data[, prop_cols],
                      st@meta.data$seurat_clusters), function(x) colMeans(x))
  save_tab(t(avg), "03_rctd_proportion_by_domain.tsv", row.names = TRUE)
}

saveRDS(st, file.path(cache_dir, "03_st_with_rctd.rds"))
message("RCTD 空间图完成：", file.path(figures_dir, "rctd"))