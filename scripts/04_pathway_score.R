# ============================================================
# 04_pathway_score.R —— KEGG 通路模块打分与空间可视化
# 对应：MAPK / HIF-1 / NF-kB / Thermogenesis / TLR / Glycolysis
#       六条通路的空间得分图
# ============================================================

source("functions/utils.R")
load_config("config/config.R")
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })

st <- readRDS(file.path(cache_dir, "01_st_clustered.rds"))
pathway_df <- read.delim(pathway_file, header = TRUE, sep = "\t")

ensure_dirs(file.path(figures_dir, "pathway"))

score_tab <- data.frame(row.names = colnames(st))

for (i in seq_along(pathway_list)) {
  pid  <- names(pathway_list)[i]
  name <- pathway_list[[i]]

  genes <- get_pathway_genes(pathway_df, pid)
  genes <- intersect(genes, rownames(st))
  if (length(genes) < 5) {
    message("通路基因过少，跳过：", name)
    next
  }

  st <- AddModuleScore(st, features = list(genes),
                       name = paste0(pid, "_score"))
  score_col <- paste0(pid, "_score1")
  new_name  <- paste0(safe_name(name), "_score")
  colnames(st@meta.data)[colnames(st@meta.data) == score_col] <- new_name
  score_tab[[new_name]] <- st@meta.data[[new_name]]

  p1 <- SpatialFeaturePlot(st, features = new_name, pt.size.factor = 1.5,
                           ncol = 3, crop = TRUE)
  ggsave(file.path(figures_dir, "pathway",
                   paste0(safe_name(name), "_spatial_score.pdf")),
         p1, width = 25, height = 25)
  ggsave(file.path(figures_dir, "pathway",
                   paste0(safe_name(name), "_spatial_score.png")),
         p1, width = 25, height = 25, dpi = fig_dpi)
}

save_tab(score_tab, "04_pathway_score_per_spot.tsv", row.names = TRUE)

# 与分组的比较（Wilcoxon）
if ("group" %in% colnames(st@meta.data)) {
  res <- lapply(colnames(score_tab), function(col) {
    d <- st@meta.data
    if (length(unique(na.omit(d$group))) < 2) return(NULL)
    data.frame(feature = col,
               p_value = wilcox.test(d[[col]] ~ d$group, exact = FALSE)$p.value)
  })
  res <- do.call(rbind, res)
  if (!is.null(res)) {
    res$p_adj <- p.adjust(res$p_value, method = "BH")
    save_tab(res, "04_pathway_group_wilcox.tsv")
  }
}

saveRDS(st, file.path(cache_dir, "04_st_pathway.rds"))
message("通路打分完成。")