# ============================================================
# 06_scenic_tf_activity.R —— SCENIC 结果解析、RSS 排序与热图
# 输入：pySCENIC 输出的 loom（Regulons / RegulonsAUC）
# 对应：CNAs 与 Myeloid 的 TF rank 图与活性热图
# ============================================================

source("functions/utils.R")
load_config("config/config.R")
suppressPackageStartupMessages({
  library(Seurat); library(SCENIC); library(SCopeLoomR); library(AUCell)
  library(ggplot2); library(ggrepel); library(cowplot); library(pheatmap)
})

#' 对一个亚群对象做主 TF 图
#' @param rds 单细胞对象路径
#' @param celltype_col 亚群列名
#' @param outdir 输出子目录
run_tf_plot <- function(rds, celltype_col, outdir) {
  obj <- readRDS(rds)
  loom <- open_loom(scenic_loom)
  on.exit(close_loom(loom))

  regulonsAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

  ensure_dirs(file.path(figures_dir, outdir))
  Idents(obj) <- obj@meta.data[[celltype_col]]

  for (grp in levels(obj$group)) {
    sub <- subset(obj, group == grp)
    sub@meta.data[[celltype_col]] <- factor(sub@meta.data[[celltype_col]])
    sub <- sub[, colnames(sub) %in% colnames(regulonsAUC)]
    if (ncol(sub) < 20) next

    sub_auc <- regulonsAUC[, match(colnames(sub), colnames(regulonsAUC))]
    sub_auc <- sub_auc[onlyNonDuplicatedExtended(rownames(sub_auc)), ]

    cell_types <- data.frame(row.names = colnames(sub),
                             celltype = sub@meta.data[[celltype_col]])
    cells_by_group <- split(rownames(cell_types), cell_types$celltype)

    # 各亚群的平均 AUC
    auc_by_group <- sapply(cells_by_group, function(cells) {
      rowMeans(getAUC(sub_auc)[, cells])
    })
    auc_scaled <- t(scale(t(auc_by_group), center = TRUE, scale = TRUE))
    auc_scaled <- na.omit(auc_scaled)

    # RSS 排序
    rss <- calcRSS(AUC = getAUC(sub_auc),
                   cellAnnotation = cell_types[colnames(sub_auc), "celltype"])
    rss <- na.omit(rss)

    rank_plots <- list()
    top_tf <- c()
    for (ct in levels(sub@meta.data[[celltype_col]])) {
      if (!ct %in% colnames(rss)) next
      d <- data.frame(TF = rownames(rss), rss = rss[, ct])
      d <- d[order(d$rss, decreasing = TRUE), ]
      d$rank <- seq_len(nrow(d))

      p <- ggplot(d, aes(x = rank, y = rss)) +
        geom_point(size = 3, shape = 16, color = "#1F77B4", alpha = 0.4) +
        geom_point(data = d[1:min(10, nrow(d)), ],
                   size = 3, color = "#DC050C") +
        geom_text_repel(data = d[1:min(10, nrow(d)), ], aes(label = TF),
                        size = 3, fontface = "italic", max.overlaps = 160,
                        box.padding = 0.2, point.padding = 0.3,
                        segment.size = 0.3) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        labs(x = "Regulons Rank", y = "RSS", title = ct)
      rank_plots[[length(rank_plots) + 1]] <- p
      top_tf <- c(top_tf, d$TF[1:min(10, nrow(d))])
    }

    if (length(rank_plots) >= 2) {
      p1 <- plot_grid(plotlist = rank_plots, ncol = min(3, length(rank_plots)))
      ggsave(file.path(figures_dir, outdir,
                       paste0(grp, "_celltype_tf_rank_plot.pdf")),
             p1, width = 20, height = 15)
    }

    top_tf <- intersect(unique(top_tf), rownames(auc_scaled))
    if (length(top_tf) >= 2) {
      pheatmap(auc_scaled[top_tf, , drop = FALSE],
               filename = file.path(figures_dir, outdir,
                                    paste0(grp, "_celltype_tf_heatmap.pdf")),
               cellwidth = 15, cellheight = 15)
    }
  }
  message("SCENIC 图完成：", file.path(figures_dir, outdir))
}

# 两个亚群分别出图
if (file.exists(sc_cna_rds)) {
  run_tf_plot(sc_cna_rds, "subtype", "scenic_epithelial")
}
if (file.exists(sc_myeloid_rds)) {
  run_tf_plot(sc_myeloid_rds, "subtype2", "scenic_myeloid")
}