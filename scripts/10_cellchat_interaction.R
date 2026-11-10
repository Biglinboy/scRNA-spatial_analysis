# ============================================================
# 10_cellchat_interaction.R —— 细胞通讯分析
# 对应：new8pair_supply_result 的 6.cellchat_interaction
# ============================================================

source("functions/utils.R")
load_config("config/config.R")
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(NMF); library(ggalluvial)
})

sc <- readRDS(file.path(cache_dir, "07_sc_annotated.rds"))
ensure_dirs(file.path(figures_dir, "cellchat"))

#' 单个分组跑 CellChat
#' @param obj Seurat 对象
#' @param grp 分组名（C 或 P）
run_cellchat <- function(obj, grp) {
  sub <- subset(obj, subset = group == grp)
  if (ncol(sub) < 200) {
    message("细胞数不足，跳过：", grp)
    return(NULL)
  }

  data_input <- GetAssayData(sub, assay = "RNA", slot = "data")
  meta       <- data.frame(labels = sub$annot, row.names = colnames(sub))
  cc <- createCellChat(object = data_input, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.human

  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)

  # 网络图
  p_net <- netVisual_circle(cc@net$count, vertex.weight = rowSums(cc@net$count),
                            weight.scale = TRUE, label.edge = FALSE,
                            title.name = paste0(grp, " interaction count"))
  save_fig(function() p_net, paste0(grp, "_celltype_network"),
           width = 10, height = 10, subdir = "cellchat")

  # 流入流出信号模式
  p_out <- try(netAnalysis_river(cc, pattern = "outgoing"), silent = TRUE)
  if (!inherits(p_out, "try-error")) {
    save_fig(function() p_out, paste0(grp, "_outgoing_river"),
             width = 12, height = 8, subdir = "cellchat")
  }
  p_in <- try(netAnalysis_river(cc, pattern = "incoming"), silent = TRUE)
  if (!inherits(p_in, "try-error")) {
    save_fig(function() p_in, paste0(grp, "_incoming_river"),
             width = 12, height = 8, subdir = "cellchat")
  }

  saveRDS(cc, file.path(cache_dir, paste0("10_cellchat_", grp, ".rds")))
  cc
}

cc_list <- list()
for (grp in group_levels) {
  cc_list[[grp]] <- run_cellchat(sc, grp)
}
message("CellChat 完成。")