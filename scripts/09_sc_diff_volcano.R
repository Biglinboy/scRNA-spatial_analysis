# ============================================================
# 09_sc_diff_volcano.R —— 配对样本差异表达与多火山图
# 对应：单个病人差异分析作图（jjVolcano 多面板）
# 说明：差异结果可由上游 diff 流程产出，本脚本负责汇总与作图
# ============================================================

source("functions/utils.R")
load_config("config/config.R")
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(scRNAtoolVis)
})

sc <- readRDS(file.path(cache_dir, "07_sc_annotated.rds"))
diff_dir <- file.path(results_dir, "diff")

ensure_dirs(file.path(figures_dir, "diff"))

for (each_sample in paired_compare) {
  # ---- 1. 汇总该配对下各细胞类型的差异表 ---------------------
  all_diff <- data.frame()
  for (ct in levels(sc$annot)) {
    f <- file.path(diff_dir, each_sample,
                   paste0(ct, "_", each_sample, "_diffgene.xls"))
    if (!file.exists(f)) next
    d <- read.delim(f, header = TRUE, sep = "\t")
    d$cluster <- ct
    all_diff <- rbind(all_diff, d)
  }
  if (nrow(all_diff) == 0) {
    message("跳过（无差异结果）：", each_sample)
    next
  }
  all_diff$cluster <- factor(all_diff$cluster, levels = levels(sc$annot))

  # ---- 2. 多面板火山图 ---------------------------------------
  p1 <- jjVolcano_pvalue(
    diffData        = all_diff,
    pvalue.cutoff   = 0.05,
    adjustP.cutoff  = 0.05,
    log2FC.cutoff   = 0.25,
    col.type        = "adjustP",
    size            = 2.8,
    celltypeSize    = 5,
    tile.col        = color_panel,
    topGeneN        = 10,
    myMarkers       = c("JUNB", "CCL2", "CCL3", "CCL3L1", "CCL4", "CCL4L2")
  )
  save_fig(p1, paste0(each_sample, "_multi_volcano"),
           width = 15, height = 10, subdir = "diff")
}

message("差异火山图完成。")