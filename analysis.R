# ============================================================
# analysis.R —— 综合分析入口
#
# 用法：
#   Rscript analysis.R              跑全部模块
#   Rscript analysis.R 04 05        只跑指定模块
#   Rscript scripts/04_pathway_score.R   单独跑某个脚本
# ============================================================

source("functions/utils.R")
load_config("config/config.R")

args <- commandArgs(trailingOnly = TRUE)

modules <- c(
  "00" = "scripts/00_setup.R",
  "01" = "scripts/01_spatial_merge_cluster.R",
  "02" = "scripts/02_rctd_deconvolution.R",
  "03" = "scripts/03_rctd_spatial_plot.R",
  "04" = "scripts/04_pathway_score.R",
  "05" = "scripts/05_genelist_spatial.R",
  "06" = "scripts/06_scenic_tf_activity.R",
  "07" = "scripts/07_sc_merge_annotate.R",
  "08" = "scripts/08_sc_subcluster.R",
  "09" = "scripts/09_sc_diff_volcano.R",
  "10" = "scripts/10_cellchat_interaction.R"
)

run_ids <- if (length(args) > 0) args else names(modules)
miss <- setdiff(run_ids, names(modules))
if (length(miss) > 0) stop("未知模块编号: ", paste(miss, collapse = ", "))

ensure_dirs(results_dir, figures_dir, cache_dir)
t0 <- Sys.time()

for (id in run_ids) {
  f <- modules[[id]]
  message("\n===== 模块 ", id, " : ", f, " =====")
  if (!file.exists(f)) {
    message("跳过（文件不存在）：", f)
    next
  }
  sys.source(f, envir = new.env(parent = globalenv()))
}

dt <- difftime(Sys.time(), t0, units = "mins")
message("\n全部完成，用时 ", round(as.numeric(dt), 1), " 分钟。")
message("结果目录：", results_dir)
message("成图目录：", figures_dir)
