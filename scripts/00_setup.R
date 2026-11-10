# ============================================================
# 00_setup.R —— 依赖检查、目录创建、随机种子
# ============================================================

source("functions/utils.R")
load_config("config/config.R")

ensure_dirs(results_dir, figures_dir, cache_dir)
set.seed(seed)
options(future.globals.maxSize = 32 * 1024^3)

message("项目根目录：", project_root)
message("数据目录：  ", data_root)

check_pkgs()

# 关键输入文件存在性检查，缺失只提示不中断
for (f in c(sc_merged_rds, st_merged_rds, pathway_file)) {
  if (!file.exists(f)) message("注意：未找到 ", f)
}

save_session_info(file.path(project_root, "sessionInfo.txt"))
message("环境就绪。")