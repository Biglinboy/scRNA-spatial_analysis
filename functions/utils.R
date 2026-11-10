# ============================================================
# utils.R —— 通用函数
# ============================================================

#' 载入全局配置
#' @param config 配置文件路径，默认 config/config.R
load_config <- function(config = "config/config.R") {
  if (!file.exists(config)) stop("找不到配置文件: ", config)
  sys.source(config, envir = globalenv())
  invisible(TRUE)
}

#' 创建输出目录
ensure_dirs <- function(...) {
  for (d in c(...)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

#' 保存图片，统一分辨率与命名
#' @param plot ggplot 对象或返回 ggplot 的函数
#' @param name 文件名（不含扩展名）
save_fig <- function(plot, name, width = fig_width, height = fig_height,
                     subdir = NULL) {
  dir <- if (is.null(subdir)) figures_dir else file.path(figures_dir, subdir)
  ensure_dirs(dir)
  if (is.function(plot)) plot <- plot()
  fmt <- if (exists("fig_format")) fig_format else "pdf"
  if (fmt %in% c("pdf", "both")) {
    ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), plot,
                    width = width, height = height, device = "pdf")
  }
  if (fmt %in% c("png", "both")) {
    ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot,
                    width = width, height = height, dpi = fig_dpi)
  }
  invisible(NULL)
}

#' 保存表格到 results/
save_tab <- function(x, name, row.names = FALSE) {
  ensure_dirs(results_dir)
  utils::write.table(x, file.path(results_dir, name),
                     quote = FALSE, sep = "\t", row.names = row.names)
  invisible(NULL)
}

#' 检查依赖
check_pkgs <- function(pkgs = required_pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    message("缺少依赖：", paste(miss, collapse = ", "))
    message("安装示例：install.packages(c(",
            paste0('"', miss, '"', collapse = ", "), "))")
  } else {
    message("依赖检查通过。")
  }
  invisible(miss)
}

#' 记录会话信息，便于复现
save_session_info <- function(path = "sessionInfo.txt") {
  writeLines(capture.output(sessionInfo()), path)
  invisible(path)
}

#' 按对照/病变给样本打分组标签
#' @param orig_id 样本 ID 向量
add_group <- function(orig_id, control_ids, case_ids,
                      levels = c("C", "P")) {
  g <- ifelse(orig_id %in% control_ids, "C",
              ifelse(orig_id %in% case_ids, "P", NA))
  factor(g, levels = levels)
}

#' 从 KEGG 基因集表中取某条通路的基因
#' @param pathway_df 两列表：Pathway / Gene
#' @param pathway_id 例如 "hsa04010"
get_pathway_genes <- function(pathway_df, pathway_id) {
  unique(pathway_df$Gene[pathway_df$Pathway == pathway_id])
}

#' 把通路名转成文件名安全的形式
safe_name <- function(x) gsub(" | / ", "_", x)
