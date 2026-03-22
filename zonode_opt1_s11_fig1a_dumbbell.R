#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "purrr", "stringr", "tibble", "tidyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    sprintf(
      "缺少必需 R 包：%s。请先安装后再运行脚本。",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

results_dir <- "/mnt/d/其他-小论文/数据库/zonode/陆宏龙/2/优化方案1/结果"
code_dir <- "/mnt/d/其他-小论文/数据库/zonode/陆宏龙/2/优化方案1/代码"

analysis_units_path <- file.path(results_dir, "zonode_opt1_s02_analysis_units.tsv")
asv_dir <- results_dir

script_output_path <- file.path(code_dir, "zonode_opt1_s11_fig1a_dumbbell.R")
source_data_path <- file.path(results_dir, "zonode_opt1_s11_fig1a_source_data.tsv")
png_path <- file.path(results_dir, "zonode_opt1_s11_fig1a_dumbbell.png")
pdf_path <- file.path(results_dir, "zonode_opt1_s11_fig1a_dumbbell.pdf")

required_analysis_cols <- c(
  "analysis_unit", "sample_id", "file_observation_id", "marker"
)
required_asv_cols <- c("asv_id", "count") 
caption_cn <- "原始 FASTQ 读段与去噪后保留读段的对数尺度审计图。颜色表示质控模式，点形状表示 marker。读段计数以对数尺度展示，仅用于质控、排序与审计，不得解释为真实物理组成比例。"
expected_markers <- c("ITS1", "ITS2")

if (!dir.exists(results_dir)) {
  stop(sprintf("结果目录不存在：%s", results_dir), call. = FALSE)
}
if (!file.exists(analysis_units_path)) {
  stop(sprintf("缺少输入文件：%s", analysis_units_path), call. = FALSE)
}
if (!capabilities("cairo")) {
  stop("当前 R 环境不支持 cairo，无法按要求输出 PDF。", call. = FALSE)
}

analysis_units <- readr::read_tsv(analysis_units_path, show_col_types = FALSE)

# 获取输入 FASTQ 质控审计数据
qc_path <- file.path(results_dir, "zonode_opt1_s03_fastq_qc.tsv")
if (file.exists(qc_path)) {
  qc_data <- readr::read_tsv(qc_path, show_col_types = FALSE) %>%
    dplyr::select(analysis_unit, read_count) %>%
    dplyr::mutate(qc_mode = "primary") 
  
  analysis_units <- analysis_units %>%
    dplyr::left_join(qc_data, by = "analysis_unit")
}

missing_analysis_cols <- setdiff(required_analysis_cols, names(analysis_units))
if (length(missing_analysis_cols) > 0) {
  stop(sprintf("analysis_units.tsv 缺少必需列：%s", paste(missing_analysis_cols, collapse = ", ")), call. = FALSE)
}
if (nrow(analysis_units) == 0) stop("analysis_units.tsv 为空。", call. = FALSE)

asv_files <- list.files(asv_dir, pattern = "^zonode_opt1_s05_.*_asv_table\\.tsv$", full.names = TRUE)
if (length(asv_files) == 0) stop(sprintf("目录中未找到任何 s05_asv_table 文件：%s", asv_dir), call. = FALSE)

read_single_asv <- function(path) {
  tbl <- readr::read_tsv(path, show_col_types = FALSE)
  if (nrow(tbl) == 0) return(NULL)
  
  analysis_unit_from_file <- stringr::str_match(basename(path), "^zonode_opt1_s05_(.*)_asv_table\\.tsv$")[, 2]
  tibble(
    analysis_unit = analysis_unit_from_file,
    retained_asv_reads = sum(tbl$count, na.rm = TRUE),
    retained_asv_n = dplyr::n_distinct(tbl$asv_id)
  )
}

asv_summary <- purrr::map_dfr(asv_files, read_single_asv)

# 显式补齐0-read样本：使用 left_join 并用 0 填充 NA
audit_wide <- analysis_units %>%
  dplyr::transmute(
    analysis_unit, sample_id, file_observation_id, marker, qc_mode,
    raw_fastq_reads = read_count
  ) %>%
  dplyr::left_join(asv_summary, by = "analysis_unit") %>%
  dplyr::mutate(
    retained_asv_reads = ifelse(is.na(retained_asv_reads), 0, retained_asv_reads),
    retained_asv_n = ifelse(is.na(retained_asv_n), 0, retained_asv_n)
  )

stage_levels <- c("raw_fastq_reads", "retained_asv_reads")
audit_long <- audit_wide %>%
  tidyr::pivot_longer(
    cols = c(raw_fastq_reads, retained_asv_reads),
    names_to = "stage",
    values_to = "count"
  ) %>%
  dplyr::mutate(
    stage = factor(stage, levels = stage_levels, labels = c("raw_fastq_reads", "retained_asv_reads")),
    stage_index = dplyr::case_when(stage == "raw_fastq_reads" ~ 1, stage == "retained_asv_reads" ~ 2),
    log10_count_plus1 = log10(count + 1)
  )

segment_data <- audit_long %>%
  dplyr::select(analysis_unit, qc_mode, stage, log10_count_plus1) %>%
  tidyr::pivot_wider(names_from = stage, values_from = log10_count_plus1)

# 脱去所有高级属性，剥离为最基础的 data.frame
segment_data_df <- as.data.frame(segment_data)
class(segment_data_df) <- "data.frame"
audit_long_df <- as.data.frame(audit_long)
class(audit_long_df) <- "data.frame"

# ====== 终极修复：直接从 ggplot2 内部提取原生的图层拼接函数，完全绕开全局 '+' 的 S7 冲突 ======
`%++%` <- utils::getFromNamespace("+.gg", "ggplot2")

# 使用自定义的 %++% 拼接所有图层
plot_obj <- ggplot2::ggplot() %++%
  ggplot2::geom_segment(
    data = segment_data_df,
    mapping = ggplot2::aes(
      x = 1, xend = 2,
      y = raw_fastq_reads, yend = retained_asv_reads,
      group = analysis_unit, colour = qc_mode
    ),
    linewidth = 0.6, alpha = 0.8, lineend = "round"
  ) %++%
  ggplot2::geom_point(
    data = audit_long_df,
    mapping = ggplot2::aes(x = stage_index, y = log10_count_plus1, colour = qc_mode, shape = marker, group = analysis_unit),
    size = 2.8, stroke = 0.5
  ) %++%
  ggplot2::scale_x_continuous(
    name = "阶段", breaks = c(1, 2),
    labels = c("原始 FASTQ 读段", "保留 ASV 读段"), limits = c(0.8, 2.2)
  ) %++%
  ggplot2::scale_y_continuous(name = "log10(count + 1)") %++%
  ggplot2::scale_colour_manual(
    name = "质控模式",
    values = c(primary = "#1B9E77", salvage = "#D95F02")
  ) %++%
  ggplot2::scale_shape_manual(
    name = "Marker",
    values = c(ITS1 = 16, ITS2 = 17)
  ) %++%
  ggplot2::labs(
    title = "Figure 1A｜原始 FASTQ 到保留读段的审计图",
    subtitle = "每条连线代表一个 analysis_unit；纵轴为 log10(count + 1)",
    caption = caption_cn
  ) %++%
  ggplot2::theme_bw(base_size = 11) %++%
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
    legend.position = "right",
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = ggplot2::element_text(hjust = 0)
  )

readr::write_tsv(audit_long_df, source_data_path)

ggplot2::ggsave(filename = png_path, plot = plot_obj, width = 8, height = 6, dpi = 300, bg = "white")
ggplot2::ggsave(filename = pdf_path, plot = plot_obj, width = 8, height = 6, device = grDevices::cairo_pdf, bg = "white")

message("Figure 1A 绘制完成。")
