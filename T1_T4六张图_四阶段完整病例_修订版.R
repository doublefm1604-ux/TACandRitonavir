# ============================================================
# T1-T4 三个指标的六张图：四阶段完整病例
# 使用前提：先运行“四阶段提取与作图_完整R代码.R”，确保环境中已有
# stage_data 和 output_dir。本脚本不重新提取原始浓度记录，只负责筛选和作图。
# ============================================================

library(tidyverse)

# 1. 基本设置 -------------------------------------------------
stage_levels <- c("T1", "T2", "T3", "T4")
strategy_levels <- c(
  "Early_Discontinuation",
  "Bridging",
  "Concomitant_Use"
)

strategy_colors <- c(
  "Early_Discontinuation" = "#00BA38",
  "Bridging" = "#619CFF",
  "Concomitant_Use" = "#F8766D"
)

metric_labels <- c(
  "OBV" = "他克莫司实测浓度（ng/mL）",
  "CL" = "表观清除率（CL/F）",
  "OBV_IPRED" = "OBV/IPRED 比值"
)

if (!exists("stage_data")) {
  stop("当前环境中不存在 stage_data。请先运行四阶段数据提取代码。")
}
if (!exists("output_dir")) output_dir <- getwd()

plot_output_dir <- file.path(output_dir, "T1_T4_四阶段完整病例作图")
dir.create(plot_output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. 标准化阶段、策略，并检查 ID-阶段唯一性 --------------------
stage_data <- stage_data %>%
  mutate(
    Stage = factor(as.character(Stage), levels = stage_levels),
    Stage_num = as.integer(Stage),
    Strategy = factor(as.character(Strategy), levels = strategy_levels)
  )

stage_duplicates <- stage_data %>%
  count(ID, Stage) %>%
  filter(n > 1)
if (nrow(stage_duplicates) > 0) {
  write.csv(
    stage_duplicates,
    file.path(plot_output_dir, "QC_重复ID阶段组合.csv"),
    row.names = FALSE, fileEncoding = "UTF-8"
  )
  stop("stage_data 中存在同一患者同一阶段的重复记录，已导出 QC 文件。")
}

# 3. 转成长数据：每个实际选中的阶段记录对应三个指标 -------------
plot_data_long_all <- stage_data %>%
  select(ID, Strategy, Stage, Stage_num, OBV, CL, OBV_IPRED) %>%
  pivot_longer(
    cols = c(OBV, CL, OBV_IPRED),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value), is.finite(Value))

# OBV/IPRED 采用 log10 纵轴时必须为正；导出非正值核查表
nonpositive_qc <- plot_data_long_all %>%
  group_by(Metric) %>%
  summarise(
    n_nonpositive = sum(Value <= 0),
    IDs_nonpositive = paste(sort(unique(ID[Value <= 0])), collapse = ", "),
    .groups = "drop"
  )
write.csv(
  nonpositive_qc,
  file.path(plot_output_dir, "绘图数据_非正值核查.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# 4. 分别筛选每个指标 T1-T4 均完整的患者 -------------------------
# OBV、CL、OBV/IPRED 分别筛选；CL 不完整不会影响另外两个指标。
complete_ids_by_metric <- plot_data_long_all %>%
  group_by(Metric, ID) %>%
  summarise(
    n_records = n(),
    n_stages = n_distinct(Stage),
    .groups = "drop"
  ) %>%
  filter(n_records == 4, n_stages == 4) %>%
  select(Metric, ID)

plot_data_long <- plot_data_long_all %>%
  inner_join(complete_ids_by_metric, by = c("Metric", "ID")) %>%
  arrange(Metric, ID, Stage_num)

# 严格验证：每个指标纳入的每名患者都必须恰好有 T1-T4 四条记录
complete_check <- plot_data_long %>%
  group_by(Metric, ID) %>%
  summarise(
    n_records = n(),
    n_stages = n_distinct(Stage),
    .groups = "drop"
  )
if (any(complete_check$n_records != 4 | complete_check$n_stages != 4)) {
  stop("完整病例筛选后仍存在不完整的 ID-Metric 组合。")
}

# 导出各指标完整病例 ID
complete_case_ids <- plot_data_long %>%
  group_by(Metric) %>%
  summarise(
    n_complete_ID = n_distinct(ID),
    IDs = paste(sort(unique(ID)), collapse = ", "),
    .groups = "drop"
  )
write.csv(
  complete_case_ids,
  file.path(plot_output_dir, "各指标_T1_T4完整病例ID.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# 导出各指标、阶段、策略的患者数；每个指标、策略四阶段人数应一致
stage_n <- plot_data_long %>%
  group_by(Metric, Stage, Strategy) %>%
  summarise(n = n_distinct(ID), .groups = "drop") %>%
  complete(
    Metric,
    Stage = factor(stage_levels, levels = stage_levels),
    Strategy = factor(strategy_levels, levels = strategy_levels),
    fill = list(n = 0)
  ) %>%
  arrange(Metric, Stage, Strategy)
write.csv(
  stage_n,
  file.path(plot_output_dir, "绘图样本量_四阶段完整病例.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# 5. 汇总中位数、IQR、P10、P90 --------------------------------
summary_stats <- plot_data_long %>%
  group_by(Metric, Stage, Stage_num, Strategy) %>%
  summarise(
    n = n_distinct(ID),
    P10 = quantile(Value, 0.10, na.rm = TRUE, names = FALSE),
    Q1 = quantile(Value, 0.25, na.rm = TRUE, names = FALSE),
    Median = median(Value, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE, names = FALSE),
    P90 = quantile(Value, 0.90, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )
write.csv(
  summary_stats,
  file.path(plot_output_dir, "三个指标_阶段策略描述统计.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# 6. 绘制六张图：三种策略不同颜色，横坐标不横向错开 -----------
for (m in names(metric_labels)) {

  dat_m <- plot_data_long %>%
    filter(Metric == m) %>%
    arrange(ID, Stage_num)
  sum_m <- summary_stats %>%
    filter(Metric == m) %>%
    arrange(Strategy, Stage_num)

  if (nrow(dat_m) == 0 || nrow(sum_m) == 0) {
    warning("指标 ", m, " 没有四阶段完整病例，跳过绘图。")
    next
  }

  # 个体轨迹图：因为是四阶段完整病例，每位患者均有四个点
  p_individual <- ggplot(
    dat_m,
    aes(x = Stage_num, y = Value, group = ID, color = Strategy)
  ) +
    geom_line(alpha = 0.35, linewidth = 0.60, na.rm = TRUE) +
    geom_point(alpha = 0.85, size = 1.9, na.rm = TRUE) +
    scale_x_continuous(
      breaks = 1:4,
      labels = stage_levels,
      limits = c(0.8, 4.2)
    ) +
    scale_color_manual(values = strategy_colors, drop = FALSE) +
    labs(
      x = "阶段",
      y = metric_labels[[m]],
      color = "用药策略",
      title = paste0(metric_labels[[m]], "：个体变化轨迹"),
      subtitle = paste0(
        "仅纳入该指标 T1-T4 均有数据的患者；n = ",
        n_distinct(dat_m$ID)
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.title = element_text(color = "black")
    )

  # IQR 图：三组共用相同的 T1/T2/T3/T4 横坐标，不做 dodge
  p_iqr <- ggplot(
    sum_m,
    aes(x = Stage_num, group = Strategy)
  ) +
    geom_ribbon(
      aes(ymin = Q1, ymax = Q3, fill = Strategy),
      alpha = 0.16,
      colour = NA
    ) +
    geom_line(
      aes(y = Median, color = Strategy),
      linewidth = 0.9
    ) +
    geom_point(
      aes(y = Median),
      shape = 15,
      size = 2.6,
      color = "black"
    ) +
    geom_line(
      aes(y = P10, color = Strategy),
      linetype = "dashed",
      linewidth = 0.45,
      alpha = 0.8
    ) +
    geom_line(
      aes(y = P90, color = Strategy),
      linetype = "dashed",
      linewidth = 0.45,
      alpha = 0.8
    ) +
    scale_x_continuous(
      breaks = 1:4,
      labels = stage_levels,
      limits = c(0.8, 4.2)
    ) +
    scale_color_manual(values = strategy_colors, drop = FALSE) +
    scale_fill_manual(values = strategy_colors, drop = FALSE) +
    labs(
      x = "阶段",
      y = metric_labels[[m]],
      color = "用药策略",
      fill = "用药策略",
      title = paste0(metric_labels[[m]], "：中位数与四分位距"),
      subtitle = "仅纳入该指标 T1-T4 均有数据的患者；阴影：IQR；黑色方块：中位数；虚线：P10和P90"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.title = element_text(color = "black")
    )

  # OBV/IPRED 范围跨度较大：只改变显示坐标，不改变原始统计值
  if (m == "OBV_IPRED") {
    p_individual <- p_individual + scale_y_log10()
    p_iqr <- p_iqr + scale_y_log10()
  }

  ggsave(
    filename = file.path(plot_output_dir, paste0(m, "_individual_complete.png")),
    plot = p_individual,
    width = 8.5, height = 6, dpi = 300
  )
  ggsave(
    filename = file.path(plot_output_dir, paste0(m, "_IQR_complete.png")),
    plot = p_iqr,
    width = 8.5, height = 6, dpi = 300
  )
}

cat("\n========== 四阶段完整病例绘图完成 ==========\n")
print(complete_case_ids, n = Inf)
cat("\n图片和统计表保存位置：\n", plot_output_dir, "\n")
