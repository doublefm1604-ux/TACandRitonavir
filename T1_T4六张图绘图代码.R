# ============================================================
# T1-T4 三个指标的六张图：个体轨迹图 + IQR汇总图
# 前提：已经运行前面的四阶段提取代码，当前环境中存在 stage_data 和 output_dir
# 默认使用每个阶段所有可用记录；可切换为四阶段完整病例集。
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. 绘图数据设置
# ------------------------------------------------------------

stage_levels <- c("T1", "T2", "T3", "T4")
strategy_levels <- c(
  "Early_Discontinuation",
  "Bridging",
  "Concomitant_Use"
)

# "available"：每个阶段使用该阶段所有可用记录（默认，T1/T2/T3/T4人数可以不同）
# "complete" ：每个指标仅纳入该指标T1-T4均有数据的患者
plot_mode <- "available"

# 沿用前面约定的策略颜色
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

metric_file_names <- c(
  "OBV" = "OBV",
  "CL" = "CL",
  "OBV_IPRED" = "OBV_IPRED"
)

if (!exists("stage_data")) {
  stop("当前环境中不存在 stage_data。请先运行四阶段数据提取代码。")
}
if (!exists("output_dir")) {
  output_dir <- getwd()
}
if (!plot_mode %in% c("available", "complete")) {
  stop("plot_mode 只能设置为 'available' 或 'complete'。")
}

# 确保阶段顺序和策略顺序正确
stage_data <- stage_data %>%
  mutate(
    Stage = factor(as.character(Stage), levels = stage_levels),
    Stage_num = as.integer(Stage),
    Strategy = factor(as.character(Strategy), levels = strategy_levels)
  )

# ------------------------------------------------------------
# 2. 转成长数据
# ------------------------------------------------------------

plot_data_long_all <- stage_data %>%
  select(ID, Strategy, Stage, Stage_num, OBV, CL, OBV_IPRED) %>%
  pivot_longer(
    cols = c(OBV, CL, OBV_IPRED),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value), is.finite(Value))

# 检查各指标是否存在非正值；OBV/IPRED使用log10坐标，必须为正值
nonpositive_qc <- plot_data_long_all %>%
  group_by(Metric) %>%
  summarise(
    n_nonpositive = sum(Value <= 0),
    IDs_nonpositive = paste(unique(ID[Value <= 0]), collapse = ", "),
    .groups = "drop"
  )
write.csv(
  nonpositive_qc,
  file.path(output_dir, "绘图数据_非正值核查.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 3. 根据绘图模式筛选数据，并导出样本量
# ------------------------------------------------------------

plot_data_long <- plot_data_long_all

if (plot_mode == "complete") {
  # 每个指标分别筛选：该患者在该指标下必须同时具备T1、T2、T3、T4
  complete_ids_by_metric <- plot_data_long_all %>%
    group_by(Metric, ID) %>%
    summarise(n_stage = n_distinct(Stage), .groups = "drop") %>%
    filter(n_stage == 4) %>%
    select(Metric, ID)

  plot_data_long <- plot_data_long_all %>%
    inner_join(complete_ids_by_metric, by = c("Metric", "ID"))
}

# 输出每个指标、阶段、策略的实际样本量
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
  file.path(output_dir, paste0("绘图样本量_", plot_mode, ".csv")),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# 如使用complete模式，同时输出每个指标的完整病例数和ID
if (plot_mode == "complete") {
  complete_case_ids <- plot_data_long %>%
    group_by(Metric) %>%
    summarise(
      n_complete_ID = n_distinct(ID),
      IDs = paste(sort(unique(ID)), collapse = ", "),
      .groups = "drop"
    )
  write.csv(
    complete_case_ids,
    file.path(output_dir, "各指标_T1_T4完整病例ID.csv"),
    row.names = FALSE, na = "", fileEncoding = "UTF-8"
  )
}

# ------------------------------------------------------------
# 4. 描述性统计：中位数、IQR、P10和P90
# ------------------------------------------------------------

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
  file.path(output_dir, paste0("三个指标_阶段策略描述统计_", plot_mode, ".csv")),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# ------------------------------------------------------------
# 5. 六张图
#    个体图：所有实际数据点均显示；只连接相邻阶段，避免跨越缺失阶段连线。
#    IQR图：彩色阴影=Q1-Q3；黑色方块=中位数；虚线=P10/P90。
# ------------------------------------------------------------

# 三组策略在每个阶段内轻微错开，避免IQR汇总图的三组数据完全重叠
strategy_offsets <- c(
  "Early_Discontinuation" = -0.16,
  "Bridging" = 0,
  "Concomitant_Use" = 0.16
)

for (m in names(metric_labels)) {

  dat_m <- plot_data_long %>%
    filter(Metric == m) %>%
    arrange(ID, Stage_num)

  sum_m <- summary_stats %>%
    filter(Metric == m) %>%
    mutate(
      x_plot = Stage_num + unname(strategy_offsets[as.character(Strategy)])
    ) %>%
    arrange(Strategy, Stage_num)

  # 如果某个指标在当前模式下没有完整数据，跳过并提示
  if (nrow(dat_m) == 0 || nrow(sum_m) == 0) {
    warning("指标 ", m, " 没有可用于绘图的数据，已跳过。")
    next
  }

  # 个体轨迹的线段：只连接T1-T2、T2-T3或T3-T4这样的相邻阶段
  individual_segments <- dat_m %>%
    group_by(ID, Strategy) %>%
    arrange(Stage_num, .by_group = TRUE) %>%
    mutate(
      xend = lead(Stage_num),
      yend = lead(Value)
    ) %>%
    filter(
      !is.na(xend),
      xend - Stage_num == 1,
      !is.na(yend)
    ) %>%
    ungroup()

  # 图A：个体散点与相邻阶段连线
  p_individual <- ggplot() +
    geom_segment(
      data = individual_segments,
      aes(
        x = Stage_num, y = Value,
        xend = xend, yend = yend,
        color = Strategy,
        group = ID
      ),
      alpha = 0.30,
      linewidth = 0.55,
      na.rm = TRUE
    ) +
    geom_point(
      data = dat_m,
      aes(x = Stage_num, y = Value, color = Strategy),
      alpha = 0.82,
      size = 1.9,
      na.rm = TRUE
    ) +
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
      subtitle = if (plot_mode == "available") {
        "各阶段使用所有可用记录；仅连接相邻且均有观测的阶段"
      } else {
        "仅纳入该指标T1–T4均有数据的患者"
      }
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.title = element_text(color = "black")
    )

  # 图B：中位数、IQR及P10/P90
  p_iqr <- ggplot(
    sum_m,
    aes(x = x_plot, group = Strategy)
  ) +
    geom_ribbon(
      aes(ymin = Q1, ymax = Q3, fill = Strategy),
      alpha = 0.17,
      colour = NA
    ) +
    geom_line(
      aes(y = Median, color = Strategy),
      linewidth = 0.85
    ) +
    geom_point(
      aes(y = Median),
      shape = 15,
      size = 2.5,
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
      limits = c(0.65, 4.35)
    ) +
    scale_color_manual(values = strategy_colors, drop = FALSE) +
    scale_fill_manual(values = strategy_colors, drop = FALSE) +
    labs(
      x = "阶段",
      y = metric_labels[[m]],
      color = "用药策略",
      fill = "用药策略",
      title = paste0(metric_labels[[m]], "：中位数与四分位距"),
      subtitle = if (plot_mode == "available") {
        "阴影：IQR；黑色方块：中位数；虚线：P10和P90；各点样本量见CSV"
      } else {
        "仅纳入该指标T1–T4均有数据的患者；阴影：IQR；虚线：P10和P90"
      }
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.title = element_text(color = "black")
    )

  # OBV/IPRED跨度较大，采用log10纵坐标；统计量仍在原始尺度上计算
  if (m == "OBV_IPRED") {
    p_individual <- p_individual + scale_y_log10()
    p_iqr <- p_iqr + scale_y_log10()
  }

  ggsave(
    filename = file.path(
      output_dir,
      paste0(metric_file_names[[m]], "_individual_", plot_mode, ".png")
    ),
    plot = p_individual,
    width = 8.5, height = 6, dpi = 300
  )

  ggsave(
    filename = file.path(
      output_dir,
      paste0(metric_file_names[[m]], "_IQR_", plot_mode, ".png")
    ),
    plot = p_iqr,
    width = 8.5, height = 6, dpi = 300
  )
}

cat("\n绘图完成。绘图模式：", plot_mode,
    "\n文件保存位置：\n", output_dir, "\n")
cat("提示：若想改用四阶段完整病例集，将 plot_mode 改为 'complete' 后重新运行本段代码。\n")
