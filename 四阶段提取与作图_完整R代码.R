# ============================================================
# 他克莫司与 NMV/r：四阶段数据提取、核查与六张图
# 阶段：T1 基线；T2 P药期；T3 停P药后/Tac重启前；T4 Tac重启后
# 注意：TIME 仅用于筛选和核查；图中横坐标仅显示 T1–T4。
# ============================================================

# -----------------------------
# 0. 准备
# -----------------------------
library(tidyverse)

# 按用户提供的路径读取数据
data_path <- "C:/Mac/Home/Documents/PPKandPBPK/P/DRUGP/1008_60IDs_bys/1008_R_LMM_60.csv"

if (!file.exists(data_path)) {
  stop("找不到数据文件，请检查 data_path：\n", data_path)
}

output_dir <- file.path(dirname(data_path), "FourStage_T1_T4_results")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# 1. 读取与变量类型处理
# -----------------------------
data <- read.csv(
  data_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", ".", "NA", "NaN")
)

required_cols <- c(
  "ID", "TIME", "AMT", "EVID", "OBV", "IPRED", "CLp", "CLover",
  "First_Pax", "Last_Pax", "CLbase", "Strategy"
)
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop("数据缺少必要列：", paste(missing_cols, collapse = ", "))
}

num_vars <- c(
  "ID", "TIME", "AMT", "EVID", "DV", "OBV", "MDV", "ADDL", "II",
  "DOSE", "DAY", "WT", "AGE", "IPRED", "Ptotal", "Combo_Dose",
  "Last_Pax", "First_Pax", "CLbase", "CLp", "CLover", "wtIPRED"
)
num_vars <- intersect(num_vars, names(data))
data <- data %>%
  mutate(across(all_of(num_vars), ~ suppressWarnings(as.numeric(as.character(.x)))))

# 核对 ID、策略标签
if (any(is.na(data$ID))) stop("ID 列存在缺失，请先核查原始数据。")
strategy_levels <- c("Early_Discontinuation", "Bridging", "Concomitant_Use")
unknown_strategy <- setdiff(unique(na.omit(data$Strategy)), strategy_levels)
if (length(unknown_strategy) > 0) {
  warning("发现未预设的 Strategy 标签：", paste(unknown_strategy, collapse = ", "))
}

# 返回每位患者某变量在数据中的第一个非缺失值
first_non_na <- function(x) {
  idx <- which(!is.na(x))[1]
  if (is.na(idx)) return(x[NA_integer_][1])
  x[[idx]]
}

# -----------------------------
# 2. 患者层面信息
#    CLbase 从每位患者全部记录中提取，不限制在浓度记录或某个 TIME。
# -----------------------------
patient_info <- data %>%
  group_by(ID) %>%
  summarise(
    First_Pax = first_non_na(First_Pax),
    Last_Pax = first_non_na(Last_Pax),
    Strategy = first_non_na(Strategy),
    CLbase = first_non_na(CLbase),
    .groups = "drop"
  )

# 检查 CLbase 是否每位患者只有一个不同的非缺失值
clbase_qc <- data %>%
  group_by(ID) %>%
  summarise(
    n_CLbase_distinct = n_distinct(CLbase[!is.na(CLbase)]),
    CLbase_values = paste(unique(CLbase[!is.na(CLbase)]), collapse = "; "),
    .groups = "drop"
  )
write.csv(
  filter(clbase_qc, n_CLbase_distinct != 1),
  file.path(output_dir, "QC_CLbase非缺失值数量异常.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# -----------------------------
# 3. 计算 Tac 重启时间
#    定义：Last_Pax 之后第一条 EVID=1 且 AMT>0 的 Tac 给药记录。
# -----------------------------
restart_info <- data %>%
  select(ID, TIME, EVID, AMT) %>%
  inner_join(patient_info %>% select(ID, Last_Pax), by = "ID") %>%
  filter(
    !is.na(TIME),
    !is.na(Last_Pax),
    EVID == 1,
    !is.na(AMT),
    AMT > 0,
    TIME > Last_Pax
  ) %>%
  group_by(ID) %>%
  summarise(Restart_TIME = min(TIME), .groups = "drop")

patient_info <- patient_info %>%
  left_join(restart_info, by = "ID")

# -----------------------------
# 4. 建立有效浓度观测表
#    仅保留 EVID=0 或 2、OBV>0 的有效浓度记录。
#    CLbase 不从观测行读取，而从 patient_info 按 ID 补入。
# -----------------------------
obs <- data %>%
  filter(
    EVID %in% c(0, 2),
    !is.na(TIME), is.finite(TIME),
    !is.na(OBV), is.finite(OBV), OBV > 0
  ) %>%
  select(ID, TIME, EVID, DV, OBV, IPRED, CLp, CLover) %>%
  left_join(patient_info, by = "ID") %>%
  mutate(
    OBV_IPRED = if_else(
      !is.na(IPRED) & is.finite(IPRED) & IPRED > 0,
      OBV / IPRED,
      NA_real_
    )
  )

# -----------------------------
# 5. 按既定规则提取 T1–T4
#    T1：First_Pax 前最后一条有效浓度；CL=患者层面 CLbase。
#    T2：First_Pax <= TIME < Last_Pax+12；OBV 最大，若并列取最早 TIME；CL=该行 CLp。
#    T3：TIME >= Last_Pax+12 且早于 Tac 重启；取最后一条有效浓度；CL=该行 CLover。
#        若无 Tac 重启记录，则不设上界，取停 P 药 12 h 后最后一条有效浓度。
#    T4：TIME > Restart_TIME；取 Tac 重启后第一条有效浓度；CL=该行 CLover。
#    与 Restart_TIME 完全相同的浓度记录不纳入 T3/T4，另行导出核查。
# -----------------------------
T1 <- obs %>%
  filter(!is.na(First_Pax), TIME < First_Pax) %>%
  group_by(ID) %>%
  mutate(n_candidates = n()) %>%
  arrange(desc(TIME), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(Stage = "T1", CL = CLbase, CL_source = "CLbase")

T2 <- obs %>%
  filter(
    !is.na(First_Pax), !is.na(Last_Pax),
    TIME >= First_Pax,
    TIME < Last_Pax + 12
  ) %>%
  group_by(ID) %>%
  mutate(n_candidates = n()) %>%
  arrange(desc(OBV), TIME, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(Stage = "T2", CL = CLp, CL_source = "CLp")

T3 <- obs %>%
  filter(
    !is.na(Last_Pax),
    TIME >= Last_Pax + 12,
    is.na(Restart_TIME) | TIME < Restart_TIME
  ) %>%
  group_by(ID) %>%
  mutate(n_candidates = n()) %>%
  arrange(desc(TIME), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(Stage = "T3", CL = CLover, CL_source = "CLover")

T4 <- obs %>%
  filter(
    !is.na(Restart_TIME),
    TIME > Restart_TIME
  ) %>%
  group_by(ID) %>%
  mutate(n_candidates = n()) %>%
  arrange(TIME, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(Stage = "T4", CL = CLover, CL_source = "CLover")

stage_levels <- c("T1", "T2", "T3", "T4")
stage_data <- bind_rows(T1, T2, T3, T4) %>%
  mutate(
    Stage = factor(Stage, levels = stage_levels),
    Stage_num = as.integer(Stage),
    Strategy = factor(Strategy, levels = strategy_levels)
  ) %>%
  arrange(ID, Stage)

# 确认每位患者每个阶段至多一条记录
stage_duplicates <- stage_data %>%
  count(ID, Stage) %>%
  filter(n > 1)
if (nrow(stage_duplicates) > 0) {
  stop("发现患者-阶段重复记录，请检查 stage_duplicates。")
}

# -----------------------------
# 6. 导出阶段选择明细与缺失核查表
# -----------------------------
selected_cols <- c(
  "ID", "Strategy", "Stage", "TIME", "EVID", "OBV", "IPRED",
  "OBV_IPRED", "CL", "CL_source", "CLbase", "CLp", "CLover",
  "First_Pax", "Last_Pax", "Restart_TIME", "n_candidates"
)
selected_records <- stage_data %>% select(any_of(selected_cols))
write.csv(
  selected_records,
  file.path(output_dir, "四阶段选中记录明细.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# 每个阶段每个策略的病例数和指标完整数
stage_qc <- stage_data %>%
  group_by(Stage, Strategy) %>%
  summarise(
    n_ID = n_distinct(ID),
    OBV_complete = sum(!is.na(OBV)),
    CL_complete = sum(!is.na(CL)),
    OBV_IPRED_complete = sum(!is.na(OBV_IPRED)),
    CL_missing_ID = paste(ID[is.na(CL)], collapse = ", "),
    .groups = "drop"
  )
write.csv(
  stage_qc,
  file.path(output_dir, "四阶段病例数与指标完整性.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# 所有患者 × 四个阶段：列出没有选出有效浓度记录的患者-阶段组合
all_patient_stage <- tidyr::expand_grid(
  ID = patient_info$ID,
  Stage = stage_levels
)
selected_patient_stage <- stage_data %>%
  distinct(ID, Stage = as.character(Stage))
missing_stage <- all_patient_stage %>%
  anti_join(selected_patient_stage, by = c("ID", "Stage")) %>%
  left_join(patient_info, by = "ID") %>%
  mutate(
    Missing_reason = case_when(
      Stage == "T1" ~ "首次 P 药前无有效浓度记录",
      Stage == "T2" ~ "P 药期定义窗口内无有效浓度记录",
      Stage == "T3" ~ "停 P 药 12 h 后至 Tac 重启前（或之后无重启时的剩余随访期）无有效浓度记录",
      Stage == "T4" & is.na(Restart_TIME) ~ "未识别到 Last_Pax 后的 Tac 重启给药记录",
      Stage == "T4" & !is.na(Restart_TIME) ~ "Tac 重启后无有效浓度记录",
      TRUE ~ "需人工核查"
    )
  ) %>%
  select(ID, Strategy, Stage, First_Pax, Last_Pax, Restart_TIME, Missing_reason)
write.csv(
  missing_stage,
  file.path(output_dir, "未能提取阶段记录的患者.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# 已选中阶段记录中 CL 缺失的具体患者/记录（T2 检查 CLp；T3/T4 检查 CLover）
missing_CL <- stage_data %>%
  filter(is.na(CL)) %>%
  select(
    ID, Strategy, Stage, TIME, OBV, IPRED, OBV_IPRED,
    CL, CL_source, CLbase, CLp, CLover, Last_Pax, Restart_TIME
  ) %>%
  arrange(Stage, ID)
write.csv(
  missing_CL,
  file.path(output_dir, "选中阶段记录_CL缺失患者.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# OBV/IPRED 极端值核查：保留全部数据，不按该表自动剔除病例。
ratio_extreme_qc <- stage_data %>%
  filter(!is.na(OBV_IPRED), is.finite(OBV_IPRED)) %>%
  arrange(desc(OBV_IPRED)) %>%
  slice_head(n = 20) %>%
  select(ID, Strategy, Stage, TIME, OBV, IPRED, OBV_IPRED)
write.csv(
  ratio_extreme_qc,
  file.path(output_dir, "OBV_IPRED最高的20条记录_核查.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# 恰好在 Tac 重启时刻采集的浓度记录：按既定规则不纳入 T3/T4，单独核查
same_time_restart_obs <- obs %>%
  filter(!is.na(Restart_TIME), TIME == Restart_TIME) %>%
  select(ID, Strategy, TIME, Restart_TIME, OBV, IPRED, CLover)
write.csv(
  same_time_restart_obs,
  file.path(output_dir, "QC_浓度时间恰好等于Tac重启时间.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# -----------------------------
# 7. 指标描述性统计：每个阶段×策略
# -----------------------------
plot_data_long <- stage_data %>%
  select(ID, Strategy, Stage, Stage_num, OBV, CL, OBV_IPRED) %>%
  pivot_longer(
    cols = c(OBV, CL, OBV_IPRED),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value), is.finite(Value))

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
  file.path(output_dir, "三个指标_阶段与策略描述性统计.csv"),
  row.names = FALSE, na = "", fileEncoding = "UTF-8"
)

# -----------------------------
# 8. 六张图：三个指标分别绘制个体轨迹图和 IQR 汇总图
# -----------------------------
strategy_colors <- c(
  "Early_Discontinuation" = "#00BA38",
  "Bridging" = "#619CFF",
  "Concomitant_Use" = "#F8766D"
)

metric_labels <- c(
  "OBV" = "Observed tacrolimus concentration (ng/mL)",
  "CL" = "Apparent clearance (CL/F)",
  "OBV_IPRED" = "OBV/IPRED ratio"
)
metric_file_names <- c(
  "OBV" = "OBV",
  "CL" = "CL",
  "OBV_IPRED" = "OBV_IPRED"
)

for (m in names(metric_labels)) {

  dat_m <- plot_data_long %>% filter(Metric == m)
  sum_m <- summary_stats %>% filter(Metric == m)

  # 图 1：个体散点/连线图。每位患者连接其实际具备的阶段观测；缺失阶段不补值。
  p_individual <- ggplot(
    dat_m,
    aes(x = Stage_num, y = Value, group = ID, color = Strategy)
  ) +
    geom_line(alpha = 0.30, linewidth = 0.55, na.rm = TRUE) +
    geom_point(alpha = 0.85, size = 1.8, na.rm = TRUE) +
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
      title = paste0(metric_labels[[m]], "：个体变化轨迹")
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(color = "black")
    )

  # 图 2：阶段×策略的中位数与 IQR；虚线表示 P10/P90。
  p_iqr <- ggplot(
    sum_m,
    aes(x = Stage_num, group = Strategy)
  ) +
    geom_ribbon(
      aes(ymin = Q1, ymax = Q3, fill = Strategy),
      alpha = 0.15,
      colour = NA
    ) +
    geom_line(
      aes(y = Median, color = Strategy),
      linewidth = 0.9
    ) +
    geom_point(
      aes(y = Median, color = Strategy),
      shape = 15,
      size = 2.8
    ) +
    geom_line(
      aes(y = P10, color = Strategy),
      linetype = "dashed",
      linewidth = 0.45,
      alpha = 0.75
    ) +
    geom_line(
      aes(y = P90, color = Strategy),
      linetype = "dashed",
      linewidth = 0.45,
      alpha = 0.75
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
      title = paste0(metric_labels[[m]], "：中位数与四分位距")
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(color = "black")
    )

  # OBV/IPRED 在本数据中的跨度较大，采用 log10 坐标显示；不改变原始值和统计量。
  if (m == "OBV_IPRED") {
    p_individual <- p_individual + scale_y_log10()
    p_iqr <- p_iqr + scale_y_log10()
  }

  ggsave(
    filename = file.path(output_dir, paste0(metric_file_names[[m]], "_individual.png")),
    plot = p_individual,
    width = 8.5, height = 6, dpi = 300
  )
  ggsave(
    filename = file.path(output_dir, paste0(metric_file_names[[m]], "_IQR.png")),
    plot = p_iqr,
    width = 8.5, height = 6, dpi = 300
  )
}

# -----------------------------
# 9. 控制台输出关键核对信息
# -----------------------------
cat("\n========== 数据概况 ==========" , "\n")
cat("总行数：", nrow(data), "\n")
cat("患者数：", n_distinct(data$ID), "\n")
cat("有效浓度记录数（EVID 0/2 且 OBV>0）：", nrow(obs), "\n")

cat("\n========== 各阶段病例数 ==========" , "\n")
print(stage_data %>% count(Stage, name = "患者数"))

cat("\n========== 各阶段×策略完整性 ==========" , "\n")
print(stage_qc, n = Inf)

cat("\n========== 已选中记录中 CL 缺失 ==========" , "\n")
if (nrow(missing_CL) == 0) {
  cat("无缺失。\n")
} else {
  print(missing_CL %>% select(ID, Strategy, Stage, TIME, CL_source, CLp, CLover), n = Inf)
}

cat("\n========== 无法提取阶段记录的患者 ==========" , "\n")
print(missing_stage, n = Inf)

cat("\n全部结果已保存至：\n", output_dir, "\n")
