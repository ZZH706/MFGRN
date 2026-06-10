# =========================
# 4组GO富集结果气泡图（x=组别, y=GO条目）
# 气泡大小 = Count
# 颜色 = -Log10(q)（蓝 -> 红）
# 输入列要求：GO, Category, Description, Count, %, Log10(P), Log10(q)
# =========================

library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(janitor)

# ---- 1) 文件路径（Windows 路径建议用 / 或者 \\）
files <- c(
  C1 = "E:/SCD/下游分析/关键转录因子/富集分析/详细信息/C1.xlsx",
  C2 = "E:/SCD/下游分析/关键转录因子/富集分析/详细信息/C2.xlsx",
  S1 = "E:/SCD/下游分析/关键转录因子/富集分析/详细信息/S1.xlsx",
  S2 = "E:/SCD/下游分析/关键转录因子/富集分析/详细信息/S2.xlsx"
)

# ---- 2) 读入 + 统一字段名
read_one <- function(path, group_name, sheet = 1) {
  read_excel(path, sheet = sheet) %>%
    janitor::clean_names() %>%  # -> go, category, description, count, percent, log10_p, log10_q
    mutate(group = group_name) %>%
    rename(
      go = go,
      category = category,
      description = description,
      count = count,
      percent = percent,
      log10_p = log10_p,
      log10_q = log10_q
    ) %>%
    mutate(
      count = as.numeric(count),
      percent = as.numeric(percent),
      log10_p = as.numeric(log10_p),
      log10_q = as.numeric(log10_q)
    )
}

df <- bind_rows(lapply(names(files), function(nm) read_one(files[[nm]], nm))) %>%
  mutate(group = factor(group, levels = c("C1", "C2", "S1", "S2"))) %>%
  mutate(neglog10q = -log10_q)   # ---- 颜色用 -Log10(q)

# ---- 3) 可选：每组取Top N条（避免GO条目太多导致图挤）
TOP_N <- 20
plot_df <- df %>%
  group_by(group) %>%
  slice_max(order_by = neglog10q, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# ---- 4) 统一Y轴条目顺序：按所有组里该条目的最大显著性（neglog10q）排序
term_order <- plot_df %>%
  group_by(description) %>%
  summarise(score = max(neglog10q, na.rm = TRUE), .groups = "drop") %>%
  arrange(score)

plot_df <- plot_df %>%
  mutate(
    description = factor(description, levels = term_order$description),
    description = fct_rev(description),
    description_wrap = str_wrap(as.character(description), width = 55)
  )

# ---- 5) 画气泡图（大小=Count，颜色=-Log10(q)，蓝->红）
p <- ggplot(plot_df, aes(x = group, y = description_wrap)) +
  geom_point(aes(size = count, color = neglog10q), alpha = 0.85) +
  scale_size(range = c(2, 10)) +
  scale_color_gradient(
    low = "blue",
    high = "red",
    na.value = "grey80"
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Count",
    color = "-Log10(q)",
    title = "GO Enrichment Bubble Plot"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5)
  )

print(p)











# =========================
# 6组差异富集结果气泡图（x=组别, y=GO条目）
# 气泡大小 = Count
# 颜色 = -Log10(q)（蓝 -> 红）
# 输入列要求：GO, Category, Description, Count, %, Log10(P), Log10(q)
# =========================

library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(janitor)

# ---- 1) 文件路径（Windows 路径建议用 / 或者 \\）
files <- c(
  `C1 vs S1 上调` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/C1 VS S1（上调）.xlsx",
  `C1 vs S1 下调` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/C1 VS S1（下调）.xlsx",
  `C2 vs S2 上调` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/C2 VS S2（上调）.xlsx",
  `C2 vs S2 下调` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/C2 VS S2（下调）.xlsx",
  `S1 vs S2 上调` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/S1 VS S2（上调）.xlsx",
  `S1 vs S2 下调` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/S1 VS S2（下调）.xlsx"
)

# ---- 2) 读入 + 统一字段名
read_one <- function(path, group_name, sheet = 1) {
  read_excel(path, sheet = sheet) %>%
    janitor::clean_names() %>% # -> go, category, description, count, percent, log10_p, log10_q
    mutate(group = group_name) %>%
    rename(
      go = go,
      category = category,
      description = description,
      count = count,
      percent = percent,
      log10_p = log10_p,
      log10_q = log10_q
    ) %>%
    mutate(
      count = as.numeric(count),
      percent = as.numeric(percent),
      log10_p = as.numeric(log10_p),
      log10_q = as.numeric(log10_q)
    )
}

df <- bind_rows(lapply(names(files), function(nm) read_one(files[[nm]], nm))) %>%
  mutate(
    group = factor(group, levels = names(files)),
    neglog10q = -log10_q
  )

# ---- 3) 可选：每组取Top N条（避免GO条目太多导致图挤）
TOP_N <- 10
plot_df <- df %>%
  group_by(group) %>%
  slice_max(order_by = neglog10q, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# ---- 4) 统一Y轴条目顺序：按所有组里该条目的最大显著性（neglog10q）排序
term_order <- plot_df %>%
  group_by(description) %>%
  summarise(score = max(neglog10q, na.rm = TRUE), .groups = "drop") %>%
  arrange(score)

plot_df <- plot_df %>%
  mutate(
    description = factor(description, levels = term_order$description),
    description = fct_rev(description),
    description_wrap = str_wrap(as.character(description), width = 60)
  )

# ---- 5) 画气泡图（大小=Count，颜色=-Log10(q)，蓝->红）
p <- ggplot(plot_df, aes(x = group, y = description_wrap)) +
  geom_point(aes(size = count, color = neglog10q), alpha = 0.85) +
  scale_size(range = c(2, 10)) +
  scale_color_gradient(low = "blue", high = "red", na.value = "grey80") +
  labs(
    x = NULL,
    y = NULL,
    size = "Count",
    color = "-Log10(q)",
    title = "Differential GO Enrichment Bubble Plot (6 Groups)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5)
  )

print(p)







# =========================
# 多模块GO富集结果气泡图（x=模块, y=GO条目）
# 气泡大小 = Count
# 颜色 = -Log10(q)（蓝 -> 红）
# 每个模块取TOP 5个GO条目
# =========================

library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(janitor)

# ---- 1) 定义文件路径（按组别分类）----
# C1组模块
c1_files <- c(
  "模块37" = "E:/SCD/下游分析/功能模块/模块/C1/模块37.xlsx",
  "模块53" = "E:/SCD/下游分析/功能模块/模块/C1/模块53.xlsx",
  "模块4"  = "E:/SCD/下游分析/功能模块/模块/C1/模块4.xlsx",
  "模块22" = "E:/SCD/下游分析/功能模块/模块/C1/模块22.xlsx",
  "模块33" = "E:/SCD/下游分析/功能模块/模块/C1/模块33.xlsx"
)

# C2组模块
c2_files <- c(
  "模块11" = "E:/SCD/下游分析/功能模块/模块/C2/模块11.xlsx",
  "模块17" = "E:/SCD/下游分析/功能模块/模块/C2/模块17.xlsx",
  "模块28" = "E:/SCD/下游分析/功能模块/模块/C2/模块28.xlsx",
  "模块33" = "E:/SCD/下游分析/功能模块/模块/C2/模块33.xlsx",
  "模块10" = "E:/SCD/下游分析/功能模块/模块/C2/模块10.xlsx"
)

# S1组模块
s1_files <- c(
  "模块6"  = "E:/SCD/下游分析/功能模块/模块/S1/模块6.xlsx",
  "模块9"  = "E:/SCD/下游分析/功能模块/模块/S1/模块9.xlsx",
  "模块44" = "E:/SCD/下游分析/功能模块/模块/S1/模块44.xlsx",
  "模块71" = "E:/SCD/下游分析/功能模块/模块/S1/模块71.xlsx",
  "模块4"  = "E:/SCD/下游分析/功能模块/模块/S1/模块4.xlsx"
)

# S2组模块
s2_files <- c(
  "模块9"  = "E:/SCD/下游分析/功能模块/模块/S2/模块9.xlsx",
  "模块27" = "E:/SCD/下游分析/功能模块/模块/S2/模块27.xlsx",
  "模块28" = "E:/SCD/下游分析/功能模块/模块/S2/模块28.xlsx",
  "模块41" = "E:/SCD/下游分析/功能模块/模块/S2/模块41.xlsx",
  "模块45" = "E:/SCD/下游分析/功能模块/模块/S2/模块45.xlsx"
)

# 合并所有文件（带组别标签）
all_files <- c(c1_files, c2_files, s1_files, s2_files)
group_labels <- c(
  rep("C1", length(c1_files)),
  rep("C2", length(c2_files)),
  rep("S1", length(s1_files)),
  rep("S2", length(s2_files))
)
names(all_files) <- paste0(group_labels, "_", names(all_files))

# ---- 2) 读取数据函数 ----
read_module_data <- function(path, module_name) {
  read_excel(path, sheet = 1) %>%
    janitor::clean_names() %>%  # 列名转小写: go, category, description, count, percent, log10_p, log10_q
    mutate(
      module = module_name,  # 完整模块名（如 "C1_模块37"）
      count = as.numeric(count),
      percent = as.numeric(percent),
      log10_p = as.numeric(log10_p),
      log10_q = as.numeric(log10_q)
    ) %>%
    rename(
      go = go,
      category = category,
      description = description,
      count = count,
      percent = percent,
      log10_p = log10_p,
      log10_q = log10_q
    )
}

# 读取所有数据
df_list <- lapply(names(all_files), function(nm) {
  read_module_data(all_files[[nm]], nm)
})
df <- bind_rows(df_list)

# 计算 -Log10(q)，过滤无效数据
df <- df %>%
  filter(!is.na(log10_q) & is.finite(log10_q)) %>%
  mutate(
    neglog10q = -log10_q,
    # 提取组别信息（用于后续排序）
    group = str_extract(module, "^[^_]+")
  )

# ---- 3) 每个模块取TOP 5条GO条目（按显著性排序）----
TOP_N_PER_MODULE <- 3  # 修改为5
plot_df <- df %>%
  group_by(module) %>%
  slice_max(order_by = neglog10q, n = TOP_N_PER_MODULE, with_ties = FALSE) %>%
  ungroup()

# ---- 4) 设置Y轴GO条目顺序（按整体显著性排序）----
term_order <- plot_df %>%
  group_by(description) %>%
  summarise(max_sig = max(neglog10q, na.rm = TRUE), .groups = "drop") %>%
  arrange(max_sig)

plot_df <- plot_df %>%
  mutate(
    description = factor(description, levels = term_order$description),
    description = fct_rev(description),  # 反转使最显著的在上方
    # 长GO描述换行（宽度40字符）
    description_wrap = str_wrap(as.character(description), width = 40)
  )

# 设置X轴模块顺序（按组别排序）
module_order <- plot_df %>%
  distinct(module) %>%
  arrange(module)  # 按字母+数字自然排序
plot_df <- plot_df %>%
  mutate(module = factor(module, levels = module_order$module))

# ---- 5) 绘制气泡图 ----
p <- ggplot(plot_df, aes(x = module, y = description_wrap)) +
  geom_point(aes(size = count, color = neglog10q), alpha = 0.8) +
  scale_size_continuous(range = c(2, 12), breaks = pretty(plot_df$count)) +
  scale_color_gradient(
    low = "blue", 
    high = "red", 
    name = "-Log10(q)",
    na.value = "grey80"
  ) +
  labs(
    x = "Module", 
    y = "GO Term",
    size = "Gene Count",
    title = "GO Enrichment Bubble Plot Across Modules (Top 5 per Module)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# 显示图形
print(p)

# 可选：保存为高分辨率图片
ggsave("E:/SCD/下游分析/功能模块/模块/GO_bubble_plot_top5.png", 
       plot = p, width = 14, height = 10, dpi = 300)

# 输出数据摘要
cat("\n===== 数据摘要 =====\n")
cat("总模块数:", n_distinct(plot_df$module), "\n")
cat("总GO条目数:", nrow(plot_df), "\n")
cat("\n各模块GO条目数（应为5）:\n")
print(table(plot_df$module))