# ============================================================================
# GO Enrichment Bubble Plot Generation
# ============================================================================
# Purpose:
# This script generates bubble plots for GO enrichment analysis results across
# multiple experimental conditions. It creates three types of plots:
# 1. Basic GO enrichment across 4 groups (C1, C2, S1, S2)
# 2. Differential GO enrichment across 6 comparison groups
# 3. Module-based GO enrichment with top terms per module
#
# Input:
#   Excel files containing GO enrichment results with columns:
#   - GO: GO term ID
#   - Category: GO category (BP, CC, MF)
#   - Description: GO term description
#   - Count: Number of genes in the term
#   - Percent: Percentage of genes
#   - Log10(P): Log10 transformed p-value
#   - Log10(q): Log10 transformed q-value (FDR adjusted)
#
# Output:
#   1. Bubble plot PDF/PNG showing GO enrichment across groups
#   2. Bubble plot for differential enrichment across comparisons
#   3. Bubble plot for module-based enrichment with top terms
#   4. High-resolution PNG image saved to specified output directory
#
# Note: All column names are standardized using janitor::clean_names()
# ============================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(janitor)

# ============================================================================
# Plot 1: Basic GO Enrichment Across 4 Groups
# ============================================================================

# 1) File paths
files <- c(
  C1 = "E:/SCD/下游分析/关键转录因子/富集分析/详细信息/C1.xlsx",
  C2 = "E:/SCD/下游分析/关键转录因子/富集分析/详细信息/C2.xlsx",
  S1 = "E:/SCD/下游分析/关键转录因子/富集分析/详细信息/S1.xlsx",
  S2 = "E:/SCD/下游分析/关键转录因子/富集分析/详细信息/S2.xlsx"
)

# 2) Read and standardize column names
read_one <- function(path, group_name, sheet = 1) {
  read_excel(path, sheet = sheet) %>%
    janitor::clean_names() %>%
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
  mutate(neglog10q = -log10_q)

# 3) Optional: Take top N per group to avoid overcrowding
TOP_N <- 20
plot_df <- df %>%
  group_by(group) %>%
  slice_max(order_by = neglog10q, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# 4) Order GO terms by maximum significance across all groups
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

# 5) Generate bubble plot
p1 <- ggplot(plot_df, aes(x = group, y = description_wrap)) +
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

print(p1)

# ============================================================================
# Plot 2: Differential GO Enrichment Across 6 Comparison Groups
# ============================================================================

# 1) File paths for differential enrichment
files <- c(
  `C1 vs S1 Up` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/C1 VS S1（上调）.xlsx",
  `C1 vs S1 Down` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/C1 VS S1（下调）.xlsx",
  `C2 vs S2 Up` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/C2 VS S2（上调）.xlsx",
  `C2 vs S2 Down` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/C2 VS S2（下调）.xlsx",
  `S1 vs S2 Up` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/S1 VS S2（上调）.xlsx",
  `S1 vs S2 Down` = "E:/SCD/下游分析/度中心/差异富集分析/详细信息/S1 VS S2（下调）.xlsx"
)

# 2) Read and standardize
df <- bind_rows(lapply(names(files), function(nm) read_one(files[[nm]], nm))) %>%
  mutate(
    group = factor(group, levels = names(files)),
    neglog10q = -log10_q
  )

# 3) Take top N per group
TOP_N <- 10
plot_df <- df %>%
  group_by(group) %>%
  slice_max(order_by = neglog10q, n = TOP_N, with_ties = FALSE) %>%
  ungroup()

# 4) Order GO terms
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

# 5) Generate differential enrichment bubble plot
p2 <- ggplot(plot_df, aes(x = group, y = description_wrap)) +
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

print(p2)

# ============================================================================
# Plot 3: Module-based GO Enrichment (Top 5 per Module)
# ============================================================================

# 1) Define module file paths by group
c1_files <- c(
  "Module37" = "E:/SCD/下游分析/功能模块/模块/C1/模块37.xlsx",
  "Module53" = "E:/SCD/下游分析/功能模块/模块/C1/模块53.xlsx",
  "Module4"  = "E:/SCD/下游分析/功能模块/模块/C1/模块4.xlsx",
  "Module22" = "E:/SCD/下游分析/功能模块/模块/C1/模块22.xlsx",
  "Module33" = "E:/SCD/下游分析/功能模块/模块/C1/模块33.xlsx"
)

c2_files <- c(
  "Module11" = "E:/SCD/下游分析/功能模块/模块/C2/模块11.xlsx",
  "Module17" = "E:/SCD/下游分析/功能模块/模块/C2/模块17.xlsx",
  "Module28" = "E:/SCD/下游分析/功能模块/模块/C2/模块28.xlsx",
  "Module33" = "E:/SCD/下游分析/功能模块/模块/C2/模块33.xlsx",
  "Module10" = "E:/SCD/下游分析/功能模块/模块/C2/模块10.xlsx"
)

s1_files <- c(
  "Module6"  = "E:/SCD/下游分析/功能模块/模块/S1/模块6.xlsx",
  "Module9"  = "E:/SCD/下游分析/功能模块/模块/S1/模块9.xlsx",
  "Module44" = "E:/SCD/下游分析/功能模块/模块/S1/模块44.xlsx",
  "Module71" = "E:/SCD/下游分析/功能模块/模块/S1/模块71.xlsx",
  "Module4"  = "E:/SCD/下游分析/功能模块/模块/S1/模块4.xlsx"
)

s2_files <- c(
  "Module9"  = "E:/SCD/下游分析/功能模块/模块/S2/模块9.xlsx",
  "Module27" = "E:/SCD/下游分析/功能模块/模块/S2/模块27.xlsx",
  "Module28" = "E:/SCD/下游分析/功能模块/模块/S2/模块28.xlsx",
  "Module41" = "E:/SCD/下游分析/功能模块/模块/S2/模块41.xlsx",
  "Module45" = "E:/SCD/下游分析/功能模块/模块/S2/模块45.xlsx"
)

# Combine all files with group labels
all_files <- c(c1_files, c2_files, s1_files, s2_files)
group_labels <- c(
  rep("C1", length(c1_files)),
  rep("C2", length(c2_files)),
  rep("S1", length(s1_files)),
  rep("S2", length(s2_files))
)
names(all_files) <- paste0(group_labels, "_", names(all_files))

# 2) Read module data function
read_module_data <- function(path, module_name) {
  read_excel(path, sheet = 1) %>%
    janitor::clean_names() %>%
    mutate(
      module = module_name,
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

# Read all data
df_list <- lapply(names(all_files), function(nm) {
  read_module_data(all_files[[nm]], nm)
})
df <- bind_rows(df_list)

# Calculate -Log10(q) and filter invalid data
df <- df %>%
  filter(!is.na(log10_q) & is.finite(log10_q)) %>%
  mutate(
    neglog10q = -log10_q,
    group = str_extract(module, "^[^_]+")
  )

# 3) Take top 5 GO terms per module
TOP_N_PER_MODULE <- 3
plot_df <- df %>%
  group_by(module) %>%
  slice_max(order_by = neglog10q, n = TOP_N_PER_MODULE, with_ties = FALSE) %>%
  ungroup()

# 4) Set GO term order by overall significance
term_order <- plot_df %>%
  group_by(description) %>%
  summarise(max_sig = max(neglog10q, na.rm = TRUE), .groups = "drop") %>%
  arrange(max_sig)

plot_df <- plot_df %>%
  mutate(
    description = factor(description, levels = term_order$description),
    description = fct_rev(description),
    description_wrap = str_wrap(as.character(description), width = 40)
  )

# Set module order by group
module_order <- plot_df %>%
  distinct(module) %>%
  arrange(module)
plot_df <- plot_df %>%
  mutate(module = factor(module, levels = module_order$module))

# 5) Generate module-based bubble plot
p3 <- ggplot(plot_df, aes(x = module, y = description_wrap)) +
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

print(p3)

# Save high-resolution image
ggsave("E:/SCD/下游分析/功能模块/模块/GO_bubble_plot_top5.png", 
       plot = p3, width = 14, height = 10, dpi = 300)

# Print data summary
cat("\n===== Data Summary =====\n")
cat("Total modules:", n_distinct(plot_df$module), "\n")
cat("Total GO terms:", nrow(plot_df), "\n")
cat("\nGO terms per module:\n")
print(table(plot_df$module))