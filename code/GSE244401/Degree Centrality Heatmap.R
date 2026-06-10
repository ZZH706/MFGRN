###############################################################
#  Degree centrality — 仅转录因子 (TF) 四网络对比
#  —— 文件自带表头版 —— 2025‑05
###############################################################

library(igraph)
library(dplyr)

# ---------- ① TF 白名单 ----------
tf_list <- read.table("E:/1数据/Homo_TF_clean.txt",
                      header = FALSE,
                      stringsAsFactors = FALSE)
all_tfs <- unique(trimws(tf_list[[1]]))   # TF gene‑symbol 列表

# ---------- ② 网络文件 ----------
network_files <- c(
  "E:/1数据/≥2/c1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/c2_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s2_TF_TARGET_cleaned.txt"
)
network_names <- c("Control_T1", "Control_T2", "SCA_T1", "SCA_T2")

degree_centrality_list <- list()

# ---------- ③ 循环计算 ----------
for (i in seq_along(network_files)) {
  
  # 1) 读取边表（文件自带列名）
  edge_data <- read.table(
    network_files[i],
    header = TRUE,        # ✅ 文件已有表头
    sep    = "\t",
    stringsAsFactors = FALSE
  )
  # 如果列名并非 TF/Target/Weight，可在这里改名：
  # colnames(edge_data)[1:3] <- c("TF","Target","Weight")
  
  # 2) 构建 igraph
  g <- graph_from_data_frame(edge_data[, 1:2], directed = TRUE)
  
  # 3) 仅取白名单且位于图中的 TF
  tf_in_graph <- intersect(all_tfs, V(g)$name)
  
  degree_values <- degree(g, mode = "all")[tf_in_graph]
  degree_values <- degree_values[!is.na(degree_values)]
  
  degree_df <- data.frame(
    Gene    = names(degree_values),
    Degree  = as.numeric(degree_values),
    Network = network_names[i],
    stringsAsFactors = FALSE
  )
  
  degree_centrality_list[[i]] <- degree_df
  
  # 4) 单网结果按 Degree 降序保存
  write.table(
    degree_df[order(-degree_df$Degree), ],
    file = paste0("degree_centrality_", network_names[i], ".txt"),
    sep  = "\t", quote = FALSE, row.names = FALSE
  )
}

# ---------- ④ 合并四网结果 ----------
all_degree_data <- bind_rows(degree_centrality_list)

# ---------- ⑤ 绘制分布竖线图 ----------
pdf("degree_centrality_TF_only.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))

for (nw in network_names) {
  
  sub <- all_degree_data %>% 
    filter(Network == nw) %>% 
    arrange(desc(Degree))
  
  plot(
    sub$Degree,
    type = "h", col = "black",
    main = paste("TF Degree –", nw),
    xlab = "TFs (sorted)",
    ylab = "Degree",
    lwd  = 1
  )
}
dev.off()

cat("\n✅ 仅 TF 的 Degree 计算完成！\n",
    "- 每网结果：degree_centrality_<Network>.txt\n",
    "- 合并表   ：all_degree_data (在 R 环境)\n",
    "- 图文件   ：degree_centrality_TF_only.pdf\n")

# 对度中心性取 log10
degree_log_matrix <- log10(degree_matrix + 1)

# 画对数热图
pdf("TF_degree_heatmap_log10.pdf", width = 8, height = 10)
pheatmap(degree_log_matrix,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 6)
dev.off()

library(tidyr)

# 把 all_degree_data 转为宽格式矩阵
degree_wide <- pivot_wider(all_degree_data,
                           names_from = Network,
                           values_from = Degree,
                           values_fill = 0)  # 没有的填 0

# 转成 matrix，去掉第一列 Gene
degree_matrix <- as.matrix(degree_wide[,-1])
rownames(degree_matrix) <- degree_wide$Gene

# 取 log10 画热图
degree_log_matrix <- log10(degree_matrix + 1)

library(pheatmap)

pdf("TF_degree_heatmap_log10_no_rownames.pdf", width = 8, height = 10)
pheatmap(degree_log_matrix,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 6,
         show_rownames = FALSE)   # 不显示左侧行名（TF）
dev.off()

library(tidyr)
library(pheatmap)

# 1. 把 all_degree_data 转为宽格式矩阵
degree_wide <- pivot_wider(all_degree_data,
                           names_from = Network,
                           values_from = Degree,
                           values_fill = 0)

degree_matrix <- as.matrix(degree_wide[,-1])
rownames(degree_matrix) <- degree_wide$Gene

# 2. 计算整体热图的 log10 矩阵
degree_log_matrix <- log10(degree_matrix + 1)

# 3. 整体热图用于提取行顺序
res_all <- pheatmap(degree_log_matrix, silent = TRUE)
row_order <- res_all$tree_row$order
ordered_genes <- rownames(degree_log_matrix)[row_order]

# 4. 定义对比组和输出文件名
compare_groups <- list(
  c("Control_T1", "SCA_T1"),
  c("Control_T2", "SCA_T2"),
  c("SCA_T1", "SCA_T2"),
  c("Control_T1", "Control_T2")
)

file_names <- c(
  "TF_degree_heatmap_ControlT1_vs_SCAT1.pdf",
  "TF_degree_heatmap_ControlT2_vs_SCAT2.pdf",
  "TF_degree_heatmap_SCA_T1_vs_T2.pdf",
  "TF_degree_heatmap_Control_T1_vs_T2.pdf"
)

# 5. 画子热图（瘦、去标题、横着列名）
for (i in seq_along(compare_groups)) {
  cols <- compare_groups[[i]]
  sub_matrix <- degree_matrix[, cols, drop = FALSE]
  sub_matrix_log <- log10(sub_matrix + 1)
  
  sub_matrix_log_ordered <- sub_matrix_log[ordered_genes, , drop = FALSE]
  
  pdf(file_names[i], width = 4, height = 8)  # 更瘦
  pheatmap(sub_matrix_log_ordered,
           scale = "none",
           cluster_rows = FALSE,
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           fontsize_row = 6,
           show_rownames = FALSE,
           show_colnames = TRUE,    # 显示列名
           angle_col = 0)           # 横向显示列名
  dev.off()
}


cat("✅ 四个更瘦、无标题的热图已完成。\n")





library(igraph)
library(dplyr)
library(scales)

# ============ 参数 ============
network_file <- "E:/1数据/≥2/c1_TF_TARGET_cleaned.txt"
network_name <- "Control_T1"

# ====== 1. 读取边表构建图 ======
edge_data <- read.table(network_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
g <- graph_from_data_frame(edge_data[, 1:2], directed = FALSE)

# ====== 2. 获取 Top10 TF（先过滤掉不在图中的 TF） ======
tf_in_graph <- intersect(all_tfs, V(g)$name)
deg <- degree(g, mode = "all")[tf_in_graph]
top10_tf <- names(sort(deg, decreasing = TRUE))[1:min(10, length(deg))]

# ====== 3. 获取 TF 及其靶基因构成子图 ======
neighbor_list <- lapply(top10_tf, function(tf) {
  neighbors(g, tf, mode = "all")$name
})
sub_nodes <- unique(c(top10_tf, unlist(neighbor_list)))
sub_nodes <- intersect(sub_nodes, V(g)$name)
sub_g <- induced_subgraph(g, vids = sub_nodes)

# ====== 4. 自定义 layout：TF 正十边形，靶基因椭圆分布 ======
# ---- TF 正十边形（半径为 r = 1.2） ----
N_tf <- length(top10_tf)
theta_tf <- seq(0, 2 * pi, length.out = N_tf + 1)[- (N_tf + 1)]
tf_coords <- cbind(cos(theta_tf), sin(theta_tf)) * 1.2

# ---- 靶基因椭圆分布 + 随机扰动 ----
target_nodes <- setdiff(V(sub_g)$name, top10_tf)
N_target <- length(target_nodes)
target_theta <- runif(N_target, 0, 2 * pi)
a <- 3.5; b <- 2.0
x <- a * cos(target_theta) + rnorm(N_target, 0, 0.4)
y <- b * sin(target_theta) + rnorm(N_target, 0, 0.4)
target_coords <- cbind(x, y)

# ---- 合并 layout：按 V(sub_g) 顺序排列 ----
layout_matrix <- matrix(NA, nrow = vcount(sub_g), ncol = 2)
layout_matrix[match(top10_tf, V(sub_g)$name), ] <- tf_coords
layout_matrix[match(target_nodes, V(sub_g)$name), ] <- target_coords

# ====== 5. 设置节点属性 ======
V(sub_g)$color <- ifelse(V(sub_g)$name %in% top10_tf, "red", "black")
V(sub_g)$size <- ifelse(V(sub_g)$name %in% top10_tf, 18, 2.5)  # 靶基因为小黑点
V(sub_g)$shape <- "circle"  # 所有节点都用合法形状
V(sub_g)$label <- ifelse(V(sub_g)$name %in% top10_tf, V(sub_g)$name, NA)
V(sub_g)$label.cex <- 0.8
V(sub_g)$label.color <- "black"
V(sub_g)$frame.color <- NA


# ====== 6. 绘图 ======
pdf(paste0("Top10_TF_polygon_layout_", network_name, ".pdf"), width = 10, height = 10)
plot(sub_g,
     layout = layout_matrix,
     vertex.label.family = "sans",
     edge.arrow.mode = 0,
     edge.color = "gray60",
     edge.width = 1,
     vertex.label.dist = 0.6,
     vertex.label.degree = pi/2,
     main = paste0("Top10 TF Network – ", network_name))
dev.off()



















###############################################################
#  Degree centrality — 仅转录因子 (TF) 四网络对比
#  —— 文件自带表头版 —— 2025‑05
###############################################################

library(igraph)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)

# ---------- ① TF 白名单 ----------
tf_list <- read.table("E:/1数据/Homo_TF_clean.txt",
                      header = FALSE,
                      stringsAsFactors = FALSE)
all_tfs <- unique(trimws(tf_list[[1]]))   # TF gene‑symbol 列表

# ---------- ② 网络文件 ----------
network_files <- c(
  "E:/1数据/≥2/c1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/c2_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s2_TF_TARGET_cleaned.txt"
)
network_names <- c("Control_T1", "Control_T2", "SCA_T1", "SCA_T2")

degree_centrality_list <- list()

# ---------- ③ 循环计算 ----------
for (i in seq_along(network_files)) {
  
  # 1) 读取边表（文件自带列名）
  edge_data <- read.table(
    network_files[i],
    header = TRUE,        # ✅ 文件已有表头
    sep    = "\t",
    stringsAsFactors = FALSE
  )
  
  # 2) 构建 igraph
  g <- graph_from_data_frame(edge_data[, 1:2], directed = TRUE)
  
  # 3) 仅取白名单且位于图中的 TF
  tf_in_graph <- intersect(all_tfs, V(g)$name)
  
  degree_values <- degree(g, mode = "all")[tf_in_graph]
  degree_values <- degree_values[!is.na(degree_values)]
  
  degree_df <- data.frame(
    Gene    = names(degree_values),
    Degree  = as.numeric(degree_values),
    Network = network_names[i],
    stringsAsFactors = FALSE
  )
  
  degree_centrality_list[[i]] <- degree_df
  
  # 4) 单网结果按 Degree 降序保存
  write.table(
    degree_df[order(-degree_df$Degree), ],
    file = paste0("degree_centrality_", network_names[i], ".txt"),
    sep  = "\t", quote = FALSE, row.names = FALSE
  )
}

# ---------- ④ 合并四网结果 ----------
all_degree_data <- bind_rows(degree_centrality_list)

# ---------- ⑤ 移除错误的 Gene 行（如 "TF"） ----------
all_degree_data <- all_degree_data %>%
  filter(Gene != "TF")

# ---------- ⑥ 画四网竖线分布图 ----------
pdf("degree_centrality_TF_only.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))

for (nw in network_names) {
  sub <- all_degree_data %>% 
    filter(Network == nw) %>% 
    arrange(desc(Degree))
  
  plot(
    sub$Degree,
    type = "h", col = "black",
    main = paste("TF Degree –", nw),
    xlab = "TFs (sorted)",
    ylab = "Degree",
    lwd  = 1
  )
}
dev.off()

# ---------- ⑦ 宽表转换 + 差值计算 ----------
degree_wide <- all_degree_data %>%
  pivot_wider(names_from = Network, values_from = Degree, values_fill = 0) %>%
  filter(Gene != "TF") %>%
  mutate(
    Diff_T1  = SCA_T1 - Control_T1,      # 疾病影响（运动前）
    Diff_T2  = SCA_T2 - Control_T2,      # 疾病影响（运动后）
    Diff_Con = Control_T2 - Control_T1,  # 运动影响（对照组）
    Diff_SCA = SCA_T2 - SCA_T1           # 运动影响（SCA组）
  )

# ---------- ⑧ 热图绘制 ----------
df_heatmap <- degree_wide[, c("Control_T1", "Control_T2", "SCA_T1", "SCA_T2")] %>%
  as.data.frame()
rownames(df_heatmap) <- degree_wide$Gene

pdf("TF_degree_heatmap_no_label.pdf", width = 10, height = 12)
pheatmap(df_heatmap,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         show_rownames = FALSE,     # ✅ 关闭行名标签
         main = "TF Degree Centrality Across Networks")
dev.off()


# ---------- ⑨ top10 差异 TF（运动对 SCA 的影响） ----------
top_diff_SCA <- degree_wide %>%
  arrange(desc(abs(Diff_SCA))) %>%
  head(50)

pdf("TF_degree_diff_SCA_T2_vs_T1_top50.pdf", width = 8, height = 6)
ggplot(top_diff_SCA, aes(x = reorder(Gene, Diff_SCA), y = Diff_SCA)) +
  geom_bar(stat = "identity", fill = "black") +
  coord_flip() +
  labs(title = "Top 50 TFs with Most Degree Change in SCA (T2 vs T1)",
       x = "TF", y = "Degree Difference")
dev.off()

# 筛选 top50 疾病状态（Control_T1 vs SCA_T1）差异 TF
top_diff_disease <- degree_wide %>%
  arrange(desc(abs(Diff_T1))) %>%  # 按差异绝对值降序
  head(50)# 绘制疾病状态差异柱状图
pdf("TF_degree_diff_disease_ControlT1_vs_SCAt1_top50.pdf", width = 8, height = 6)
ggplot(top_diff_disease, aes(x = reorder(Gene, Diff_T1), y = Diff_T1)) +
  geom_bar(stat = "identity", fill = "black") +
  coord_flip() +
  labs(
    title = "Top 50 TFs with Most Degree Change in Disease State (Control_T1 vs SCA_T1)",
    x = "TF", 
    y = "Degree Difference (SCA_T1 - Control_T1)"
  )
dev.off()

# ---------- ⑨-c 运动对 SCA 的影响可视化（SCA_T1 vs SCA_T2） ----------
df_sca_exercise <- degree_wide %>%
  select(Gene, SCA_T1, SCA_T2) %>%
  pivot_longer(cols = c(SCA_T1, SCA_T2),
               names_to = "Group", values_to = "Degree")

pdf("TF_degree_SCA_T1_vs_T2_Boxplot.pdf", width = 6, height = 6)
ggplot(df_sca_exercise, aes(x = Group, y = Degree, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5) +
  scale_fill_manual(values = c("SCA_T1" = "#e41a1c", "SCA_T2" = "#377eb8")) +
  labs(title = "TF Degree Centrality: SCA_T1 vs SCA_T2",
       x = "", y = "Degree Centrality") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
dev.off()

# ---------- ⑨-b 疾病状态比较可视化（Control_T1 vs SCA_T1） ----------
df_disease <- degree_wide %>%
  select(Gene, Control_T1, SCA_T1) %>%
  pivot_longer(cols = c(Control_T1, SCA_T1),
               names_to = "Group", values_to = "Degree")

pdf("TF_degree_disease_state_Comparison_Boxplot.pdf", width = 6, height = 6)
ggplot(df_disease, aes(x = Group, y = Degree, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5) +
  scale_fill_manual(values = c("Control_T1" = "#1f77b4", "SCA_T1" = "#d62728")) +
  labs(title = "TF Degree Centrality: Control vs SCA (T1)",
       x = "", y = "Degree Centrality") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
dev.off()


# ---------- ⑩ 非参数检验 ----------
cat("Wilcoxon 检验：疾病状态（Control_T1 vs SCA_T1）\n")
print(wilcox.test(degree_wide$Control_T1, degree_wide$SCA_T1, paired = FALSE))

cat("Wilcoxon 检验：运动对 SCA 影响（SCA_T1 vs SCA_T2）\n")
print(wilcox.test(degree_wide$SCA_T1, degree_wide$SCA_T2, paired = TRUE))

# ---------- ⑪ 保存完整表 ----------
write.table(degree_wide,
            "TF_degree_centrality_wide_table.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n✅ 全流程完成！输出文件包括：\n",
    "- 各网络单独中心性排序表\n",
    "- 合并宽表：TF_degree_centrality_wide_table.txt\n",
    "- 四网热图：TF_degree_heatmap.pdf\n",
    "- 运动对 SCA 影响最大的 TF 图：TF_degree_diff_SCA_T2_vs_T1_top10.pdf\n",
    "- Wilcoxon 检验结果打印在控制台\n")



###############################################################
#  Degree centrality — 仅转录因子 (TF) 四网络对比
#  —— 文件自带表头版 —— 2025‑05
#  —— 包含对比分析 ——
###############################################################

library(igraph)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)

# ---------- ① TF 白名单 ----------
tf_list <- read.table("E:/1数据/Homo_TF_clean.txt",
                      header = FALSE,
                      stringsAsFactors = FALSE)
all_tfs <- unique(trimws(tf_list[[1]]))   # TF gene‑symbol 列表

# ---------- ② 网络文件 ----------
network_files <- c(
  "E:/1数据/≥2/c1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/c2_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s2_TF_TARGET_cleaned.txt"
)
network_names <- c("Control_T1", "Control_T2", "SCA_T1", "SCA_T2")

degree_centrality_list <- list()
edge_data_list <- list()  # 存储各网络的边数据

# ---------- ③ 循环计算 ----------
for (i in seq_along(network_files)) {
  
  # 1) 读取边表（文件自带列名）
  edge_data <- read.table(
    network_files[i],
    header = TRUE,        # ✅ 文件已有表头
    sep    = "\t",
    stringsAsFactors = FALSE
  )
  edge_data_list[[i]] <- edge_data
  
  # 2) 构建 igraph
  g <- graph_from_data_frame(edge_data[, 1:2], directed = TRUE)
  
  # 3) 仅取白名单且位于图中的 TF
  tf_in_graph <- intersect(all_tfs, V(g)$name)
  
  degree_values <- degree(g, mode = "all")[tf_in_graph]
  degree_values <- degree_values[!is.na(degree_values)]
  
  degree_df <- data.frame(
    Gene    = names(degree_values),
    Degree  = as.numeric(degree_values),
    Network = network_names[i],
    stringsAsFactors = FALSE
  )
  
  degree_centrality_list[[i]] <- degree_df
  
  # 4) 单网结果按 Degree 降序保存
  write.table(
    degree_df[order(-degree_df$Degree), ],
    file = paste0("degree_centrality_", network_names[i], ".txt"),
    sep  = "\t", quote = FALSE, row.names = FALSE
  )
}

# ---------- ④ 合并四网结果 ----------
all_degree_data <- bind_rows(degree_centrality_list)

# ---------- ⑤ 绘制分布竖线图 ----------
pdf("degree_centrality_TF_only.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))

for (nw in network_names) {
  
  sub <- all_degree_data %>% 
    filter(Network == nw) %>% 
    arrange(desc(Degree))
  
  plot(
    sub$Degree,
    type = "h", col = "black",
    main = paste("TF Degree –", nw),
    xlab = "TFs (sorted)",
    ylab = "Degree",
    lwd  = 1
  )
}
dev.off()

# ---------- ⑥ 转换为宽表并计算差异 ----------
degree_wide <- all_degree_data %>%
  pivot_wider(names_from = Network, values_from = Degree, values_fill = 0) %>%
  mutate(
    # 疾病影响
    Disease_Effect_T1 = SCA_T1 - Control_T1,  # 运动前疾病影响
    Disease_Effect_T2 = SCA_T2 - Control_T2,  # 运动后疾病影响
    
    # 运动影响
    Exercise_Control = Control_T2 - Control_T1,  # 对照组运动影响
    Exercise_SCA = SCA_T2 - SCA_T1,             # 患病组运动影响
    
    # 交互效应
    Interaction = (SCA_T2 - SCA_T1) - (Control_T2 - Control_T1)
  )

# ---------- ⑦ 保存完整宽表 ----------
write.table(degree_wide,
            "TF_degree_centrality_wide_table.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- ⑧ 疾病状态影响分析 ----------
# 筛选疾病影响（运动前）Top50差异的TF
top_diff_disease_t1 <- degree_wide %>%
  arrange(desc(abs(Disease_Effect_T1))) %>%
  head(50)

# 绘制疾病状态影响柱状图
pdf("diff_disease_t1_top50.pdf", width = 8, height = 6)
ggplot(top_diff_disease_t1, aes(x = reorder(Gene, Disease_Effect_T1), y = Disease_Effect_T1)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  coord_flip() +
  labs(
    title = "Top 50 TFs: SCA_T1 vs Control_T1 (Disease Effect at T1)",
    x = "Transcription Factor (TF)",
    y = "Degree Difference (SCA_T1 - Control_T1)"
  )
dev.off()

# 疾病状态影响的统计检验
wilcox_disease_t1 <- wilcox.test(degree_wide$SCA_T1, degree_wide$Control_T1, paired = FALSE)
cat("Wilcoxon 检验结果（SCA_T1 vs Control_T1）：\n")
print(wilcox_disease_t1)

# ---------- ⑨ 运动对患病组影响分析 ----------
# 筛选运动对患病组影响Top50差异的TF
top_diff_exercise_sca <- degree_wide %>%
  arrange(desc(abs(Exercise_SCA))) %>%
  head(50)

# 绘制运动对患病组影响柱状图
pdf("diff_exercise_sca_top50.pdf", width = 8, height = 6)
ggplot(top_diff_exercise_sca, aes(x = reorder(Gene, Exercise_SCA), y = Exercise_SCA)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 50 TFs: SCA_T2 vs SCA_T1 (Exercise Effect on SCA)",
    x = "Transcription Factor (TF)",
    y = "Degree Difference (SCA_T2 - SCA_T1)"
  )
dev.off()

# 运动对患病组影响的统计检验
wilcox_exercise_sca <- wilcox.test(degree_wide$SCA_T1, degree_wide$SCA_T2, paired = TRUE)
cat("Wilcoxon 检验结果（SCA_T2 vs SCA_T1，运动对患病组的影响）：\n")
print(wilcox_exercise_sca)

# ---------- ⑩ 交互效应分析 ----------
# 筛选交互效应Top50差异的TF
top_interaction_tfs <- degree_wide %>%
  arrange(desc(abs(Interaction))) %>%
  head(50)

# 绘制交互效应柱状图
pdf("interaction_effect_top50.pdf", width = 8, height = 6)
ggplot(top_interaction_tfs, aes(x = reorder(Gene, Interaction), y = Interaction)) +
  geom_bar(stat = "identity", fill = "purple") +
  coord_flip() +
  labs(
    title = "Top 50 TFs with Significant Interaction Effect (Exercise*Disease)",
    x = "Transcription Factor (TF)",
    y = "Interaction Effect"
  )
dev.off()
library(tibble) 
# ---------- ⑪ 多组差异热图 ----------
# 提取用于热图的差异列
heatmap_data <- degree_wide %>%
  select(Gene, Control_T1, Control_T2, SCA_T1, SCA_T2, 
         Disease_Effect_T1, Disease_Effect_T2, 
         Exercise_Control, Exercise_SCA, Interaction) %>%
  column_to_rownames(var = "Gene")

# 绘制热图
pdf("diff_heatmap.pdf", width = 12, height = 10)
pheatmap(
  heatmap_data, 
  scale = "row",  # 按行标准化
  clustering_distance_rows = "euclidean", 
  clustering_method = "ward.D2",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "TF Degree Differences Across Comparisons",
  show_rownames = FALSE  # 基因太多时可关闭行名显示
)
dev.off()

# ---------- ⑫ 计算四组网络的整体统计量 ----------
network_stats <- data.frame(
  Network = network_names,
  Nodes = numeric(4),
  Edges = numeric(4),
  Density = numeric(4),
  AvgDegree = numeric(4)
)

for(i in 1:4) {
  g <- graph_from_data_frame(edge_data_list[[i]][,1:2], directed=TRUE)
  network_stats$Nodes[i] <- vcount(g)
  network_stats$Edges[i] <- ecount(g)
  network_stats$Density[i] <- edge_density(g)
  network_stats$AvgDegree[i] <- mean(degree(g))
}

# 可视化网络统计量对比
pdf("network_stats_comparison.pdf", width = 10, height = 8)
ggplot(network_stats %>% 
         tidyr::pivot_longer(cols = -Network, names_to = "Metric", values_to = "Value"), 
       aes(x = Network, y = Value, fill = Network)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Metric, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# ---------- ⑬ 输出结果摘要 ----------
cat("\n✅ 全流程完成！输出文件包括：\n",
    "- 各网络单独中心性排序表：degree_centrality_<Network>.txt\n",
    "- 合并宽表：TF_degree_centrality_wide_table.txt\n",
    "- 四网度中心性分布图：degree_centrality_TF_only.pdf\n",
    "- 疾病状态影响分析：diff_disease_t1_top50.pdf\n",
    "- 运动对患病组影响分析：diff_exercise_sca_top50.pdf\n",
    "- 交互效应分析：interaction_effect_top50.pdf\n",
    "- 差异热图：diff_heatmap.pdf\n",
    "- 网络统计量对比图：network_stats_comparison.pdf\n",
    "- Wilcoxon 检验结果打印在控制台\n")












library(igraph)
library(dplyr)

# ---------- ① TF 白名单 ----------
tf_list <- read.table("E:/1数据/Homo_TF_clean.txt",
                      header = FALSE,
                      stringsAsFactors = FALSE)
all_tfs <- unique(trimws(tf_list[[1]]))   # TF gene‑symbol 列表


# ---------- ② 网络文件 ----------
network_files <- c(
  "E:/1数据/≥2/c1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/c2_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s2_TF_TARGET_cleaned.txt"
)
network_names <- c("Control_T1", "Control_T2", "SCA_T1", "SCA_T2")

degree_centrality_list <- list()

# ---------- ③ 循环计算每个网络 ----------
for (i in seq_along(network_files)) {
  
  # 1) 读取边表
  edge_data <- read.table(
    network_files[i],
    header = TRUE,
    sep    = "\t",
    stringsAsFactors = FALSE
  )
  
  # 2) 构建 igraph
  g <- graph_from_data_frame(edge_data[, 1:2], directed = TRUE)
  
  # 3) 仅取白名单且位于图中的 TF
  tf_in_graph <- intersect(all_tfs, V(g)$name)
  
  degree_values <- degree(g, mode = "all")[tf_in_graph]
  degree_values <- degree_values[!is.na(degree_values)]
  
  degree_df <- data.frame(
    Gene    = names(degree_values),
    Degree  = as.numeric(degree_values),
    Network = network_names[i],
    stringsAsFactors = FALSE
  )
  
  degree_centrality_list[[i]] <- degree_df
  
  # 4) 单网结果保存
  write.table(
    degree_df[order(-degree_df$Degree), ],
    file = paste0("degree_centrality_", network_names[i], ".txt"),
    sep  = "\t", quote = FALSE, row.names = FALSE
  )
}


# ---------- ④ 合并四网结果 ----------
all_degree_data <- bind_rows(degree_centrality_list)


# ---------- ⑤ 差异分析：前200个TF ----------
get_top200_diff <- function(df1, df2) {
  merged <- merge(df1, df2, by = "Gene", suffixes = c(".x", ".y"))
  merged <- merged %>%
    mutate(Diff = abs(Degree.x - Degree.y)) %>%
    arrange(desc(Diff))
  top200 <- head(merged, 200)
  return(top200)
}


# ---------- ⑥ 取不同组合的子集 ----------
control_t1 <- all_degree_data %>% filter(Network == "Control_T1")
control_t2 <- all_degree_data %>% filter(Network == "Control_T2")
sca_t1     <- all_degree_data %>% filter(Network == "SCA_T1")
sca_t2     <- all_degree_data %>% filter(Network == "SCA_T2")


# ---------- ⑦ 差异前200 TF计算 ----------
diff_c1_c2 <- get_top200_diff(control_t1, control_t2)
diff_s1_s2 <- get_top200_diff(sca_t1, sca_t2)
diff_c1_s1 <- get_top200_diff(control_t1, sca_t1)
diff_c2_s2 <- get_top200_diff(control_t2, sca_t2)


# ---------- ⑧ 保存结果 ----------
write.table(diff_c1_c2, file = "degree_diff_Control_T1_vs_Control_T2_top200.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diff_s1_s2, file = "degree_diff_SCA_T1_vs_SCA_T2_top200.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diff_c1_s1, file = "degree_diff_Control_T1_vs_SCA_T1_top200.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(diff_c2_s2, file = "degree_diff_Control_T2_vs_SCA_T2_top200.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# ---------- ⑨ 控制台提醒 ----------
cat("\n✅ 计算完成！已输出以下文件：\n",
    "- Control_T1 vs Control_T2：degree_diff_Control_T1_vs_Control_T2_top200.txt\n",
    "- SCA_T1 vs SCA_T2        ：degree_diff_SCA_T1_vs_SCA_T2_top200.txt\n",
    "- Control_T1 vs SCA_T1    ：degree_diff_Control_T1_vs_SCA_T1_top200.txt\n",
    "- Control_T2 vs SCA_T2    ：degree_diff_Control_T2_vs_SCA_T2_top200.txt\n")





library(igraph)
library(dplyr)
library(tidyr)

# ---------- ① 读取 TF 白名单 ----------
tf_list <- read.table("E:/1数据/Homo_TF_clean.txt",
                      header = FALSE,
                      stringsAsFactors = FALSE)
all_tfs <- unique(trimws(tf_list[[1]]))

# ---------- ② 读取四个网络 ----------
network_files <- c(
  "E:/1数据/≥2/c1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/c2_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s2_TF_TARGET_cleaned.txt"
)
network_names <- c("Control_T1", "Control_T2", "SCA_T1", "SCA_T2")

degree_centrality_list <- list()

for (i in seq_along(network_files)) {
  edge_data <- read.table(network_files[i],
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE)
  
  g <- graph_from_data_frame(edge_data[, 1:2], directed = TRUE)
  
  tf_in_graph <- intersect(all_tfs, V(g)$name)
  
  degree_values <- degree(g, mode = "all")[tf_in_graph]
  degree_values <- degree_values[!is.na(degree_values)]
  
  degree_df <- data.frame(
    Gene    = names(degree_values),
    Degree  = as.numeric(degree_values),
    Network = network_names[i],
    stringsAsFactors = FALSE
  )
  
  degree_centrality_list[[i]] <- degree_df
}

all_degree_data <- bind_rows(degree_centrality_list)

# ---------- ③ 找疾病相关变化最大的 200 个 TF ----------
control_t1 <- all_degree_data %>% filter(Network == "Control_T1")
sca_t1     <- all_degree_data %>% filter(Network == "SCA_T1")

merged <- merge(control_t1, sca_t1, by = "Gene", suffixes = c(".Control_T1", ".SCA_T1"))
merged <- merged %>%
  mutate(Diff = abs(Degree.Control_T1 - Degree.SCA_T1)) %>%
  arrange(desc(Diff))

top200_tfs <- head(merged$Gene, 200)

# ---------- ④ 提取这 200 个 TF 在四个组别的 Degree ----------
degree_wide <- all_degree_data %>%
  filter(Gene %in% top200_tfs) %>%
  select(Gene, Network, Degree) %>%
  pivot_wider(names_from = Network, values_from = Degree)

# ---------- ⑤ 判断是否恢复（SCA_T2 是否更接近 Control_T1，且减少幅度 > 50%） ----------
degree_wide <- degree_wide %>%
  mutate(
    Dist_before = abs(SCA_T1 - Control_T1),
    Dist_after  = abs(SCA_T2 - Control_T1),
    Recovery = ifelse(Dist_after < Dist_before, "Yes", "No"),
    Effective_Recovery = ifelse((Dist_before - Dist_after) / Dist_before > 0.5, "Yes", "No")
  )

# ---------- ⑥ 保存为 CSV ----------
write.csv(degree_wide, file = "TF_degree_trend_SCA_recovery.csv", row.names = FALSE)

cat("\n✅ 已输出文件：TF_degree_trend_SCA_recovery.csv\n")




library(data.table)
library(igraph)
library(dplyr)
library(tidyr)

## =============== 1️⃣ Poisson 方法筛选关键 TF ====================

edge_files <- c(
  c1 = "E:/1数据/≥2/evaluation_results-c1/two_or_more_overlap_c1_TF_TARGET.txt",
  c2 = "E:/1数据/≥2/evaluation_results-c2/two_or_more_overlap_c2_TF_TARGET.txt",
  s1 = "E:/1数据/≥2/evaluation_results-s1/two_or_more_overlap_s1_TF_TARGET.txt",
  s2 = "E:/1数据/≥2/evaluation_results-s2/two_or_more_overlap_s2_TF_TARGET.txt"
)

label_map <- c(
  c1 = "对照组运动前",
  c2 = "对照组运动后",
  s1 = "患病组运动前",
  s2 = "患病组运动后"
)

out_dir <- "E:/1数据/泊松分布筛选关键转录因子"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

alpha     <- 0.05
tail_prob <- 0.999

tf_list <- fread("E:/1数据/Homo_TF_clean.txt", header = TRUE, sep = "\t", colClasses = "character")[, TF]

poisson_tf_all <- list()

for (tag in names(edge_files)) {
  ef <- edge_files[[tag]]
  edges <- fread(ef, header = TRUE, sep = "\t", colClasses = "character")
  setnames(edges, 1:2, c("TF", "Target"))
  
  g <- graph_from_data_frame(edges, directed = TRUE)
  
  n <- vcount(g)
  m <- ecount(g)
  lambda <- m / n
  
  deg_out <- degree(g, mode = "out")
  deg_tf  <- deg_out[names(deg_out) %in% tf_list]
  
  p_raw <- 1 - ppois(deg_tf - 1, lambda = lambda)
  q_val <- p.adjust(p_raw, method = "BH")
  
  deg_cut <- qpois(tail_prob, lambda = lambda)
  sig_idx <- which(q_val < alpha & deg_tf >= deg_cut)
  
  net_tag <- sub("\\.txt$", "", basename(ef))
  
  res_all <- data.frame(
    Gene    = names(deg_tf),
    Degree  = as.integer(deg_tf),
    p_value = signif(p_raw,  3),
    q_value = signif(q_val,  3),
    stringsAsFactors = FALSE
  )
  res_sig <- res_all[sig_idx, ]
  
  fwrite(res_all, file = file.path(out_dir, paste0(net_tag, "_poisson_all.txt")), sep  = "\t", quote = FALSE)
  fwrite(res_sig, file = file.path(out_dir, paste0(net_tag, "_poisson_TF.txt")), sep  = "\t", quote = FALSE)
  
  cat(sprintf("\n[%s] λ = %.2f, 度阈 = %d, 关键 TF = %d\n", label_map[tag], lambda, deg_cut, nrow(res_sig)))
  
  poisson_tf_all[[tag]] <- res_sig$Gene
}

# 合并四个网络的 Poisson 筛选结果
important_tf <- unique(unlist(poisson_tf_all))

cat("\n🎯 筛选后保留的转录因子数量：", length(important_tf), "\n")


## =============== 2️⃣ 这批 TF 再做恢复分析 ====================

network_files <- c(
  "E:/1数据/≥2/c1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/c2_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s2_TF_TARGET_cleaned.txt"
)
network_names <- c("Control_T1", "Control_T2", "SCA_T1", "SCA_T2")

degree_centrality_list <- list()

for (i in seq_along(network_files)) {
  edge_data <- fread(network_files[i], header = TRUE, sep = "\t", colClasses = "character")
  
  g <- graph_from_data_frame(edge_data[, 1:2], directed = TRUE)
  
  tf_in_graph <- intersect(important_tf, V(g)$name)
  
  degree_values <- degree(g, mode = "all")[tf_in_graph]
  degree_values <- degree_values[!is.na(degree_values)]
  
  degree_df <- data.frame(
    Gene    = names(degree_values),
    Degree  = as.numeric(degree_values),
    Network = network_names[i],
    stringsAsFactors = FALSE
  )
  
  degree_centrality_list[[i]] <- degree_df
}

all_degree_data <- bind_rows(degree_centrality_list)

# ----------------- 转为宽表
degree_wide <- all_degree_data %>%
  select(Gene, Network, Degree) %>%
  pivot_wider(names_from = Network, values_from = Degree)


# ----------------- 恢复判断，50% 阈值
degree_wide <- degree_wide %>%
  mutate(
    Dist_before = abs(SCA_T1 - Control_T1),
    Dist_after  = abs(SCA_T2 - Control_T1),
    Recovery = ifelse(Dist_after < Dist_before, "Yes", "No"),
    Effective_Recovery = ifelse((Dist_before - Dist_after) / Dist_before > 0.5, "Yes", "No")
  )


# ----------------- 保存
write.csv(degree_wide, file = file.path(out_dir, "TF_degree_trend_SCA_recovery.csv"), row.names = FALSE)
cat("\n✅ 已输出最终文件：", file.path(out_dir, "TF_degree_trend_SCA_recovery.csv"), "\n")
