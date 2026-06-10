# -----------------------------------------------------------------------------
# 完整脚本：计算所有方法组合的 AUROC/AUPR (GSE117221数据集)
# 包括：
#   1. 单个方法（3个）
#   2. 两两组合交集（3个）
#   3. 三个方法交集（1个）
#   4. 两个及以上方法（并集/平均，1个）
# 总计：8个网络
# -----------------------------------------------------------------------------

# —— 0. 清理环境 & 加载包 ——
rm(list = ls()); gc()

if (!requireNamespace("data.table", quietly=TRUE))   install.packages("data.table")
if (!requireNamespace("dplyr", quietly=TRUE))        install.packages("dplyr")
if (!requireNamespace("pROC", quietly=TRUE))         install.packages("pROC")
if (!requireNamespace("PRROC", quietly=TRUE))        install.packages("PRROC")
if (!requireNamespace("ggplot2", quietly=TRUE))      install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly=TRUE))        install.packages("tidyr")

library(data.table)
library(dplyr)
library(pROC)
library(PRROC)
library(ggplot2)
library(tidyr)

# —— 1. 定义清洗函数 ——
clean_name <- function(x) {
  x %>% trimws() %>% tolower() %>% gsub("[^a-z0-9]", "", .)
}

# —— 2. 输入：网络文件（3 个 top10pct.tsv）——
net_paths <- list(
  deepfgrn = "E:/SCD/其他数据/GSE117221/grn/deepfgrn_H_top10pct.tsv",
  `3dcema` = "E:/SCD/其他数据/GSE117221/grn/3dcema_H_top10pct.tsv",
  deepsem  = "E:/SCD/其他数据/GSE117221/grn/deepsem_H_top10pct.tsv"
)

# 输出根目录
out_root <- "E:/SCD/其他数据/GSE117221/AUROC"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# —— 3. 读取三个方法预测文件 & 构建 dt_list ——
dt_list <- lapply(names(net_paths), function(meth) {
  dt <- fread(net_paths[[meth]], sep = "\t", header = TRUE)
  
  # 兼容不同的权重列名
  if (!("EdgeWeight" %in% names(dt))) {
    cand <- intersect(names(dt), c("weight", "Weight", "score", "Score", "edgeweight", "Edge_Weight"))
    if (length(cand) >= 1) setnames(dt, cand[1], "EdgeWeight")
  }
  
  dt[, TF     := clean_name(TF)]
  dt[, Target := clean_name(Target)]
  dt[, pair   := paste0(TF, "_", Target)]
  
  # 若同一方法内同一 pair 多行，取均值聚合
  dt <- dt[, .(
    EdgeWeight = mean(EdgeWeight, na.rm = TRUE)
  ), by = .(TF, Target, pair)]
  
  # 添加方法名称列
  dt[, Method := meth]
  dt
})
names(dt_list) <- names(net_paths)

# —— 4. 读取金标准数据 ——
# 注意：GSE117221的金标准文件路径（请根据实际情况修改）
chip_file <- "E:/SCD/其他数据/gse133181_单细胞/GSE133181_RAW/chip/Integrated_ChIP_Top200.txt"

# 检查文件是否存在
if (!file.exists(chip_file)) {
  stop(paste("金标准文件不存在:", chip_file, 
             "\n请更新正确的金标准文件路径！"))
}

cat("读取金标准文件:", chip_file, "\n")
chip_dt <- fread(chip_file, sep = "\t", header = TRUE)

# 检查列名
if (!all(c("TF", "Target") %in% names(chip_dt))) {
  stop("ChIP文件必须包含 'TF' 和 'Target' 两列")
}

# 清洗并生成 pair
chip_dt[, TF_clean := clean_name(TF)]
chip_dt[, Target_clean := clean_name(Target)]
chip_dt[, pair := paste0(TF_clean, "_", Target_clean)]

# 获取所有金标准pair
all_chip_pairs <- unique(chip_dt$pair)

cat("金标准总边数:", nrow(chip_dt), "\n")
cat("去重后金标准边数:", length(all_chip_pairs), "\n")
cat("涉及TF数量:", length(unique(chip_dt$TF_clean)), "\n\n")

# —— 5. 定义评估函数 ——
evaluate_network <- function(network_dt, network_name, gold_pairs, out_dir) {
  cat("\n评估网络:", network_name, "\n")
  cat("  边数量:", nrow(network_dt), "\n")
  
  # 准备数据
  df <- copy(network_dt)[
    , `:=`(
      label = as.integer(pair %in% gold_pairs),
      score = EdgeWeight
    )
  ]
  
  # 统计
  pos_count <- sum(df$label)
  neg_count <- sum(1 - df$label)
  cat("  正样本数:", pos_count, "(", round(pos_count/nrow(df)*100, 2), "%)\n")
  cat("  负样本数:", neg_count, "(", round(neg_count/nrow(df)*100, 2), "%)\n")
  
  # 检查是否可以计算
  if (length(unique(df$label)) < 2) {
    cat("  ⚠️ 单一标签，无法计算 ROC/PR\n")
    return(NULL)
  }
  
  # 计算 AUROC
  roc_obj <- roc(df$label, df$score, quiet = TRUE, direction = "<")
  auroc <- auc(roc_obj)
  
  # 计算 AUPR
  pr_obj <- tryCatch({
    pr.curve(
      scores.class0 = df$score[df$label == 1],
      scores.class1 = df$score[df$label == 0],
      curve = TRUE
    )
  }, error = function(e) {
    cat("  警告: PR计算出错，尝试替代方法\n")
    return(NULL)
  })
  
  if (!is.null(pr_obj) && is.list(pr_obj) && !is.null(pr_obj$auc.integral)) {
    aupr <- pr_obj$auc.integral
  } else {
    # 使用pROC计算近似AUPR
    pr_roc <- roc(df$label, df$score, quiet = TRUE)
    coords_df <- coords(pr_roc, "all", ret = c("precision", "recall"))
    coords_df <- na.omit(coords_df)
    if (nrow(coords_df) > 1) {
      aupr <- sum(diff(coords_df$recall) * 
                    (head(coords_df$precision, -1) + tail(coords_df$precision, -1)) / 2)
    } else {
      aupr <- NA
    }
  }
  
  # 创建网络专属目录
  net_dir <- file.path(out_dir, gsub("[^a-zA-Z0-9]", "_", network_name))
  dir.create(net_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 绘制 ROC 曲线
  png(file.path(net_dir, "ROC_curve.png"), 800, 600)
  plot(roc_obj, 
       main = sprintf("%s\nAUROC = %.4f", network_name, auroc),
       col = "blue", lwd = 2, legacy.axes = TRUE)
  abline(a = 0, b = 1, lty = 2, col = "gray")
  legend("bottomright", legend = c("ROC曲线", "随机分类器"), 
         col = c("blue", "gray"), lty = c(1, 2), lwd = c(2, 1))
  dev.off()
  
  # 绘制 PR 曲线（如果成功）
  if (!is.null(pr_obj) && is.list(pr_obj) && !is.null(pr_obj$curve)) {
    png(file.path(net_dir, "PR_curve.png"), 800, 600)
    plot(pr_obj$curve[, 1], pr_obj$curve[, 2],
         type = "l", lwd = 2, col = "red",
         xlab = "Recall (灵敏度)", ylab = "Precision (精确率)",
         main = sprintf("%s\nAUPR = %.4f", network_name, aupr))
    grid()
    dev.off()
  }
  
  # 保存结果
  result <- data.table(
    Network = network_name,
    Type = ifelse(grepl("_vs_", network_name), "Intersection",
                  ifelse(grepl("TwoOrMore", network_name), "Union_2plus",
                         ifelse(grepl("AllThree", network_name), "Intersection_3", "Single"))),
    N_Edges = nrow(df),
    Pos_Edges = pos_count,
    Neg_Edges = neg_count,
    Pos_Ratio = round(pos_count/nrow(df)*100, 4),
    AUROC = round(as.numeric(auroc), 4),
    AUPR = round(aupr, 4)
  )
  
  # 保存网络文件
  fwrite(network_dt, file.path(net_dir, "network_edges.tsv"), sep = "\t", quote = FALSE)
  
  # 保存结果
  fwrite(result, file.path(net_dir, "evaluation_results.tsv"), sep = "\t", quote = FALSE)
  
  cat("  ✓ AUROC =", round(auroc, 4), ", AUPR =", round(aupr, 4), "\n")
  
  return(result)
}

# —— 6. 构建所有网络 ——

cat("\n", rep("=", 80), "\n", sep = "")
cat("开始构建所有网络组合\n")
cat(rep("=", 80), "\n", sep = "")

all_networks <- list()
method_names <- names(net_paths)

# 6.1 单个方法网络
cat("\n1. 构建单个方法网络...\n")
for (meth in method_names) {
  network_dt <- dt_list[[meth]][, .(TF, Target, pair, EdgeWeight)]
  all_networks[[meth]] <- network_dt
  cat("  ✓", meth, ":", nrow(network_dt), "条边\n")
}

# 6.2 两两组合交集
cat("\n2. 构建两两组合交集网络...\n")
combinations <- combn(method_names, 2, simplify = FALSE)
for (comb in combinations) {
  comb_name <- paste(comb, collapse = "_vs_")
  
  # 获取两个方法的pair集合
  pairs1 <- dt_list[[comb[1]]]$pair
  pairs2 <- dt_list[[comb[2]]]$pair
  
  # 计算交集
  intersect_pairs <- intersect(pairs1, pairs2)
  
  if (length(intersect_pairs) > 0) {
    # 提取交集边的信息，权重取两个方法的平均值
    dt1_sub <- dt_list[[comb[1]]][pair %in% intersect_pairs, .(pair, TF, Target, Weight1 = EdgeWeight)]
    dt2_sub <- dt_list[[comb[2]]][pair %in% intersect_pairs, .(pair, Weight2 = EdgeWeight)]
    
    network_dt <- merge(dt1_sub, dt2_sub, by = "pair")
    network_dt[, EdgeWeight := (Weight1 + Weight2) / 2]
    network_dt <- network_dt[, .(TF, Target, pair, EdgeWeight)]
    
    all_networks[[comb_name]] <- network_dt
    cat("  ✓", comb_name, ":", nrow(network_dt), "条边\n")
  } else {
    cat("  ✗", comb_name, ": 交集为空\n")
  }
}

# 6.3 三个方法交集
cat("\n3. 构建三个方法交集网络...\n")
pairs_all <- lapply(dt_list, function(x) x$pair)
common_pairs <- Reduce(intersect, pairs_all)

if (length(common_pairs) > 0) {
  # 提取三个方法的权重并取平均
  network_dt <- data.table(pair = common_pairs)
  for (meth in method_names) {
    dt_sub <- dt_list[[meth]][pair %in% common_pairs, .(pair, Weight = EdgeWeight)]
    setnames(dt_sub, "Weight", paste0("Weight_", meth))
    network_dt <- merge(network_dt, dt_sub, by = "pair", all.x = TRUE)
  }
  
  # 计算平均权重
  weight_cols <- paste0("Weight_", method_names)
  network_dt[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = weight_cols]
  
  # 添加TF和Target信息
  network_dt <- merge(network_dt, dt_list[[1]][, .(pair, TF, Target)], by = "pair")
  network_dt <- network_dt[, .(TF, Target, pair, EdgeWeight)]
  
  all_networks[["AllThree_Intersection"]] <- network_dt
  cat("  ✓ AllThree_Intersection:", nrow(network_dt), "条边\n")
} else {
  cat("  ✗ AllThree_Intersection: 交集为空\n")
}

# 6.4 两个及以上方法（保留至少两个方法预测的边，权重取平均）
cat("\n4. 构建两个及以上方法网络...\n")

# 合并所有方法的权重
merged_dt <- data.table(pair = unique(unlist(lapply(dt_list, function(x) x$pair))))

for (meth in method_names) {
  dt_sub <- dt_list[[meth]][, .(pair, Weight = EdgeWeight)]
  setnames(dt_sub, "Weight", paste0("Weight_", meth))
  merged_dt <- merge(merged_dt, dt_sub, by = "pair", all.x = TRUE)
}

# 计算每个pair出现在几个方法中
weight_cols <- paste0("Weight_", method_names)
merged_dt[, count_present := rowSums(!is.na(.SD)), .SDcols = weight_cols]

# 筛选出现次数 >= 2 的边
two_or_more <- merged_dt[count_present >= 2]

if (nrow(two_or_more) > 0) {
  # 计算平均权重
  two_or_more[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = weight_cols]
  
  # 添加TF和Target信息
  two_or_more <- merge(two_or_more, dt_list[[1]][, .(pair, TF, Target)], by = "pair")
  network_dt <- two_or_more[, .(TF, Target, pair, EdgeWeight)]
  
  all_networks[["TwoOrMore_Union"]] <- network_dt
  cat("  ✓ TwoOrMore_Union:", nrow(network_dt), "条边\n")
} else {
  cat("  ✗ TwoOrMore_Union: 没有边\n")
}

# —— 7. 对所有网络进行评估 ——
cat("\n", rep("=", 80), "\n", sep = "")
cat("开始评估所有网络\n")
cat(rep("=", 80), "\n", sep = "")

all_results <- list()

for (net_name in names(all_networks)) {
  result <- evaluate_network(all_networks[[net_name]], net_name, all_chip_pairs, out_root)
  if (!is.null(result)) {
    all_results[[net_name]] <- result
  }
}

# —— 8. 汇总所有结果 ——
cat("\n", rep("=", 80), "\n", sep = "")
cat("汇总结果\n")
cat(rep("=", 80), "\n", sep = "")

if (length(all_results) > 0) {
  summary_all <- rbindlist(all_results)
  
  # 按AUROC排序
  setorder(summary_all, -AUROC)
  
  # 保存汇总结果
  fwrite(summary_all, file.path(out_root, "ALL_NETWORKS_SUMMARY.tsv"), 
         sep = "\t", quote = FALSE)
  
  # 打印汇总表格
  cat("\n网络性能汇总表:\n")
  print(summary_all)
  
  # —— 9. 可视化 ——
  
  # 9.1 AUROC和AUPR柱状图
  summary_plot <- summary_all %>%
    pivot_longer(cols = c(AUROC, AUPR), names_to = "Metric", values_to = "Value")
  
  p1 <- ggplot(summary_plot, aes(x = reorder(Network, Value), y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "GSE117221: 所有网络性能比较",
         x = "网络", y = "分数") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "performance_comparison.png"), p1, width = 12, height = 7)
  
  # 9.2 按类型分组比较
  p2 <- ggplot(summary_all, aes(x = Type, y = AUROC, fill = Type)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    theme_minimal() +
    labs(title = "GSE117221: 不同网络类型的AUROC比较",
         x = "网络类型", y = "AUROC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "aurocs_by_type.png"), p2, width = 8, height = 6)
  
  # 9.3 AUROC vs 网络大小的关系
  p3 <- ggplot(summary_all, aes(x = N_Edges, y = AUROC, color = Type, label = Network)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
    theme_minimal() +
    labs(title = "GSE117221: AUROC vs 网络大小",
         x = "网络边数", y = "AUROC") +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "aurocs_vs_network_size.png"), p3, width = 10, height = 7)
  
  # 9.4 热图展示各方法的重叠情况（可选）
  # 创建重叠矩阵
  overlap_matrix <- matrix(0, nrow = length(method_names), ncol = length(method_names))
  rownames(overlap_matrix) <- method_names
  colnames(overlap_matrix) <- method_names
  
  for (i in 1:length(method_names)) {
    for (j in 1:length(method_names)) {
      if (i <= j) {
        pairs_i <- dt_list[[method_names[i]]]$pair
        pairs_j <- dt_list[[method_names[j]]]$pair
        overlap <- length(intersect(pairs_i, pairs_j))
        total <- length(union(pairs_i, pairs_j))
        overlap_matrix[i, j] <- overlap / total * 100
        overlap_matrix[j, i] <- overlap / total * 100
      }
    }
  }
  
  # 保存重叠矩阵
  fwrite(as.data.table(overlap_matrix, keep.rownames = "Method"), 
         file.path(out_root, "method_overlap_matrix.tsv"), sep = "\t")
  
  # 9.5 打印最佳网络信息
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("最佳网络 (按AUROC):\n")
  best_network <- summary_all[which.max(AUROC)]
  print(best_network)
  
  cat("\n最佳网络 (按AUPR):\n")
  best_pr_network <- summary_all[which.max(AUPR)]
  print(best_pr_network)
  
  cat("\n网络大小与性能:\n")
  cat("  最小网络:", summary_all[which.min(N_Edges), .(Network, N_Edges, AUROC)] %>% 
        mutate(across(where(is.numeric), ~round(., 4))), "\n")
  cat("  最大网络:", summary_all[which.max(N_Edges), .(Network, N_Edges, AUROC)] %>% 
        mutate(across(where(is.numeric), ~round(., 4))), "\n")
  
  cat("\n✔️ 所有评估完成！\n")
  cat("结果保存在:", out_root, "\n")
  cat("   - 汇总文件: ALL_NETWORKS_SUMMARY.tsv\n")
  cat("   - 各网络详细结果在对应子目录中\n")
  cat("   - 可视化图表: performance_comparison.png, aurocs_by_type.png, aurocs_vs_network_size.png\n")
  
} else {
  cat("\n⚠️ 没有成功评估的网络。\n")
}

# —— 10. 生成报告 ——
report_file <- file.path(out_root, "README.txt")
writeLines(c(
  "========================================",
  "网络性能评估报告 - GSE117221数据集",
  "========================================",
  paste("生成时间:", Sys.time()),
  "",
  "数据集信息:",
  "  - GSE117221 (H系列样本)",
  "  - 基因调控网络预测",
  "",
  "评估的网络类型:",
  "1. 单个方法网络 (3个): deepfgrn, 3dcema, deepsem",
  "2. 两两组合交集 (3个): deepfgrn_vs_3dcema, deepfgrn_vs_deepsem, 3dcema_vs_deepsem",
  "3. 三个方法交集 (1个): AllThree_Intersection",
  "4. 两个及以上方法 (1个): TwoOrMore_Union",
  "",
  "输入文件:",
  "  - deepfgrn: E:/SCD/其他数据/GSE117221/grn/deepfgrn_H_top10pct.tsv",
  "  - 3dcema: E:/SCD/其他数据/GSE117221/grn/3dcema_H_top10pct.tsv",
  "  - deepsem: E:/SCD/其他数据/GSE117221/grn/deepsem_H_top10pct.tsv",
  "",
  "金标准来源: Integrated_ChIP_Top200.txt",
  paste("金标准边数:", length(all_chip_pairs)),
  "",
  "输出文件:",
  "- ALL_NETWORKS_SUMMARY.tsv: 所有网络的汇总结果",
  "- performance_comparison.png: 性能比较柱状图",
  "- aurocs_by_type.png: 按网络类型的AUROC比较",
  "- aurocs_vs_network_size.png: AUROC与网络大小的关系",
  "- method_overlap_matrix.tsv: 方法间重叠矩阵",
  "",
  "各网络详细结果保存在对应的子目录中:",
  "  - ROC曲线和PR曲线",
  "  - 网络边列表",
  "  - 详细评估结果",
  "",
  paste("预测网络总数:", length(all_networks)),
  paste("成功评估网络数:", length(all_results))
), con = report_file)

cat("\n报告已保存至:", report_file, "\n")
























# -----------------------------------------------------------------------------
# 完整脚本：计算所有方法组合的 AUROC/AUPR (GSE117221 TI数据集)
# 包括：
#   1. 单个方法（3个）
#   2. 两两组合交集（3个）
#   3. 三个方法交集（1个）
#   4. 两个及以上方法（并集/平均，1个）
# 总计：8个网络
# -----------------------------------------------------------------------------

# —— 0. 清理环境 & 加载包 ——
rm(list = ls()); gc()

if (!requireNamespace("data.table", quietly=TRUE))   install.packages("data.table")
if (!requireNamespace("dplyr", quietly=TRUE))        install.packages("dplyr")
if (!requireNamespace("pROC", quietly=TRUE))         install.packages("pROC")
if (!requireNamespace("PRROC", quietly=TRUE))        install.packages("PRROC")
if (!requireNamespace("ggplot2", quietly=TRUE))      install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly=TRUE))        install.packages("tidyr")

library(data.table)
library(dplyr)
library(pROC)
library(PRROC)
library(ggplot2)
library(tidyr)

# —— 1. 定义清洗函数 ——
clean_name <- function(x) {
  x %>% trimws() %>% tolower() %>% gsub("[^a-z0-9]", "", .)
}

# —— 2. 输入：网络文件（3 个 top10pct.tsv）——
net_paths <- list(
  deepfgrn = "E:/SCD/其他数据/GSE117221/grn/deepfgrn_TI_top10pct.tsv",
  `3dcema` = "E:/SCD/其他数据/GSE117221/grn/3dcema_TI_top10pct.tsv",
  deepsem  = "E:/SCD/其他数据/GSE117221/grn/deepsem_TI_top10pct.tsv"
)

# 输出根目录
out_root <- "E:/SCD/其他数据/GSE117221/AUROC_TI"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# —— 3. 读取三个方法预测文件 & 构建 dt_list ——
dt_list <- lapply(names(net_paths), function(meth) {
  dt <- fread(net_paths[[meth]], sep = "\t", header = TRUE)
  
  # 兼容不同的权重列名
  if (!("EdgeWeight" %in% names(dt))) {
    cand <- intersect(names(dt), c("weight", "Weight", "score", "Score", "edgeweight", "Edge_Weight"))
    if (length(cand) >= 1) setnames(dt, cand[1], "EdgeWeight")
  }
  
  dt[, TF     := clean_name(TF)]
  dt[, Target := clean_name(Target)]
  dt[, pair   := paste0(TF, "_", Target)]
  
  # 若同一方法内同一 pair 多行，取均值聚合
  dt <- dt[, .(
    EdgeWeight = mean(EdgeWeight, na.rm = TRUE)
  ), by = .(TF, Target, pair)]
  
  # 添加方法名称列
  dt[, Method := meth]
  dt
})
names(dt_list) <- names(net_paths)

# —— 4. 读取金标准数据 ——
# 注意：GSE117221的金标准文件路径（请根据实际情况修改）
chip_file <- "E:/SCD/其他数据/gse133181_单细胞/GSE133181_RAW/chip/Integrated_ChIP_Top200.txt"

# 检查文件是否存在
if (!file.exists(chip_file)) {
  stop(paste("金标准文件不存在:", chip_file, 
             "\n请更新正确的金标准文件路径！"))
}

cat("读取金标准文件:", chip_file, "\n")
chip_dt <- fread(chip_file, sep = "\t", header = TRUE)

# 检查列名
if (!all(c("TF", "Target") %in% names(chip_dt))) {
  stop("ChIP文件必须包含 'TF' 和 'Target' 两列")
}

# 清洗并生成 pair
chip_dt[, TF_clean := clean_name(TF)]
chip_dt[, Target_clean := clean_name(Target)]
chip_dt[, pair := paste0(TF_clean, "_", Target_clean)]

# 获取所有金标准pair
all_chip_pairs <- unique(chip_dt$pair)

cat("金标准总边数:", nrow(chip_dt), "\n")
cat("去重后金标准边数:", length(all_chip_pairs), "\n")
cat("涉及TF数量:", length(unique(chip_dt$TF_clean)), "\n\n")

# —— 5. 定义评估函数 ——
evaluate_network <- function(network_dt, network_name, gold_pairs, out_dir) {
  cat("\n评估网络:", network_name, "\n")
  cat("  边数量:", nrow(network_dt), "\n")
  
  # 准备数据
  df <- copy(network_dt)[
    , `:=`(
      label = as.integer(pair %in% gold_pairs),
      score = EdgeWeight
    )
  ]
  
  # 统计
  pos_count <- sum(df$label)
  neg_count <- sum(1 - df$label)
  cat("  正样本数:", pos_count, "(", round(pos_count/nrow(df)*100, 2), "%)\n")
  cat("  负样本数:", neg_count, "(", round(neg_count/nrow(df)*100, 2), "%)\n")
  
  # 检查是否可以计算
  if (length(unique(df$label)) < 2) {
    cat("  ⚠️ 单一标签，无法计算 ROC/PR\n")
    return(NULL)
  }
  
  # 计算 AUROC
  roc_obj <- roc(df$label, df$score, quiet = TRUE, direction = "<")
  auroc <- auc(roc_obj)
  
  # 计算 AUPR
  pr_obj <- tryCatch({
    pr.curve(
      scores.class0 = df$score[df$label == 1],
      scores.class1 = df$score[df$label == 0],
      curve = TRUE
    )
  }, error = function(e) {
    cat("  警告: PR计算出错，尝试替代方法\n")
    return(NULL)
  })
  
  if (!is.null(pr_obj) && is.list(pr_obj) && !is.null(pr_obj$auc.integral)) {
    aupr <- pr_obj$auc.integral
  } else {
    # 使用pROC计算近似AUPR
    pr_roc <- roc(df$label, df$score, quiet = TRUE)
    coords_df <- coords(pr_roc, "all", ret = c("precision", "recall"))
    coords_df <- na.omit(coords_df)
    if (nrow(coords_df) > 1) {
      aupr <- sum(diff(coords_df$recall) * 
                    (head(coords_df$precision, -1) + tail(coords_df$precision, -1)) / 2)
    } else {
      aupr <- NA
    }
  }
  
  # 创建网络专属目录
  net_dir <- file.path(out_dir, gsub("[^a-zA-Z0-9]", "_", network_name))
  dir.create(net_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 绘制 ROC 曲线
  png(file.path(net_dir, "ROC_curve.png"), 800, 600)
  plot(roc_obj, 
       main = sprintf("%s\nAUROC = %.4f", network_name, auroc),
       col = "blue", lwd = 2, legacy.axes = TRUE)
  abline(a = 0, b = 1, lty = 2, col = "gray")
  legend("bottomright", legend = c("ROC曲线", "随机分类器"), 
         col = c("blue", "gray"), lty = c(1, 2), lwd = c(2, 1))
  dev.off()
  
  # 绘制 PR 曲线（如果成功）
  if (!is.null(pr_obj) && is.list(pr_obj) && !is.null(pr_obj$curve)) {
    png(file.path(net_dir, "PR_curve.png"), 800, 600)
    plot(pr_obj$curve[, 1], pr_obj$curve[, 2],
         type = "l", lwd = 2, col = "red",
         xlab = "Recall (灵敏度)", ylab = "Precision (精确率)",
         main = sprintf("%s\nAUPR = %.4f", network_name, aupr))
    grid()
    dev.off()
  }
  
  # 保存结果
  result <- data.table(
    Network = network_name,
    Type = ifelse(grepl("_vs_", network_name), "Intersection",
                  ifelse(grepl("TwoOrMore", network_name), "Union_2plus",
                         ifelse(grepl("AllThree", network_name), "Intersection_3", "Single"))),
    N_Edges = nrow(df),
    Pos_Edges = pos_count,
    Neg_Edges = neg_count,
    Pos_Ratio = round(pos_count/nrow(df)*100, 4),
    AUROC = round(as.numeric(auroc), 4),
    AUPR = round(aupr, 4)
  )
  
  # 保存网络文件
  fwrite(network_dt, file.path(net_dir, "network_edges.tsv"), sep = "\t", quote = FALSE)
  
  # 保存结果
  fwrite(result, file.path(net_dir, "evaluation_results.tsv"), sep = "\t", quote = FALSE)
  
  cat("  ✓ AUROC =", round(auroc, 4), ", AUPR =", round(aupr, 4), "\n")
  
  return(result)
}

# —— 6. 构建所有网络 ——

cat("\n", rep("=", 80), "\n", sep = "")
cat("开始构建所有网络组合\n")
cat(rep("=", 80), "\n", sep = "")

all_networks <- list()
method_names <- names(net_paths)

# 6.1 单个方法网络
cat("\n1. 构建单个方法网络...\n")
for (meth in method_names) {
  network_dt <- dt_list[[meth]][, .(TF, Target, pair, EdgeWeight)]
  all_networks[[meth]] <- network_dt
  cat("  ✓", meth, ":", nrow(network_dt), "条边\n")
}

# 6.2 两两组合交集
cat("\n2. 构建两两组合交集网络...\n")
combinations <- combn(method_names, 2, simplify = FALSE)
for (comb in combinations) {
  comb_name <- paste(comb, collapse = "_vs_")
  
  # 获取两个方法的pair集合
  pairs1 <- dt_list[[comb[1]]]$pair
  pairs2 <- dt_list[[comb[2]]]$pair
  
  # 计算交集
  intersect_pairs <- intersect(pairs1, pairs2)
  
  if (length(intersect_pairs) > 0) {
    # 提取交集边的信息，权重取两个方法的平均值
    dt1_sub <- dt_list[[comb[1]]][pair %in% intersect_pairs, .(pair, TF, Target, Weight1 = EdgeWeight)]
    dt2_sub <- dt_list[[comb[2]]][pair %in% intersect_pairs, .(pair, Weight2 = EdgeWeight)]
    
    network_dt <- merge(dt1_sub, dt2_sub, by = "pair")
    network_dt[, EdgeWeight := (Weight1 + Weight2) / 2]
    network_dt <- network_dt[, .(TF, Target, pair, EdgeWeight)]
    
    all_networks[[comb_name]] <- network_dt
    cat("  ✓", comb_name, ":", nrow(network_dt), "条边\n")
  } else {
    cat("  ✗", comb_name, ": 交集为空\n")
  }
}

# 6.3 三个方法交集
cat("\n3. 构建三个方法交集网络...\n")
pairs_all <- lapply(dt_list, function(x) x$pair)
common_pairs <- Reduce(intersect, pairs_all)

if (length(common_pairs) > 0) {
  # 提取三个方法的权重并取平均
  network_dt <- data.table(pair = common_pairs)
  for (meth in method_names) {
    dt_sub <- dt_list[[meth]][pair %in% common_pairs, .(pair, Weight = EdgeWeight)]
    setnames(dt_sub, "Weight", paste0("Weight_", meth))
    network_dt <- merge(network_dt, dt_sub, by = "pair", all.x = TRUE)
  }
  
  # 计算平均权重
  weight_cols <- paste0("Weight_", method_names)
  network_dt[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = weight_cols]
  
  # 添加TF和Target信息
  network_dt <- merge(network_dt, dt_list[[1]][, .(pair, TF, Target)], by = "pair")
  network_dt <- network_dt[, .(TF, Target, pair, EdgeWeight)]
  
  all_networks[["AllThree_Intersection"]] <- network_dt
  cat("  ✓ AllThree_Intersection:", nrow(network_dt), "条边\n")
} else {
  cat("  ✗ AllThree_Intersection: 交集为空\n")
}

# 6.4 两个及以上方法（保留至少两个方法预测的边，权重取平均）
cat("\n4. 构建两个及以上方法网络...\n")

# 合并所有方法的权重
merged_dt <- data.table(pair = unique(unlist(lapply(dt_list, function(x) x$pair))))

for (meth in method_names) {
  dt_sub <- dt_list[[meth]][, .(pair, Weight = EdgeWeight)]
  setnames(dt_sub, "Weight", paste0("Weight_", meth))
  merged_dt <- merge(merged_dt, dt_sub, by = "pair", all.x = TRUE)
}

# 计算每个pair出现在几个方法中
weight_cols <- paste0("Weight_", method_names)
merged_dt[, count_present := rowSums(!is.na(.SD)), .SDcols = weight_cols]

# 筛选出现次数 >= 2 的边
two_or_more <- merged_dt[count_present >= 2]

if (nrow(two_or_more) > 0) {
  # 计算平均权重
  two_or_more[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = weight_cols]
  
  # 添加TF和Target信息
  two_or_more <- merge(two_or_more, dt_list[[1]][, .(pair, TF, Target)], by = "pair")
  network_dt <- two_or_more[, .(TF, Target, pair, EdgeWeight)]
  
  all_networks[["TwoOrMore_Union"]] <- network_dt
  cat("  ✓ TwoOrMore_Union:", nrow(network_dt), "条边\n")
} else {
  cat("  ✗ TwoOrMore_Union: 没有边\n")
}

# —— 7. 对所有网络进行评估 ——
cat("\n", rep("=", 80), "\n", sep = "")
cat("开始评估所有网络\n")
cat(rep("=", 80), "\n", sep = "")

all_results <- list()

for (net_name in names(all_networks)) {
  result <- evaluate_network(all_networks[[net_name]], net_name, all_chip_pairs, out_root)
  if (!is.null(result)) {
    all_results[[net_name]] <- result
  }
}

# —— 8. 汇总所有结果 ——
cat("\n", rep("=", 80), "\n", sep = "")
cat("汇总结果\n")
cat(rep("=", 80), "\n", sep = "")

if (length(all_results) > 0) {
  summary_all <- rbindlist(all_results)
  
  # 按AUROC排序
  setorder(summary_all, -AUROC)
  
  # 保存汇总结果
  fwrite(summary_all, file.path(out_root, "ALL_NETWORKS_SUMMARY.tsv"), 
         sep = "\t", quote = FALSE)
  
  # 打印汇总表格
  cat("\n网络性能汇总表:\n")
  print(summary_all)
  
  # —— 9. 可视化 ——
  
  # 9.1 AUROC和AUPR柱状图
  summary_plot <- summary_all %>%
    pivot_longer(cols = c(AUROC, AUPR), names_to = "Metric", values_to = "Value")
  
  p1 <- ggplot(summary_plot, aes(x = reorder(Network, Value), y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "GSE117221 TI数据集: 所有网络性能比较",
         x = "网络", y = "分数") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "performance_comparison.png"), p1, width = 12, height = 7)
  
  # 9.2 按类型分组比较
  p2 <- ggplot(summary_all, aes(x = Type, y = AUROC, fill = Type)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    theme_minimal() +
    labs(title = "GSE117221 TI数据集: 不同网络类型的AUROC比较",
         x = "网络类型", y = "AUROC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "aurocs_by_type.png"), p2, width = 8, height = 6)
  
  # 9.3 AUROC vs 网络大小的关系
  p3 <- ggplot(summary_all, aes(x = N_Edges, y = AUROC, color = Type, label = Network)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
    theme_minimal() +
    labs(title = "GSE117221 TI数据集: AUROC vs 网络大小",
         x = "网络边数", y = "AUROC") +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "aurocs_vs_network_size.png"), p3, width = 10, height = 7)
  
  # 9.4 热图展示各方法的重叠情况（可选）
  # 创建重叠矩阵
  overlap_matrix <- matrix(0, nrow = length(method_names), ncol = length(method_names))
  rownames(overlap_matrix) <- method_names
  colnames(overlap_matrix) <- method_names
  
  for (i in 1:length(method_names)) {
    for (j in 1:length(method_names)) {
      if (i <= j) {
        pairs_i <- dt_list[[method_names[i]]]$pair
        pairs_j <- dt_list[[method_names[j]]]$pair
        overlap <- length(intersect(pairs_i, pairs_j))
        total <- length(union(pairs_i, pairs_j))
        overlap_matrix[i, j] <- overlap / total * 100
        overlap_matrix[j, i] <- overlap / total * 100
      }
    }
  }
  
  # 保存重叠矩阵
  fwrite(as.data.table(overlap_matrix, keep.rownames = "Method"), 
         file.path(out_root, "method_overlap_matrix.tsv"), sep = "\t")
  
  # 9.5 打印最佳网络信息
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("最佳网络 (按AUROC):\n")
  best_network <- summary_all[which.max(AUROC)]
  print(best_network)
  
  cat("\n最佳网络 (按AUPR):\n")
  best_pr_network <- summary_all[which.max(AUPR)]
  print(best_pr_network)
  
  cat("\n网络大小与性能:\n")
  cat("  最小网络:", summary_all[which.min(N_Edges), .(Network, N_Edges, AUROC)] %>% 
        mutate(across(where(is.numeric), ~round(., 4))), "\n")
  cat("  最大网络:", summary_all[which.max(N_Edges), .(Network, N_Edges, AUROC)] %>% 
        mutate(across(where(is.numeric), ~round(., 4))), "\n")
  
  cat("\n✔️ 所有评估完成！\n")
  cat("结果保存在:", out_root, "\n")
  cat("   - 汇总文件: ALL_NETWORKS_SUMMARY.tsv\n")
  cat("   - 各网络详细结果在对应子目录中\n")
  cat("   - 可视化图表: performance_comparison.png, aurocs_by_type.png, aurocs_vs_network_size.png\n")
  
} else {
  cat("\n⚠️ 没有成功评估的网络。\n")
}

# —— 10. 生成报告 ——
report_file <- file.path(out_root, "README.txt")
writeLines(c(
  "========================================",
  "网络性能评估报告 - GSE117221 TI数据集",
  "========================================",
  paste("生成时间:", Sys.time()),
  "",
  "数据集信息:",
  "  - GSE117221 (TI系列样本)",
  "  - 基因调控网络预测",
  "",
  "评估的网络类型:",
  "1. 单个方法网络 (3个): deepfgrn, 3dcema, deepsem",
  "2. 两两组合交集 (3个): deepfgrn_vs_3dcema, deepfgrn_vs_deepsem, 3dcema_vs_deepsem",
  "3. 三个方法交集 (1个): AllThree_Intersection",
  "4. 两个及以上方法 (1个): TwoOrMore_Union",
  "",
  "输入文件:",
  "  - deepfgrn: E:/SCD/其他数据/GSE117221/grn/deepfgrn_TI_top10pct.tsv",
  "  - 3dcema: E:/SCD/其他数据/GSE117221/grn/3dcema_TI_top10pct.tsv",
  "  - deepsem: E:/SCD/其他数据/GSE117221/grn/deepsem_TI_top10pct.tsv",
  "",
  "金标准来源: Integrated_ChIP_Top200.txt",
  paste("金标准边数:", length(all_chip_pairs)),
  "",
  "输出文件:",
  "- ALL_NETWORKS_SUMMARY.tsv: 所有网络的汇总结果",
  "- performance_comparison.png: 性能比较柱状图",
  "- aurocs_by_type.png: 按网络类型的AUROC比较",
  "- aurocs_vs_network_size.png: AUROC与网络大小的关系",
  "- method_overlap_matrix.tsv: 方法间重叠矩阵",
  "",
  "各网络详细结果保存在对应的子目录中:",
  "  - ROC曲线和PR曲线",
  "  - 网络边列表",
  "  - 详细评估结果",
  "",
  paste("预测网络总数:", length(all_networks)),
  paste("成功评估网络数:", length(all_results))
), con = report_file)

cat("\n报告已保存至:", report_file, "\n")


























# -----------------------------------------------------------------------------
# 完整脚本：计算所有方法组合的 AUROC/AUPR (GSE117221 TM数据集)
# 包括：
#   1. 单个方法（3个）
#   2. 两两组合交集（3个）
#   3. 三个方法交集（1个）
#   4. 两个及以上方法（并集/平均，1个）
# 总计：8个网络
# -----------------------------------------------------------------------------

# —— 0. 清理环境 & 加载包 ——
rm(list = ls()); gc()

if (!requireNamespace("data.table", quietly=TRUE))   install.packages("data.table")
if (!requireNamespace("dplyr", quietly=TRUE))        install.packages("dplyr")
if (!requireNamespace("pROC", quietly=TRUE))         install.packages("pROC")
if (!requireNamespace("PRROC", quietly=TRUE))        install.packages("PRROC")
if (!requireNamespace("ggplot2", quietly=TRUE))      install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly=TRUE))        install.packages("tidyr")

library(data.table)
library(dplyr)
library(pROC)
library(PRROC)
library(ggplot2)
library(tidyr)

# —— 1. 定义清洗函数 ——
clean_name <- function(x) {
  x %>% trimws() %>% tolower() %>% gsub("[^a-z0-9]", "", .)
}

# —— 2. 输入：网络文件（3 个 top10pct.tsv）——
net_paths <- list(
  deepfgrn = "E:/SCD/其他数据/GSE117221/grn/deepfgrn_TM_top10pct.tsv",
  `3dcema` = "E:/SCD/其他数据/GSE117221/grn/3dcema_TM_top10pct.tsv",
  deepsem  = "E:/SCD/其他数据/GSE117221/grn/deepsem_TM_top10pct.tsv"
)

# 输出根目录
out_root <- "E:/SCD/其他数据/GSE117221/AUROC_TM"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# —— 3. 读取三个方法预测文件 & 构建 dt_list ——
dt_list <- lapply(names(net_paths), function(meth) {
  dt <- fread(net_paths[[meth]], sep = "\t", header = TRUE)
  
  # 兼容不同的权重列名
  if (!("EdgeWeight" %in% names(dt))) {
    cand <- intersect(names(dt), c("weight", "Weight", "score", "Score", "edgeweight", "Edge_Weight"))
    if (length(cand) >= 1) setnames(dt, cand[1], "EdgeWeight")
  }
  
  dt[, TF     := clean_name(TF)]
  dt[, Target := clean_name(Target)]
  dt[, pair   := paste0(TF, "_", Target)]
  
  # 若同一方法内同一 pair 多行，取均值聚合
  dt <- dt[, .(
    EdgeWeight = mean(EdgeWeight, na.rm = TRUE)
  ), by = .(TF, Target, pair)]
  
  # 添加方法名称列
  dt[, Method := meth]
  dt
})
names(dt_list) <- names(net_paths)

# —— 4. 读取金标准数据 ——
# 注意：GSE117221的金标准文件路径（请根据实际情况修改）
chip_file <- "E:/SCD/其他数据/gse133181_单细胞/GSE133181_RAW/chip/Integrated_ChIP_Top200.txt"

# 检查文件是否存在
if (!file.exists(chip_file)) {
  stop(paste("金标准文件不存在:", chip_file, 
             "\n请更新正确的金标准文件路径！"))
}

cat("读取金标准文件:", chip_file, "\n")
chip_dt <- fread(chip_file, sep = "\t", header = TRUE)

# 检查列名
if (!all(c("TF", "Target") %in% names(chip_dt))) {
  stop("ChIP文件必须包含 'TF' 和 'Target' 两列")
}

# 清洗并生成 pair
chip_dt[, TF_clean := clean_name(TF)]
chip_dt[, Target_clean := clean_name(Target)]
chip_dt[, pair := paste0(TF_clean, "_", Target_clean)]

# 获取所有金标准pair
all_chip_pairs <- unique(chip_dt$pair)

cat("金标准总边数:", nrow(chip_dt), "\n")
cat("去重后金标准边数:", length(all_chip_pairs), "\n")
cat("涉及TF数量:", length(unique(chip_dt$TF_clean)), "\n\n")

# —— 5. 定义评估函数 ——
evaluate_network <- function(network_dt, network_name, gold_pairs, out_dir) {
  cat("\n评估网络:", network_name, "\n")
  cat("  边数量:", nrow(network_dt), "\n")
  
  # 准备数据
  df <- copy(network_dt)[
    , `:=`(
      label = as.integer(pair %in% gold_pairs),
      score = EdgeWeight
    )
  ]
  
  # 统计
  pos_count <- sum(df$label)
  neg_count <- sum(1 - df$label)
  cat("  正样本数:", pos_count, "(", round(pos_count/nrow(df)*100, 2), "%)\n")
  cat("  负样本数:", neg_count, "(", round(neg_count/nrow(df)*100, 2), "%)\n")
  
  # 检查是否可以计算
  if (length(unique(df$label)) < 2) {
    cat("  ⚠️ 单一标签，无法计算 ROC/PR\n")
    return(NULL)
  }
  
  # 计算 AUROC
  roc_obj <- roc(df$label, df$score, quiet = TRUE, direction = "<")
  auroc <- auc(roc_obj)
  
  # 计算 AUPR
  pr_obj <- tryCatch({
    pr.curve(
      scores.class0 = df$score[df$label == 1],
      scores.class1 = df$score[df$label == 0],
      curve = TRUE
    )
  }, error = function(e) {
    cat("  警告: PR计算出错，尝试替代方法\n")
    return(NULL)
  })
  
  if (!is.null(pr_obj) && is.list(pr_obj) && !is.null(pr_obj$auc.integral)) {
    aupr <- pr_obj$auc.integral
  } else {
    # 使用pROC计算近似AUPR
    pr_roc <- roc(df$label, df$score, quiet = TRUE)
    coords_df <- coords(pr_roc, "all", ret = c("precision", "recall"))
    coords_df <- na.omit(coords_df)
    if (nrow(coords_df) > 1) {
      aupr <- sum(diff(coords_df$recall) * 
                    (head(coords_df$precision, -1) + tail(coords_df$precision, -1)) / 2)
    } else {
      aupr <- NA
    }
  }
  
  # 创建网络专属目录
  net_dir <- file.path(out_dir, gsub("[^a-zA-Z0-9]", "_", network_name))
  dir.create(net_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 绘制 ROC 曲线
  png(file.path(net_dir, "ROC_curve.png"), 800, 600)
  plot(roc_obj, 
       main = sprintf("%s\nAUROC = %.4f", network_name, auroc),
       col = "blue", lwd = 2, legacy.axes = TRUE)
  abline(a = 0, b = 1, lty = 2, col = "gray")
  legend("bottomright", legend = c("ROC曲线", "随机分类器"), 
         col = c("blue", "gray"), lty = c(1, 2), lwd = c(2, 1))
  dev.off()
  
  # 绘制 PR 曲线（如果成功）
  if (!is.null(pr_obj) && is.list(pr_obj) && !is.null(pr_obj$curve)) {
    png(file.path(net_dir, "PR_curve.png"), 800, 600)
    plot(pr_obj$curve[, 1], pr_obj$curve[, 2],
         type = "l", lwd = 2, col = "red",
         xlab = "Recall (灵敏度)", ylab = "Precision (精确率)",
         main = sprintf("%s\nAUPR = %.4f", network_name, aupr))
    grid()
    dev.off()
  }
  
  # 保存结果
  result <- data.table(
    Network = network_name,
    Type = ifelse(grepl("_vs_", network_name), "Intersection",
                  ifelse(grepl("TwoOrMore", network_name), "Union_2plus",
                         ifelse(grepl("AllThree", network_name), "Intersection_3", "Single"))),
    N_Edges = nrow(df),
    Pos_Edges = pos_count,
    Neg_Edges = neg_count,
    Pos_Ratio = round(pos_count/nrow(df)*100, 4),
    AUROC = round(as.numeric(auroc), 4),
    AUPR = round(aupr, 4)
  )
  
  # 保存网络文件
  fwrite(network_dt, file.path(net_dir, "network_edges.tsv"), sep = "\t", quote = FALSE)
  
  # 保存结果
  fwrite(result, file.path(net_dir, "evaluation_results.tsv"), sep = "\t", quote = FALSE)
  
  cat("  ✓ AUROC =", round(auroc, 4), ", AUPR =", round(aupr, 4), "\n")
  
  return(result)
}

# —— 6. 构建所有网络 ——

cat("\n", rep("=", 80), "\n", sep = "")
cat("开始构建所有网络组合\n")
cat(rep("=", 80), "\n", sep = "")

all_networks <- list()
method_names <- names(net_paths)

# 6.1 单个方法网络
cat("\n1. 构建单个方法网络...\n")
for (meth in method_names) {
  network_dt <- dt_list[[meth]][, .(TF, Target, pair, EdgeWeight)]
  all_networks[[meth]] <- network_dt
  cat("  ✓", meth, ":", nrow(network_dt), "条边\n")
}

# 6.2 两两组合交集
cat("\n2. 构建两两组合交集网络...\n")
combinations <- combn(method_names, 2, simplify = FALSE)
for (comb in combinations) {
  comb_name <- paste(comb, collapse = "_vs_")
  
  # 获取两个方法的pair集合
  pairs1 <- dt_list[[comb[1]]]$pair
  pairs2 <- dt_list[[comb[2]]]$pair
  
  # 计算交集
  intersect_pairs <- intersect(pairs1, pairs2)
  
  if (length(intersect_pairs) > 0) {
    # 提取交集边的信息，权重取两个方法的平均值
    dt1_sub <- dt_list[[comb[1]]][pair %in% intersect_pairs, .(pair, TF, Target, Weight1 = EdgeWeight)]
    dt2_sub <- dt_list[[comb[2]]][pair %in% intersect_pairs, .(pair, Weight2 = EdgeWeight)]
    
    network_dt <- merge(dt1_sub, dt2_sub, by = "pair")
    network_dt[, EdgeWeight := (Weight1 + Weight2) / 2]
    network_dt <- network_dt[, .(TF, Target, pair, EdgeWeight)]
    
    all_networks[[comb_name]] <- network_dt
    cat("  ✓", comb_name, ":", nrow(network_dt), "条边\n")
  } else {
    cat("  ✗", comb_name, ": 交集为空\n")
  }
}

# 6.3 三个方法交集
cat("\n3. 构建三个方法交集网络...\n")
pairs_all <- lapply(dt_list, function(x) x$pair)
common_pairs <- Reduce(intersect, pairs_all)

if (length(common_pairs) > 0) {
  # 提取三个方法的权重并取平均
  network_dt <- data.table(pair = common_pairs)
  for (meth in method_names) {
    dt_sub <- dt_list[[meth]][pair %in% common_pairs, .(pair, Weight = EdgeWeight)]
    setnames(dt_sub, "Weight", paste0("Weight_", meth))
    network_dt <- merge(network_dt, dt_sub, by = "pair", all.x = TRUE)
  }
  
  # 计算平均权重
  weight_cols <- paste0("Weight_", method_names)
  network_dt[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = weight_cols]
  
  # 添加TF和Target信息
  network_dt <- merge(network_dt, dt_list[[1]][, .(pair, TF, Target)], by = "pair")
  network_dt <- network_dt[, .(TF, Target, pair, EdgeWeight)]
  
  all_networks[["AllThree_Intersection"]] <- network_dt
  cat("  ✓ AllThree_Intersection:", nrow(network_dt), "条边\n")
} else {
  cat("  ✗ AllThree_Intersection: 交集为空\n")
}

# 6.4 两个及以上方法（保留至少两个方法预测的边，权重取平均）
cat("\n4. 构建两个及以上方法网络...\n")

# 合并所有方法的权重
merged_dt <- data.table(pair = unique(unlist(lapply(dt_list, function(x) x$pair))))

for (meth in method_names) {
  dt_sub <- dt_list[[meth]][, .(pair, Weight = EdgeWeight)]
  setnames(dt_sub, "Weight", paste0("Weight_", meth))
  merged_dt <- merge(merged_dt, dt_sub, by = "pair", all.x = TRUE)
}

# 计算每个pair出现在几个方法中
weight_cols <- paste0("Weight_", method_names)
merged_dt[, count_present := rowSums(!is.na(.SD)), .SDcols = weight_cols]

# 筛选出现次数 >= 2 的边
two_or_more <- merged_dt[count_present >= 2]

if (nrow(two_or_more) > 0) {
  # 计算平均权重
  two_or_more[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = weight_cols]
  
  # 添加TF和Target信息
  two_or_more <- merge(two_or_more, dt_list[[1]][, .(pair, TF, Target)], by = "pair")
  network_dt <- two_or_more[, .(TF, Target, pair, EdgeWeight)]
  
  all_networks[["TwoOrMore_Union"]] <- network_dt
  cat("  ✓ TwoOrMore_Union:", nrow(network_dt), "条边\n")
} else {
  cat("  ✗ TwoOrMore_Union: 没有边\n")
}

# —— 7. 对所有网络进行评估 ——
cat("\n", rep("=", 80), "\n", sep = "")
cat("开始评估所有网络\n")
cat(rep("=", 80), "\n", sep = "")

all_results <- list()

for (net_name in names(all_networks)) {
  result <- evaluate_network(all_networks[[net_name]], net_name, all_chip_pairs, out_root)
  if (!is.null(result)) {
    all_results[[net_name]] <- result
  }
}

# —— 8. 汇总所有结果 ——
cat("\n", rep("=", 80), "\n", sep = "")
cat("汇总结果\n")
cat(rep("=", 80), "\n", sep = "")

if (length(all_results) > 0) {
  summary_all <- rbindlist(all_results)
  
  # 按AUROC排序
  setorder(summary_all, -AUROC)
  
  # 保存汇总结果
  fwrite(summary_all, file.path(out_root, "ALL_NETWORKS_SUMMARY.tsv"), 
         sep = "\t", quote = FALSE)
  
  # 打印汇总表格
  cat("\n网络性能汇总表:\n")
  print(summary_all)
  
  # —— 9. 可视化 ——
  
  # 9.1 AUROC和AUPR柱状图
  summary_plot <- summary_all %>%
    pivot_longer(cols = c(AUROC, AUPR), names_to = "Metric", values_to = "Value")
  
  p1 <- ggplot(summary_plot, aes(x = reorder(Network, Value), y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "GSE117221 TM数据集: 所有网络性能比较",
         x = "网络", y = "分数") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "performance_comparison.png"), p1, width = 12, height = 7)
  
  # 9.2 按类型分组比较
  p2 <- ggplot(summary_all, aes(x = Type, y = AUROC, fill = Type)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    theme_minimal() +
    labs(title = "GSE117221 TM数据集: 不同网络类型的AUROC比较",
         x = "网络类型", y = "AUROC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "aurocs_by_type.png"), p2, width = 8, height = 6)
  
  # 9.3 AUROC vs 网络大小的关系
  p3 <- ggplot(summary_all, aes(x = N_Edges, y = AUROC, color = Type, label = Network)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
    theme_minimal() +
    labs(title = "GSE117221 TM数据集: AUROC vs 网络大小",
         x = "网络边数", y = "AUROC") +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "aurocs_vs_network_size.png"), p3, width = 10, height = 7)
  
  # 9.4 热图展示各方法的重叠情况（可选）
  # 创建重叠矩阵
  overlap_matrix <- matrix(0, nrow = length(method_names), ncol = length(method_names))
  rownames(overlap_matrix) <- method_names
  colnames(overlap_matrix) <- method_names
  
  for (i in 1:length(method_names)) {
    for (j in 1:length(method_names)) {
      if (i <= j) {
        pairs_i <- dt_list[[method_names[i]]]$pair
        pairs_j <- dt_list[[method_names[j]]]$pair
        overlap <- length(intersect(pairs_i, pairs_j))
        total <- length(union(pairs_i, pairs_j))
        overlap_matrix[i, j] <- overlap / total * 100
        overlap_matrix[j, i] <- overlap / total * 100
      }
    }
  }
  
  # 保存重叠矩阵
  fwrite(as.data.table(overlap_matrix, keep.rownames = "Method"), 
         file.path(out_root, "method_overlap_matrix.tsv"), sep = "\t")
  
  # 9.5 打印最佳网络信息
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("最佳网络 (按AUROC):\n")
  best_network <- summary_all[which.max(AUROC)]
  print(best_network)
  
  cat("\n最佳网络 (按AUPR):\n")
  best_pr_network <- summary_all[which.max(AUPR)]
  print(best_pr_network)
  
  cat("\n网络大小与性能:\n")
  cat("  最小网络:", summary_all[which.min(N_Edges), .(Network, N_Edges, AUROC)] %>% 
        mutate(across(where(is.numeric), ~round(., 4))), "\n")
  cat("  最大网络:", summary_all[which.max(N_Edges), .(Network, N_Edges, AUROC)] %>% 
        mutate(across(where(is.numeric), ~round(., 4))), "\n")
  
  cat("\n✔️ 所有评估完成！\n")
  cat("结果保存在:", out_root, "\n")
  cat("   - 汇总文件: ALL_NETWORKS_SUMMARY.tsv\n")
  cat("   - 各网络详细结果在对应子目录中\n")
  cat("   - 可视化图表: performance_comparison.png, aurocs_by_type.png, aurocs_vs_network_size.png\n")
  
} else {
  cat("\n⚠️ 没有成功评估的网络。\n")
}

# —— 10. 生成报告 ——
report_file <- file.path(out_root, "README.txt")
writeLines(c(
  "========================================",
  "网络性能评估报告 - GSE117221 TM数据集",
  "========================================",
  paste("生成时间:", Sys.time()),
  "",
  "数据集信息:",
  "  - GSE117221 (TM系列样本)",
  "  - 基因调控网络预测",
  "",
  "评估的网络类型:",
  "1. 单个方法网络 (3个): deepfgrn, 3dcema, deepsem",
  "2. 两两组合交集 (3个): deepfgrn_vs_3dcema, deepfgrn_vs_deepsem, 3dcema_vs_deepsem",
  "3. 三个方法交集 (1个): AllThree_Intersection",
  "4. 两个及以上方法 (1个): TwoOrMore_Union",
  "",
  "输入文件:",
  "  - deepfgrn: E:/SCD/其他数据/GSE117221/grn/deepfgrn_TM_top10pct.tsv",
  "  - 3dcema: E:/SCD/其他数据/GSE117221/grn/3dcema_TM_top10pct.tsv",
  "  - deepsem: E:/SCD/其他数据/GSE117221/grn/deepsem_TM_top10pct.tsv",
  "",
  "金标准来源: Integrated_ChIP_Top200.txt",
  paste("金标准边数:", length(all_chip_pairs)),
  "",
  "输出文件:",
  "- ALL_NETWORKS_SUMMARY.tsv: 所有网络的汇总结果",
  "- performance_comparison.png: 性能比较柱状图",
  "- aurocs_by_type.png: 按网络类型的AUROC比较",
  "- aurocs_vs_network_size.png: AUROC与网络大小的关系",
  "- method_overlap_matrix.tsv: 方法间重叠矩阵",
  "",
  "各网络详细结果保存在对应的子目录中:",
  "  - ROC曲线和PR曲线",
  "  - 网络边列表",
  "  - 详细评估结果",
  "",
  paste("预测网络总数:", length(all_networks)),
  paste("成功评估网络数:", length(all_results))
), con = report_file)

cat("\n报告已保存至:", report_file, "\n")