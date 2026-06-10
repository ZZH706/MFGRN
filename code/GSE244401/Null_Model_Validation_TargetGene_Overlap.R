# -----------------------------------------------------------------------------
# 提取 4 个转录因子的 ChIP-seq 对应靶基因
# 输入：BED 文件
# 输出：每个 TF 一个靶基因列表文件，保存在 E:/SCD/数据/CHIP-seq数据
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ChIPseeker)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

# -----------------------------
# 1) 参数
# -----------------------------
tf_list <- list(
  list(name = "BCL11A", bed = "E:/SCD/数据/CHIP-seq数据/GSE103445_HUDEP2_BCL11A_optimal_idr_peaks.bed"),
  list(name = "GATA1",  bed = "E:/SCD/数据/CHIP-seq数据/GSM970257_GATA1-F_peaks.bed"),
  list(name = "TAL1",   bed = "E:/SCD/数据/CHIP-seq数据/GSM1816083_TAL1-F5.peak.bed"),
  list(name = "NFE2",   bed = "E:/SCD/数据/CHIP-seq数据/GSM1816086_NFE2-F5.peak.bed")
)

out_dir <- "E:/SCD/数据/CHIP-seq数据"

# 规范化路径
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
tf_list <- lapply(tf_list, function(x) {
  x$bed <- normalizePath(x$bed, winslash = "/", mustWork = FALSE)
  x
})

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2) 名称清洗函数
# -----------------------------
clean_name <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(x)
  x <- gsub("[^A-Z0-9]", "", x)
  x
}

# -----------------------------
# 3) 从 BED 注释到靶基因
# -----------------------------
get_chip_targets <- function(bed_path, tss_region = c(-10000, 10000)) {
  bed_raw <- fread(bed_path, header = FALSE)[, 1:3]
  colnames(bed_raw) <- c("chr", "start", "end")
  
  bed_clean <- bed_raw %>%
    filter(!is.na(chr), !is.na(start), !is.na(end)) %>%
    filter(suppressWarnings(!is.na(as.numeric(start)))) %>%
    filter(suppressWarnings(!is.na(as.numeric(end)))) %>%
    mutate(
      start = as.numeric(start),
      end   = as.numeric(end)
    ) %>%
    filter(start <= end)
  
  if (nrow(bed_clean) == 0) {
    stop("清洗后 BED 文件为空：", bed_path)
  }
  
  tmp_bed <- tempfile(fileext = ".bed")
  fwrite(bed_clean, tmp_bed, sep = "\t", col.names = FALSE, quote = FALSE)
  
  peakAnno <- annotatePeak(
    tmp_bed,
    tssRegion = tss_region,
    TxDb      = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb    = "org.Hs.eg.db"
  )
  
  anno_df <- as.data.frame(peakAnno)
  genes <- anno_df$SYMBOL
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- unique(clean_name(genes))
  
  return(genes)
}

# -----------------------------
# 4) 主流程：提取每个 TF 的靶基因并保存
# -----------------------------
summary_list <- list()

for (tf in tf_list) {
  cat("\n开始处理：", tf$name, "\n")
  
  targets <- get_chip_targets(tf$bed, tss_region = c(-2000, 2000))
  
  target_df <- data.frame(
    TF = tf$name,
    Target = targets,
    stringsAsFactors = FALSE
  )
  
  out_file <- file.path(out_dir, paste0(tf$name, "_ChIP_targets.tsv"))
  fwrite(target_df, out_file, sep = "\t", quote = FALSE)
  
  cat("靶基因数：", nrow(target_df), "\n")
  cat("结果已保存：", out_file, "\n")
  
  summary_list[[tf$name]] <- data.frame(
    TF = tf$name,
    Target_Count = nrow(target_df),
    Output_File = out_file,
    stringsAsFactors = FALSE
  )
}

# 保存汇总文件
summary_df <- do.call(rbind, summary_list)
fwrite(summary_df,
       file = file.path(out_dir, "TF_ChIP_targets_summary.tsv"),
       sep = "\t", quote = FALSE)

cat("\n全部完成。\n")\





























# -----------------------------------------------------------------------------
# 为每个表达文件分别构造 TF-Target 调控对
# 输入文件：第一列为基因名，后面为样本表达值
# 输出文件：每个文件各自一个 TF-Target 对文件
# 同时输出：每个文件中每个 TF 的调控对数量统计
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# 1. 输入文件
expr_files <- c(
  "E:/SCD/数据/构建网络的数据/data_c1_no0_transposed.csv",
  "E:/SCD/数据/构建网络的数据/data_c2_no0_transposed.csv",
  "E:/SCD/数据/构建网络的数据/data_s1_no0_transposed.csv",
  "E:/SCD/数据/构建网络的数据/data_s2_no0_transposed.csv"
)

# 2. 转录因子
tf_vec <- c("BCL11A", "GATA1", "TAL1")

# 3. 输出目录
out_dir <- "E:/SCD/下游分析/零模型"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 4. 用于汇总所有文件的统计结果
all_summary <- list()

# 5. 逐个文件处理
for (file in expr_files) {
  cat("正在处理文件：", file, "\n")
  
  # 读取表达矩阵
  df <- fread(file, header = TRUE, data.table = FALSE)
  
  # 提取第一列基因名
  genes <- df[[1]]
  genes <- trimws(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- unique(genes)
  
  cat("提取到基因数：", length(genes), "\n")
  
  # 构造 TF-Target 对
  tf_target_df <- expand.grid(
    TF = tf_vec,
    Target = genes,
    stringsAsFactors = FALSE
  )
  
  # 生成输出文件名
  base_name <- tools::file_path_sans_ext(basename(file))
  out_file <- file.path(out_dir, paste0(base_name, "_TF_Target_pairs.tsv"))
  
  # 保存 TF-Target 对
  fwrite(tf_target_df, out_file, sep = "\t", quote = FALSE)
  
  # 统计每个 TF 的调控对数量
  tf_count_df <- as.data.frame(table(tf_target_df$TF), stringsAsFactors = FALSE)
  colnames(tf_count_df) <- c("TF", "Regulatory_Pair_Count")
  tf_count_df$File <- base_name
  tf_count_df <- tf_count_df[, c("File", "TF", "Regulatory_Pair_Count")]
  
  # 保存当前文件的统计结果
  count_file <- file.path(out_dir, paste0(base_name, "_TF_pair_count.tsv"))
  fwrite(tf_count_df, count_file, sep = "\t", quote = FALSE)
  
  # 加入总汇总
  all_summary[[base_name]] <- tf_count_df
  
  cat("生成调控对总数：", nrow(tf_target_df), "\n")
  cat("每个转录因子的调控对数量：\n")
  print(tf_count_df)
  cat("TF-Target对已保存：", out_file, "\n")
  cat("统计结果已保存：", count_file, "\n\n")
}

# 6. 保存所有文件的总汇总表
all_summary_df <- do.call(rbind, all_summary)
summary_file <- file.path(out_dir, "All_files_TF_pair_count_summary.tsv")
fwrite(all_summary_df, summary_file, sep = "\t", quote = FALSE)

cat("全部完成。\n")
cat("所有文件汇总统计已保存：", summary_file, "\n")

















# -----------------------------------------------------------------------------
# 4组数据 × 3个TF 随机采样100次，与ChIP金标准做重叠，并计算平均P值
# 注意：抽样使用的是 TF-Target 对文件，不是 TF_pair_count 文件
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

set.seed(123)  # 保证结果可重复

# -----------------------------
# 1. 输入文件
# -----------------------------
pair_files <- list(
  c1 = "E:/SCD/下游分析/零模型/data_c1_no0_transposed_TF_Target_pairs.tsv",
  c2 = "E:/SCD/下游分析/零模型/data_c2_no0_transposed_TF_Target_pairs.tsv",
  s1 = "E:/SCD/下游分析/零模型/data_s1_no0_transposed_TF_Target_pairs.tsv",
  s2 = "E:/SCD/下游分析/零模型/data_s2_no0_transposed_TF_Target_pairs.tsv"
)

chip_files <- list(
  BCL11A = "E:/SCD/数据/CHIP-seq数据/BCL11A_ChIP_targets.tsv",
  GATA1  = "E:/SCD/数据/CHIP-seq数据/GATA1_ChIP_targets.tsv",
  TAL1   = "E:/SCD/数据/CHIP-seq数据/TAL1_ChIP_targets.tsv"
)

# 每组总边数
total_genes_list <- list(
  c1 = 573029,
  c2 = 530811,
  s1 = 577357,
  s2 = 574243
)

# 每组每个TF的采样数
sample_sizes <- list(
  c1 = c(BCL11A = 928,  GATA1 = 767,  TAL1 = 629),
  c2 = c(BCL11A = 2764, GATA1 = 1227, TAL1 = 1694),
  s1 = c(BCL11A = 413,  GATA1 = 1257, TAL1 = 718),
  s2 = c(BCL11A = 1437, GATA1 = 2436, TAL1 = 437)
)

# 重复次数
n_iter <- 1000

# 输出目录
out_dir <- "E:/SCD/下游分析/零模型"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. 读取金标准靶基因
# -----------------------------
chip_targets <- list()

for (tf in names(chip_files)) {
  chip_df <- fread(chip_files[[tf]], header = TRUE, data.table = FALSE)
  
  # 默认取 Target 列
  if (!"Target" %in% colnames(chip_df)) {
    stop(paste("金标准文件缺少 Target 列：", chip_files[[tf]]))
  }
  
  chip_targets[[tf]] <- unique(trimws(as.character(chip_df$Target)))
}

# -----------------------------
# 3. 单次抽样 + Fisher检验
#    按你提供的计算方式执行
# -----------------------------
run_one_test <- function(pair_df, tf_name, sample_n, chip_target_set, total_genes) {
  
  tf_df <- pair_df %>%
    filter(TF == tf_name)
  
  if (nrow(tf_df) == 0) {
    return(data.frame(
      TF = tf_name,
      Sample_Size = sample_n,
      Overlap = NA_integer_,
      CHIP_Targets = length(chip_target_set),
      Predicted_Targets = sample_n,
      Total_Genes = total_genes,
      P_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  # 去重，避免重复Target
  tf_targets <- unique(trimws(as.character(tf_df$Target)))
  
  if (length(tf_targets) < sample_n) {
    warning(paste0(tf_name, " 可供抽样的靶基因数不足：可用 ", length(tf_targets), 
                   "，需要 ", sample_n))
    return(data.frame(
      TF = tf_name,
      Sample_Size = sample_n,
      Overlap = NA_integer_,
      CHIP_Targets = length(chip_target_set),
      Predicted_Targets = sample_n,
      Total_Genes = total_genes,
      P_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  # 随机抽样
  sampled_targets <- sample(tf_targets, size = sample_n, replace = FALSE)
  
  # 重叠数
  a <- length(intersect(sampled_targets, chip_target_set))
  b <- length(chip_target_set)
  c <- sample_n
  d <- total_genes
  
  # 按你给的方式构造列联表
  contingency_table <- matrix(c(a, b, c, d), nrow = 2, ncol = 2, byrow = TRUE)
  colnames(contingency_table) <- c("Overlap", "No_Overlap")
  rownames(contingency_table) <- c("Target_Genes", "Non_Target_Genes")
  
  fisher_result <- fisher.test(contingency_table)
  
  data.frame(
    TF = tf_name,
    Sample_Size = sample_n,
    Overlap = a,
    CHIP_Targets = b,
    Predicted_Targets = c,
    Total_Genes = d,
    P_value = fisher_result$p.value,
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# 4. 主循环
# -----------------------------
all_iter_results <- list()
summary_results <- list()

for (grp in names(pair_files)) {
  cat("正在处理组别：", grp, "\n")
  
  pair_df <- fread(pair_files[[grp]], header = TRUE, data.table = FALSE)
  
  if (!all(c("TF", "Target") %in% colnames(pair_df))) {
    stop(paste("文件缺少 TF 或 Target 列：", pair_files[[grp]]))
  }
  
  pair_df$TF <- trimws(as.character(pair_df$TF))
  pair_df$Target <- trimws(as.character(pair_df$Target))
  
  group_iter_res <- list()
  
  for (tf in c("BCL11A", "GATA1", "TAL1")) {
    cat("  TF：", tf, "\n")
    
    tf_iter_res <- lapply(1:n_iter, function(i) {
      res <- run_one_test(
        pair_df = pair_df,
        tf_name = tf,
        sample_n = sample_sizes[[grp]][tf],
        chip_target_set = chip_targets[[tf]],
        total_genes = total_genes_list[[grp]]
      )
      res$Group <- grp
      res$Iteration <- i
      res
    })
    
    tf_iter_res <- do.call(rbind, tf_iter_res)
    group_iter_res[[tf]] <- tf_iter_res
    
    # 计算平均P值
    summary_results[[paste(grp, tf, sep = "_")]] <- data.frame(
      Group = grp,
      TF = tf,
      Sample_Size = sample_sizes[[grp]][tf],
      Mean_P_value = mean(tf_iter_res$P_value, na.rm = TRUE),
      Mean_Overlap = mean(tf_iter_res$Overlap, na.rm = TRUE),
      CHIP_Targets = length(chip_targets[[tf]]),
      Total_Genes = total_genes_list[[grp]],
      stringsAsFactors = FALSE
    )
  }
  
  group_iter_res <- do.call(rbind, group_iter_res)
  all_iter_results[[grp]] <- group_iter_res
}

# -----------------------------
# 5. 合并结果并保存
# -----------------------------
all_iter_results_df <- do.call(rbind, all_iter_results)
summary_results_df <- do.call(rbind, summary_results)

# 保存每次迭代结果
fwrite(
  all_iter_results_df,
  file.path(out_dir, "Random_sampling_100times_all_iteration_results.tsv"),
  sep = "\t",
  quote = FALSE
)

# 保存12组平均P值结果
fwrite(
  summary_results_df,
  file.path(out_dir, "Random_sampling_100times_mean_pvalues.tsv"),
  sep = "\t",
  quote = FALSE
)

cat("\n全部完成。\n")
cat("每次迭代结果：", file.path(out_dir, "Random_sampling_100times_all_iteration_results.tsv"), "\n")
cat("12组平均P值结果：", file.path(out_dir, "Random_sampling_100times_mean_pvalues.tsv"), "\n")





















library(data.table)
library(dplyr)
library(openxlsx)

# 1. 读取100次随机抽样结果
random_file <- "E:/SCD/下游分析/零模型/Random_sampling_100times_all_iteration_results.tsv"
random_df <- fread(random_file, data.table = FALSE)

# 2. 录入真实预测网络的重叠数
observed_df <- data.frame(
  Group = c("c1","c1","c1",
            "c2","c2","c2",
            "s1","s1","s1",
            "s2","s2","s2"),
  TF = c("BCL11A","GATA1","TAL1",
         "BCL11A","GATA1","TAL1",
         "BCL11A","GATA1","TAL1",
         "BCL11A","GATA1","TAL1"),
  Observed_Overlap = c(267,162,159,
                       896,249,498,
                       169,274,189,
                       552,469,126),
  stringsAsFactors = FALSE
)

# 3. 计算随机分布统计量
random_summary <- random_df %>%
  group_by(Group, TF) %>%
  summarise(
    Random_Mean = mean(Overlap, na.rm = TRUE),
    Random_SD   = sd(Overlap, na.rm = TRUE),
    Random_Min  = min(Overlap, na.rm = TRUE),
    Random_Max  = max(Overlap, na.rm = TRUE),
    N = sum(!is.na(Overlap)),
    .groups = "drop"
  )

# 4. 合并真实值与随机统计
result_df <- observed_df %>%
  left_join(random_summary, by = c("Group", "TF"))

# 5. 计算经验p值、Z-score、Fold enrichment
empirical_p_fun <- function(obs, rand_vec) {
  rand_vec <- rand_vec[!is.na(rand_vec)]
  (sum(rand_vec >= obs) + 1) / (length(rand_vec) + 1)
}

result_list <- list()

for (i in 1:nrow(observed_df)) {
  grp <- observed_df$Group[i]
  tf  <- observed_df$TF[i]
  obs <- observed_df$Observed_Overlap[i]
  
  rand_vec <- random_df %>%
    filter(Group == grp, TF == tf) %>%
    pull(Overlap)
  
  rand_vec <- rand_vec[!is.na(rand_vec)]
  
  rand_mean <- mean(rand_vec)
  rand_sd   <- sd(rand_vec)
  
  emp_p <- empirical_p_fun(obs, rand_vec)
  z_val <- ifelse(rand_sd > 0, (obs - rand_mean) / rand_sd, NA)
  fold  <- ifelse(rand_mean > 0, obs / rand_mean, NA)
  
  result_list[[i]] <- data.frame(
    Group = grp,
    TF = tf,
    Observed_Overlap = obs,
    Random_Mean = rand_mean,
    Random_SD = rand_sd,
    Random_Min = min(rand_vec),
    Random_Max = max(rand_vec),
    Empirical_P = emp_p,
    Z_score = z_val,
    Fold_Enrichment = fold,
    stringsAsFactors = FALSE
  )
}

final_result <- do.call(rbind, result_list)

# 6. 多重检验校正
final_result$Empirical_FDR <- p.adjust(final_result$Empirical_P, method = "BH")

# 7. 保存为 TSV
out_tsv <- "E:/SCD/下游分析/零模型/Observed_vs_Random_significance.tsv"
fwrite(final_result, out_tsv, sep = "\t", quote = FALSE)

# 8. 保存为 Excel
out_xlsx <- "E:/SCD/下游分析/零模型/Observed_vs_Random_significance.xlsx"
write.xlsx(final_result, out_xlsx, rowNames = FALSE)

# 9. 输出结果
print(final_result)
cat("TSV结果已保存到：", out_tsv, "\n")
cat("Excel结果已保存到：", out_xlsx, "\n")


library(data.table)
library(dplyr)
library(openxlsx)
library(ggplot2)

# 1. 读取100次随机抽样结果
random_file <- "E:/SCD/下游分析/零模型/Random_sampling_100times_all_iteration_results.tsv"
random_df <- fread(random_file, data.table = FALSE)

# 2. 录入真实预测网络的重叠数
observed_df <- data.frame(
  Group = c("c1","c1","c1",
            "c2","c2","c2",
            "s1","s1","s1",
            "s2","s2","s2"),
  TF = c("BCL11A","GATA1","TAL1",
         "BCL11A","GATA1","TAL1",
         "BCL11A","GATA1","TAL1",
         "BCL11A","GATA1","TAL1"),
  Observed_Overlap = c(267,162,159,
                       896,249,498,
                       169,274,189,
                       552,469,126),
  stringsAsFactors = FALSE
)

# 3. 计算随机分布统计量
random_summary <- random_df %>%
  group_by(Group, TF) %>%
  summarise(
    Random_Mean = mean(Overlap, na.rm = TRUE),
    Random_SD   = sd(Overlap, na.rm = TRUE),
    Random_Min  = min(Overlap, na.rm = TRUE),
    Random_Max  = max(Overlap, na.rm = TRUE),
    N = sum(!is.na(Overlap)),
    .groups = "drop"
  )

# 4. 合并真实值与随机统计
final_result <- observed_df %>%
  left_join(random_summary, by = c("Group", "TF"))

# 5. 计算经验p值、Z-score、Fold enrichment
empirical_p_fun <- function(obs, rand_vec) {
  rand_vec <- rand_vec[!is.na(rand_vec)]
  (sum(rand_vec >= obs) + 1) / (length(rand_vec) + 1)
}

for (i in 1:nrow(final_result)) {
  grp <- final_result$Group[i]
  tf  <- final_result$TF[i]
  obs <- final_result$Observed_Overlap[i]
  
  rand_vec <- random_df %>%
    filter(Group == grp, TF == tf) %>%
    pull(Overlap)
  
  rand_vec <- rand_vec[!is.na(rand_vec)]
  
  rand_mean <- mean(rand_vec)
  rand_sd   <- sd(rand_vec)
  
  final_result$Empirical_P[i]     <- empirical_p_fun(obs, rand_vec)
  final_result$Z_score[i]         <- ifelse(rand_sd > 0, (obs - rand_mean) / rand_sd, NA)
  final_result$Fold_Enrichment[i] <- ifelse(rand_mean > 0, obs / rand_mean, NA)
}

# 6. 多重检验校正
final_result$Empirical_FDR <- p.adjust(final_result$Empirical_P, method = "BH")

# 7. 设置顺序
final_result$Group <- factor(final_result$Group, levels = c("c1", "c2", "s1", "s2"))
final_result$TF <- factor(final_result$TF, levels = c("BCL11A", "GATA1", "TAL1"))

# 组合标签
final_result$Label <- paste(final_result$Group, final_result$TF, sep = "_")
final_result$Label <- factor(
  final_result$Label,
  levels = c("c1_BCL11A","c1_GATA1","c1_TAL1",
             "c2_BCL11A","c2_GATA1","c2_TAL1",
             "s1_BCL11A","s1_GATA1","s1_TAL1",
             "s2_BCL11A","s2_GATA1","s2_TAL1")
)

# 增加数值型横坐标，便于画阴影带
final_result$x <- 1:nrow(final_result)

# 8. 保存为 TSV
out_tsv <- "E:/SCD/下游分析/零模型/Observed_vs_Random_significance.tsv"
fwrite(final_result, out_tsv, sep = "\t", quote = FALSE)

# 9. 保存为 Excel
out_xlsx <- "E:/SCD/下游分析/零模型/Observed_vs_Random_significance.xlsx"
write.xlsx(final_result, out_xlsx, rowNames = FALSE)

# 10. 绘制两条折线图：
# 一条 = 预测网络重叠数
# 一条 = 随机网络平均重叠数

p_line <- ggplot(final_result, aes(x = x)) +
  geom_line(aes(y = Observed_Overlap, color = "Predicted Network"), linewidth = 1) +
  geom_point(aes(y = Observed_Overlap, color = "Predicted Network"), size = 2.5) +
  geom_line(aes(y = Random_Mean, color = "Random Mean"), linewidth = 1) +
  geom_point(aes(y = Random_Mean, color = "Random Mean"), size = 2.5) +
  scale_x_continuous(
    breaks = final_result$x,
    labels = final_result$Label
  ) +
  scale_color_manual(
    values = c(
      "Predicted Network" = "C30D23",
      "Random Mean" = "1D2088"
    )
  ) +
  labs(
    title = "Observed overlap vs random mean overlap",
    x = "Group_TF",
    y = "Overlap count",
    color = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 11. 保存图片
out_line <- "E:/SCD/下游分析/零模型/Observed_vs_RandomMean_lineplot.png"
ggsave(out_line, p_line, width = 10, height = 5.5, dpi = 300)

# 12. 输出图片
print(p_line)
cat("折线图已保存到：", out_line, "\n")