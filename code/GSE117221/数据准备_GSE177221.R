# ================================
# GSE117221 完整流程（按平均表达量筛前10000个基因）
# 修正版：
# 1. 读取三个合并后的 counts 矩阵时，稳健剔除 summary 行
# 2. 先过滤低表达，再按 log2(CPM+1) 平均表达量筛前10000基因
# 3. 保存原始转置矩阵、top10000转置矩阵、top10000基因列表、平均表达量表、专用先验调控对
# 4. 从 top10000 转置矩阵恢复不转置矩阵
# 5. 构建三个数据专用 TF 列表
# 6. 基于不转置矩阵构建基因ID文件、专用先验数据集（CSV）、ID文件（TSV）
# ================================

# ---------- 1. 加载包 ----------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("edgeR", quietly = TRUE)) install.packages("edgeR")

library(data.table)
library(edgeR)

# ---------- 2. 路径设置 ----------
data_dir <- "E:/SCD/其他数据/GSE117221"

# 三个合并后的表达矩阵
file_H  <- file.path(data_dir, "GSE117221_H_counts_merged.txt")
file_TI <- file.path(data_dir, "GSE117221_TI_counts_merged.txt")
file_TM <- file.path(data_dir, "GSE117221_TM_counts_merged.txt")

# 第一套先验调控对文件（两列）
prior_file_1 <- "E:/SCD/数据/构建网络的数据/label_c1.csv"

# TF 列表
tf_file <- "E:/SCD/数据/构建网络的数据/human_tf_list.txt"

# 第二套先验信息文件（三列：Gene1, Gene2, Type）
prior_file_2 <- "C:/Users/Administrator/Desktop/new_GRN_Lung_GEN_counts_genename_c1.csv"

# ---------- 3. 定义 summary 行识别函数 ----------
# 比精确匹配更稳健：支持空格、大小写、下划线、双下划线、Unassigned 等各种格式
is_bad_feature <- function(x) {
  x2 <- trimws(as.character(x))
  x2_low <- tolower(x2)
  
  bad_exact <- c(
    "assigned",
    "alignment_not_unique",
    "no_feature",
    "ambiguous",
    "__no_feature",
    "__ambiguous",
    "__too_low_aqual",
    "__not_aligned",
    "__alignment_not_unique"
  )
  
  bad_pattern <- grepl("^unassigned_", x2_low) |
    grepl("^__", x2_low) |
    grepl("alignment[_ ]?not[_ ]?unique", x2_low) |
    grepl("no[_ ]?feature", x2_low) |
    grepl("^ambiguous$", x2_low)
  
  bad_exact_match <- x2_low %in% bad_exact
  
  return(bad_exact_match | bad_pattern)
}

# ---------- 4. 清洗表达矩阵函数 ----------
clean_expr_df <- function(df, dataset_name = "dataset") {
  colnames(df)[1] <- "Gene"
  
  # 转字符并去前后空格
  df$Gene <- trimws(as.character(df$Gene))
  
  # 去掉空值
  df <- df[!is.na(df$Gene) & df$Gene != "", , drop = FALSE]
  
  # 找出 summary 行
  bad_idx <- is_bad_feature(df$Gene)
  removed_genes <- unique(df$Gene[bad_idx])
  
  if (length(removed_genes) > 0) {
    cat("\n", dataset_name, "中移除的非基因行：\n", sep = "")
    print(removed_genes)
  } else {
    cat("\n", dataset_name, "中未检测到 summary 行。\n", sep = "")
  }
  
  # 去掉 summary 行
  df <- df[!bad_idx, , drop = FALSE]
  
  # 去重
  dup_genes <- unique(df$Gene[duplicated(df$Gene)])
  if (length(dup_genes) > 0) {
    cat(dataset_name, "中检测到重复基因，保留第一条，重复基因示例：\n")
    print(head(dup_genes, 20))
  }
  df <- df[!duplicated(df$Gene), , drop = FALSE]
  
  rownames(df) <- NULL
  return(df)
}

# =========================================================
# Part 1：按平均表达量筛前10000个基因，并保存相关文件
# =========================================================

# ---------- 5. 读取表达矩阵（修正版） ----------
read_expr_matrix <- function(file_path, dataset_name = "dataset") {
  df <- fread(file_path, data.table = FALSE)
  df <- clean_expr_df(df, dataset_name = dataset_name)
  
  rownames(df) <- df$Gene
  df$Gene <- NULL
  
  mat <- as.matrix(df)
  mode(mat) <- "numeric"
  
  return(mat)   # 行=gene, 列=sample
}

expr_H  <- read_expr_matrix(file_H,  dataset_name = "H merged counts")
expr_TI <- read_expr_matrix(file_TI, dataset_name = "TI merged counts")
expr_TM <- read_expr_matrix(file_TM, dataset_name = "TM merged counts")

# ---------- 6. 定义基因筛选函数 ----------
select_top_genes <- function(count_mat, top_n = 10000, cpm_cutoff = 1, min_prop = 0.2) {
  min_samples <- ceiling(ncol(count_mat) * min_prop)
  
  cpm_mat <- edgeR::cpm(count_mat)
  
  keep <- rowSums(cpm_mat > cpm_cutoff) >= min_samples
  filtered_mat <- count_mat[keep, , drop = FALSE]
  
  log_cpm <- log2(edgeR::cpm(filtered_mat) + 1)
  
  gene_mean <- rowMeans(log_cpm)
  
  top_n_use <- min(top_n, length(gene_mean))
  top_idx <- order(gene_mean, decreasing = TRUE)[1:top_n_use]
  top_genes <- rownames(log_cpm)[top_idx]
  
  return(list(
    filtered_count_mat = filtered_mat,
    log_cpm_mat = log_cpm,
    top_genes = top_genes,
    top_count_mat = filtered_mat[top_genes, , drop = FALSE],
    top_log_cpm_mat = log_cpm[top_genes, , drop = FALSE],
    gene_mean_table = data.frame(
      Gene = rownames(log_cpm),
      MeanExpression = gene_mean,
      stringsAsFactors = FALSE
    )[order(gene_mean, decreasing = TRUE), ]
  ))
}

res_H  <- select_top_genes(expr_H,  top_n = 10000)
res_TI <- select_top_genes(expr_TI, top_n = 10000)
res_TM <- select_top_genes(expr_TM, top_n = 10000)

# ---------- 7. 转置矩阵 ----------
transpose_expr <- function(mat) {
  tmat <- as.data.frame(t(mat), check.names = FALSE)
  tmat <- cbind(Sample = rownames(tmat), tmat)
  rownames(tmat) <- NULL
  return(tmat)
}

expr_H_t  <- transpose_expr(expr_H)
expr_TI_t <- transpose_expr(expr_TI)
expr_TM_t <- transpose_expr(expr_TM)

top_H_t  <- transpose_expr(res_H$top_count_mat)
top_TI_t <- transpose_expr(res_TI$top_count_mat)
top_TM_t <- transpose_expr(res_TM$top_count_mat)

# ---------- 8. 读取第一套先验调控对 ----------
prior_df_1 <- fread(prior_file_1, data.table = FALSE)
colnames(prior_df_1)[1:2] <- c("Gene1", "Gene2")
prior_df_1$Gene1 <- trimws(as.character(prior_df_1$Gene1))
prior_df_1$Gene2 <- trimws(as.character(prior_df_1$Gene2))
prior_df_1 <- prior_df_1[!is.na(prior_df_1$Gene1) & !is.na(prior_df_1$Gene2), , drop = FALSE]
prior_df_1 <- unique(prior_df_1)

# ---------- 9. 构建三个数据专用先验调控对 ----------
build_prior_for_dataset <- function(prior_df, gene_set) {
  prior_sub <- prior_df[
    prior_df$Gene1 %in% gene_set & prior_df$Gene2 %in% gene_set,
    ,
    drop = FALSE
  ]
  prior_sub <- unique(prior_sub)
  return(prior_sub)
}

prior_H_1  <- build_prior_for_dataset(prior_df_1, res_H$top_genes)
prior_TI_1 <- build_prior_for_dataset(prior_df_1, res_TI$top_genes)
prior_TM_1 <- build_prior_for_dataset(prior_df_1, res_TM$top_genes)

# ---------- 10. 保存结果 ----------
write.csv(expr_H_t,  file.path(data_dir, "GSE117221_H_transposed.csv"),  row.names = FALSE)
write.csv(expr_TI_t, file.path(data_dir, "GSE117221_TI_transposed.csv"), row.names = FALSE)
write.csv(expr_TM_t, file.path(data_dir, "GSE117221_TM_transposed.csv"), row.names = FALSE)

write.csv(top_H_t,  file.path(data_dir, "GSE117221_H_top10000_transposed.csv"),  row.names = FALSE)
write.csv(top_TI_t, file.path(data_dir, "GSE117221_TI_top10000_transposed.csv"), row.names = FALSE)
write.csv(top_TM_t, file.path(data_dir, "GSE117221_TM_top10000_transposed.csv"), row.names = FALSE)

write.csv(data.frame(Gene = res_H$top_genes),  file.path(data_dir, "GSE117221_H_top10000_genes.csv"),  row.names = FALSE)
write.csv(data.frame(Gene = res_TI$top_genes), file.path(data_dir, "GSE117221_TI_top10000_genes.csv"), row.names = FALSE)
write.csv(data.frame(Gene = res_TM$top_genes), file.path(data_dir, "GSE117221_TM_top10000_genes.csv"), row.names = FALSE)

write.csv(res_H$gene_mean_table,  file.path(data_dir, "GSE117221_H_gene_mean_expression.csv"),  row.names = FALSE)
write.csv(res_TI$gene_mean_table, file.path(data_dir, "GSE117221_TI_gene_mean_expression.csv"), row.names = FALSE)
write.csv(res_TM$gene_mean_table, file.path(data_dir, "GSE117221_TM_gene_mean_expression.csv"), row.names = FALSE)

write.csv(prior_H_1,  file.path(data_dir, "GSE117221_H_top10000_prior_pairs.csv"),  row.names = FALSE)
write.csv(prior_TI_1, file.path(data_dir, "GSE117221_TI_top10000_prior_pairs.csv"), row.names = FALSE)
write.csv(prior_TM_1, file.path(data_dir, "GSE117221_TM_top10000_prior_pairs.csv"), row.names = FALSE)

# ---------- 11. 输出结果概况 ----------
cat("\nPart 1 处理完成！\n\n")

cat("H组：\n")
cat("清理summary行后基因数 =", nrow(expr_H), "\n")
cat("低表达过滤后基因数 =", nrow(res_H$filtered_count_mat), "\n")
cat("Top基因数 =", length(res_H$top_genes), "\n")
cat("专用先验调控对数 =", nrow(prior_H_1), "\n\n")

cat("TI组：\n")
cat("清理summary行后基因数 =", nrow(expr_TI), "\n")
cat("低表达过滤后基因数 =", nrow(res_TI$filtered_count_mat), "\n")
cat("Top基因数 =", length(res_TI$top_genes), "\n")
cat("专用先验调控对数 =", nrow(prior_TI_1), "\n\n")

cat("TM组：\n")
cat("清理summary行后基因数 =", nrow(expr_TM), "\n")
cat("低表达过滤后基因数 =", nrow(res_TM$filtered_count_mat), "\n")
cat("Top基因数 =", length(res_TM$top_genes), "\n")
cat("专用先验调控对数 =", nrow(prior_TM_1), "\n\n")

cat("Part 1 所有文件已保存到：\n", data_dir, "\n\n")


# =========================================================
# Part 2：从 top10000_transposed.csv 恢复不转置矩阵，并构建专用 TF 列表
# =========================================================

# ---------- 12. 文件路径 ----------
file_H_t  <- file.path(data_dir, "GSE117221_H_top10000_transposed.csv")
file_TI_t <- file.path(data_dir, "GSE117221_TI_top10000_transposed.csv")
file_TM_t <- file.path(data_dir, "GSE117221_TM_top10000_transposed.csv")

# ---------- 13. 读取转置后的表达矩阵 ----------
read_transposed_expr <- function(file_path) {
  df <- fread(file_path, data.table = FALSE, check.names = FALSE)
  colnames(df)[1] <- "Sample"
  return(df)
}

expr_H_t2  <- read_transposed_expr(file_H_t)
expr_TI_t2 <- read_transposed_expr(file_TI_t)
expr_TM_t2 <- read_transposed_expr(file_TM_t)

# ---------- 14. 恢复为“不转置”格式 ----------
restore_non_transposed <- function(df_t) {
  sample_names <- df_t$Sample
  expr_only <- df_t[, -1, drop = FALSE]
  
  mat <- t(as.matrix(expr_only))
  colnames(mat) <- sample_names
  
  out_df <- data.frame(
    Gene = rownames(mat),
    mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  rownames(out_df) <- NULL
  return(out_df)
}

expr_H_notrans  <- restore_non_transposed(expr_H_t2)
expr_TI_notrans <- restore_non_transposed(expr_TI_t2)
expr_TM_notrans <- restore_non_transposed(expr_TM_t2)

# 再清洗一次，确保写出的不转置矩阵也绝对干净
expr_H_notrans  <- clean_expr_df(expr_H_notrans,  dataset_name = "H top10000 notransposed")
expr_TI_notrans <- clean_expr_df(expr_TI_notrans, dataset_name = "TI top10000 notransposed")
expr_TM_notrans <- clean_expr_df(expr_TM_notrans, dataset_name = "TM top10000 notransposed")

# ---------- 15. 保存恢复后的不转置矩阵 ----------
write.csv(expr_H_notrans,
          file.path(data_dir, "GSE117221_H_top10000_notransposed.csv"),
          row.names = FALSE)

write.csv(expr_TI_notrans,
          file.path(data_dir, "GSE117221_TI_top10000_notransposed.csv"),
          row.names = FALSE)

write.csv(expr_TM_notrans,
          file.path(data_dir, "GSE117221_TM_top10000_notransposed.csv"),
          row.names = FALSE)

# ---------- 16. 读取 TF 列表 ----------
tf_df <- fread(tf_file, data.table = FALSE)
colnames(tf_df)[1] <- "TF"
tf_df$TF <- trimws(as.character(tf_df$TF))
tf_df <- tf_df[!is.na(tf_df$TF) & tf_df$TF != "", , drop = FALSE]
tf_df <- unique(tf_df)
tf_list <- tf_df$TF

# ---------- 17. 提取三个数据各自的 top10000 基因 ----------
genes_H  <- expr_H_notrans$Gene
genes_TI <- expr_TI_notrans$Gene
genes_TM <- expr_TM_notrans$Gene

# ---------- 18. 构建三个数据专用 TF 列表 ----------
tf_H  <- data.frame(TF = intersect(tf_list, genes_H),  stringsAsFactors = FALSE)
tf_TI <- data.frame(TF = intersect(tf_list, genes_TI), stringsAsFactors = FALSE)
tf_TM <- data.frame(TF = intersect(tf_list, genes_TM), stringsAsFactors = FALSE)

tf_H  <- tf_H[order(tf_H$TF), , drop = FALSE]
tf_TI <- tf_TI[order(tf_TI$TF), , drop = FALSE]
tf_TM <- tf_TM[order(tf_TM$TF), , drop = FALSE]

# ---------- 19. 保存三个专用 TF 列表 ----------
write.table(tf_H,
            file.path(data_dir, "GSE117221_H_top10000_TF_list.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(tf_TI,
            file.path(data_dir, "GSE117221_TI_top10000_TF_list.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(tf_TM,
            file.path(data_dir, "GSE117221_TM_top10000_TF_list.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ---------- 20. 输出结果统计 ----------
cat("Part 2 处理完成！\n\n")

cat("恢复后的不转置表达矩阵：\n")
cat("H  :", dim(expr_H_notrans)[1],  "基因 x", dim(expr_H_notrans)[2] - 1, "样本\n")
cat("TI :", dim(expr_TI_notrans)[1], "基因 x", dim(expr_TI_notrans)[2] - 1, "样本\n")
cat("TM :", dim(expr_TM_notrans)[1], "基因 x", dim(expr_TM_notrans)[2] - 1, "样本\n\n")

cat("专用 TF 数量：\n")
cat("H  :", nrow(tf_H), "\n")
cat("TI :", nrow(tf_TI), "\n")
cat("TM :", nrow(tf_TM), "\n\n")

cat("Part 2 输出文件位置：\n", data_dir, "\n\n")


# =========================================================
# Part 3：基于 top10000_notransposed 构建基因ID文件、专用先验数据集、ID文件
# =========================================================

# ---------- 21. 三个不转置表达矩阵文件路径 ----------
file_H_notrans  <- file.path(data_dir, "GSE117221_H_top10000_notransposed.csv")
file_TI_notrans <- file.path(data_dir, "GSE117221_TI_top10000_notransposed.csv")
file_TM_notrans <- file.path(data_dir, "GSE117221_TM_top10000_notransposed.csv")

# ---------- 22. 读取表达矩阵（修正版） ----------
read_expr_file <- function(file_path, dataset_name = "dataset") {
  df <- fread(file_path, data.table = FALSE, check.names = FALSE)
  df <- clean_expr_df(df, dataset_name = dataset_name)
  return(df)
}

expr_H3  <- read_expr_file(file_H_notrans,  dataset_name = "H final notransposed")
expr_TI3 <- read_expr_file(file_TI_notrans, dataset_name = "TI final notransposed")
expr_TM3 <- read_expr_file(file_TM_notrans, dataset_name = "TM final notransposed")

# ---------- 23. 构建基因ID文件 ----------
build_gene_id_table <- function(expr_df) {
  gene_df <- data.frame(
    name = expr_df$Gene,
    ids = 0:(nrow(expr_df) - 1),
    stringsAsFactors = FALSE
  )
  return(gene_df)
}

gene_id_H  <- build_gene_id_table(expr_H3)
gene_id_TI <- build_gene_id_table(expr_TI3)
gene_id_TM <- build_gene_id_table(expr_TM3)

write.table(gene_id_H,
            file.path(data_dir, "GSE117221_H_gene_ids.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(gene_id_TI,
            file.path(data_dir, "GSE117221_TI_gene_ids.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(gene_id_TM,
            file.path(data_dir, "GSE117221_TM_gene_ids.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ---------- 24. 读取第二套先验文件 ----------
prior_df_2 <- fread(prior_file_2, data.table = FALSE, header = FALSE)

if (ncol(prior_df_2) < 3) {
  stop("先验文件列数少于3列，请检查文件格式。")
}

prior_df_2 <- prior_df_2[, 1:3, drop = FALSE]
colnames(prior_df_2) <- c("Gene1", "Gene2", "Type")

prior_df_2$Gene1 <- trimws(as.character(prior_df_2$Gene1))
prior_df_2$Gene2 <- trimws(as.character(prior_df_2$Gene2))
prior_df_2$Type  <- trimws(as.character(prior_df_2$Type))

prior_df_2 <- prior_df_2[
  !is.na(prior_df_2$Gene1) & prior_df_2$Gene1 != "" &
    !is.na(prior_df_2$Gene2) & prior_df_2$Gene2 != "",
  ,
  drop = FALSE
]

prior_df_2 <- unique(prior_df_2)

# ---------- 25. 构建三个数据专用先验数据集 ----------
build_dataset_prior <- function(prior_df, gene_vector) {
  out <- prior_df[
    prior_df$Gene1 %in% gene_vector & prior_df$Gene2 %in% gene_vector,
    ,
    drop = FALSE
  ]
  out <- unique(out)
  rownames(out) <- NULL
  return(out)
}

prior_H_2  <- build_dataset_prior(prior_df_2, expr_H3$Gene)
prior_TI_2 <- build_dataset_prior(prior_df_2, expr_TI3$Gene)
prior_TM_2 <- build_dataset_prior(prior_df_2, expr_TM3$Gene)

write.csv(prior_H_2,
          file.path(data_dir, "GSE117221_H_prior_dataset.csv"),
          row.names = FALSE)

write.csv(prior_TI_2,
          file.path(data_dir, "GSE117221_TI_prior_dataset.csv"),
          row.names = FALSE)

write.csv(prior_TM_2,
          file.path(data_dir, "GSE117221_TM_prior_dataset.csv"),
          row.names = FALSE)

# ---------- 26. 将专用先验数据集映射为基因ID对 ----------
convert_prior_to_ids <- function(prior_sub, gene_id_df) {
  gene_to_id <- setNames(gene_id_df$ids, gene_id_df$name)
  
  out_id <- data.frame(
    Gene1_ID = unname(gene_to_id[prior_sub$Gene1]),
    Gene2_ID = unname(gene_to_id[prior_sub$Gene2]),
    stringsAsFactors = FALSE
  )
  
  out_id <- out_id[
    !is.na(out_id$Gene1_ID) & !is.na(out_id$Gene2_ID),
    ,
    drop = FALSE
  ]
  
  rownames(out_id) <- NULL
  return(out_id)
}

prior_H_ids  <- convert_prior_to_ids(prior_H_2,  gene_id_H)
prior_TI_ids <- convert_prior_to_ids(prior_TI_2, gene_id_TI)
prior_TM_ids <- convert_prior_to_ids(prior_TM_2, gene_id_TM)

write.table(prior_H_ids,
            file.path(data_dir, "GSE117221_H_prior_ids.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(prior_TI_ids,
            file.path(data_dir, "GSE117221_TI_prior_ids.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(prior_TM_ids,
            file.path(data_dir, "GSE117221_TM_prior_ids.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ---------- 27. 生成带Type的ID版先验文件 ----------
convert_prior_to_ids_with_type <- function(prior_sub, gene_id_df) {
  gene_to_id <- setNames(gene_id_df$ids, gene_id_df$name)
  
  out2 <- data.frame(
    Gene1_ID = unname(gene_to_id[prior_sub$Gene1]),
    Gene2_ID = unname(gene_to_id[prior_sub$Gene2]),
    Type = prior_sub$Type,
    stringsAsFactors = FALSE
  )
  
  out2 <- out2[
    !is.na(out2$Gene1_ID) & !is.na(out2$Gene2_ID),
    ,
    drop = FALSE
  ]
  
  rownames(out2) <- NULL
  return(out2)
}

prior_H_ids_type  <- convert_prior_to_ids_with_type(prior_H_2,  gene_id_H)
prior_TI_ids_type <- convert_prior_to_ids_with_type(prior_TI_2, gene_id_TI)
prior_TM_ids_type <- convert_prior_to_ids_with_type(prior_TM_2, gene_id_TM)

write.table(prior_H_ids_type,
            file.path(data_dir, "GSE117221_H_prior_ids_with_type.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(prior_TI_ids_type,
            file.path(data_dir, "GSE117221_TI_prior_ids_with_type.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(prior_TM_ids_type,
            file.path(data_dir, "GSE117221_TM_prior_ids_with_type.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ---------- 28. 核对文件 ----------
build_check_table <- function(prior_sub, gene_id_df) {
  gene_to_id <- setNames(gene_id_df$ids, gene_id_df$name)
  
  check_df <- data.frame(
    Gene1 = prior_sub$Gene1,
    Gene2 = prior_sub$Gene2,
    Type = prior_sub$Type,
    Gene1_ID = unname(gene_to_id[prior_sub$Gene1]),
    Gene2_ID = unname(gene_to_id[prior_sub$Gene2]),
    stringsAsFactors = FALSE
  )
  
  check_df <- check_df[
    !is.na(check_df$Gene1_ID) & !is.na(check_df$Gene2_ID),
    ,
    drop = FALSE
  ]
  
  rownames(check_df) <- NULL
  return(check_df)
}

check_H  <- build_check_table(prior_H_2,  gene_id_H)
check_TI <- build_check_table(prior_TI_2, gene_id_TI)
check_TM <- build_check_table(prior_TM_2, gene_id_TM)

write.table(check_H,
            file.path(data_dir, "GSE117221_H_prior_check_table.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(check_TI,
            file.path(data_dir, "GSE117221_TI_prior_check_table.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(check_TM,
            file.path(data_dir, "GSE117221_TM_prior_check_table.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ---------- 29. 输出结果统计 ----------
cat("Part 3 处理完成！\n\n")

cat("三个基因ID文件（TXT）：\n")
cat("H  :", nrow(gene_id_H),  "个基因\n")
cat("TI :", nrow(gene_id_TI), "个基因\n")
cat("TM :", nrow(gene_id_TM), "个基因\n\n")

cat("三个专用先验数据集（CSV）：\n")
cat("H  :", nrow(prior_H_2),  "条\n")
cat("TI :", nrow(prior_TI_2), "条\n")
cat("TM :", nrow(prior_TM_2), "条\n\n")

cat("三个ID版先验文件（TSV）：\n")
cat("H  :", nrow(prior_H_ids),  "条\n")
cat("TI :", nrow(prior_TI_ids), "条\n")
cat("TM :", nrow(prior_TM_ids), "条\n\n")

cat("所有结果已保存到：\n", data_dir, "\n")