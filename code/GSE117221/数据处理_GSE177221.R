# ================================
# GSE117221 单样本 counts 文件合并为 3 组矩阵
# H / TI / TM
# ================================

# 1. 设置工作目录
data_dir <- "E:/SCD/其他数据/GSE117221"

# 2. 获取所有 counts 文件
files <- list.files(
  path = data_dir,
  pattern = "_counts\\.txt$",
  full.names = TRUE
)

# 查看文件数
cat("共检测到文件数：", length(files), "\n")

# 3. 定义一个函数：读取单个 counts 文件
# 假设每个文件两列：gene 和 count，没有表头
read_count_file <- function(file_path) {
  df <- read.table(
    file_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  
  # 若文件不是tab分隔，也可能被读成1列，这里尝试自动处理
  if (ncol(df) == 1) {
    df <- read.table(
      file_path,
      header = FALSE,
      sep = "",
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )
  }
  
  if (ncol(df) < 2) {
    stop(paste("文件格式有问题，无法读取两列数据：", file_path))
  }
  
  df <- df[, 1:2]
  colnames(df) <- c("Gene", "Count")
  
  # 样本名：去掉路径和后缀
  sample_name <- basename(file_path)
  sample_name <- sub("_counts\\.txt$", "", sample_name)
  
  colnames(df)[2] <- sample_name
  return(df)
}

# 4. 读取所有文件
count_list <- lapply(files, read_count_file)

# 5. 提取样本名
sample_names <- sapply(count_list, function(x) colnames(x)[2])

# 6. 根据文件名分组
# 注意顺序：先判 TM，再判 TI/TITD，再判 H，避免误匹配
group_labels <- ifelse(
  grepl("-TM_", sample_names),
  "TM",
  ifelse(
    grepl("-TI_|-TITD_", sample_names),
    "TI",
    ifelse(
      grepl("-H_", sample_names),
      "H",
      NA
    )
  )
)

# 查看哪些样本没被成功分组
group_info <- data.frame(
  sample = sample_names,
  group = group_labels,
  stringsAsFactors = FALSE
)

cat("样本分组情况：\n")
print(group_info)

if (any(is.na(group_labels))) {
  cat("\n以下样本未成功分组，请检查文件名：\n")
  print(group_info[is.na(group_labels), ])
}

# 7. 编写合并函数
merge_counts <- function(df_list) {
  merged_df <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), df_list)
  
  # NA 填 0
  merged_df[is.na(merged_df)] <- 0
  
  # 保证 count 为数值型
  for (i in 2:ncol(merged_df)) {
    merged_df[[i]] <- as.numeric(merged_df[[i]])
  }
  
  return(merged_df)
}

# 8. 分组提取
H_list  <- count_list[group_labels == "H"]
TI_list <- count_list[group_labels == "TI"]
TM_list <- count_list[group_labels == "TM"]

# 9. 合并
H_merged  <- merge_counts(H_list)
TI_merged <- merge_counts(TI_list)
TM_merged <- merge_counts(TM_list)

# 10. 可选：按基因名排序
H_merged  <- H_merged[order(H_merged$Gene), ]
TI_merged <- TI_merged[order(TI_merged$Gene), ]
TM_merged <- TM_merged[order(TM_merged$Gene), ]

# 11. 保存结果
write.table(
  H_merged,
  file = file.path(data_dir, "GSE117221_H_counts_merged.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  TI_merged,
  file = file.path(data_dir, "GSE117221_TI_counts_merged.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  TM_merged,
  file = file.path(data_dir, "GSE117221_TM_counts_merged.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\n合并完成！输出文件为：\n")
cat(file.path(data_dir, "GSE117221_H_counts_merged.txt"), "\n")
cat(file.path(data_dir, "GSE117221_TI_counts_merged.txt"), "\n")
cat(file.path(data_dir, "GSE117221_TM_counts_merged.txt"), "\n")

# 12. 查看维度
cat("\n各组矩阵维度：\n")
cat("H 组：", dim(H_merged)[1], "基因 x", dim(H_merged)[2] - 1, "样本\n")
cat("TI组：", dim(TI_merged)[1], "基因 x", dim(TI_merged)[2] - 1, "样本\n")
cat("TM组：", dim(TM_merged)[1], "基因 x", dim(TM_merged)[2] - 1, "样本\n")

















# ================================
# GSE117221 三个表达矩阵：
# 1. 转置
# 2. 筛选前5000基因
# 3. 基于筛选基因构建三个数据专用先验调控对
# 4. 所有结果保存为CSV
# ================================

# ---------- 1. 加载包 ----------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("edgeR", quietly = TRUE)) install.packages("edgeR")
if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")

library(data.table)
library(edgeR)
library(matrixStats)

# ---------- 2. 路径设置 ----------
data_dir <- "E:/SCD/其他数据/GSE117221"

# 三个合并后的表达矩阵
file_H  <- file.path(data_dir, "GSE117221_H_counts_merged.txt")
file_TI <- file.path(data_dir, "GSE117221_TI_counts_merged.txt")
file_TM <- file.path(data_dir, "GSE117221_TM_counts_merged.txt")

# 先验调控对文件
prior_file <- "E:/SCD/数据/构建网络的数据/label_c1.csv"

# ---------- 3. 读取表达矩阵 ----------
read_expr_matrix <- function(file_path) {
  df <- fread(file_path, data.table = FALSE)
  
  # 第一列改名为 Gene
  colnames(df)[1] <- "Gene"
  
  # 去重：如果基因名重复，保留第一条
  df <- df[!duplicated(df$Gene), ]
  
  # 设置行为基因名
  rownames(df) <- df$Gene
  df$Gene <- NULL
  
  # 转为数值矩阵
  mat <- as.matrix(df)
  mode(mat) <- "numeric"
  
  return(mat)   # 行=gene, 列=sample
}

expr_H  <- read_expr_matrix(file_H)
expr_TI <- read_expr_matrix(file_TI)
expr_TM <- read_expr_matrix(file_TM)

# ---------- 4. 定义基因筛选函数 ----------
# 思路：
# (1) 先按 CPM 过滤低表达
# (2) 再按 log2(CPM+1) 方差筛前5000
select_top_genes <- function(count_mat, top_n = 5000, cpm_cutoff = 1, min_prop = 0.2) {
  # count_mat: 行=gene, 列=sample
  
  # 计算每个基因至少在多少个样本中表达
  min_samples <- ceiling(ncol(count_mat) * min_prop)
  
  # CPM
  cpm_mat <- edgeR::cpm(count_mat)
  
  # 低表达过滤
  keep <- rowSums(cpm_mat > cpm_cutoff) >= min_samples
  filtered_mat <- count_mat[keep, , drop = FALSE]
  
  # 再计算 log2(CPM+1)
  log_cpm <- log2(edgeR::cpm(filtered_mat) + 1)
  
  # 计算方差
  gene_var <- matrixStats::rowVars(log_cpm)
  
  # 若过滤后不足5000个，则全部保留
  top_n_use <- min(top_n, length(gene_var))
  
  # 取方差最大的前 top_n_use 个基因
  top_idx <- order(gene_var, decreasing = TRUE)[1:top_n_use]
  top_genes <- rownames(log_cpm)[top_idx]
  
  # 返回结果
  return(list(
    filtered_count_mat = filtered_mat,
    log_cpm_mat = log_cpm,
    top_genes = top_genes,
    top_count_mat = filtered_mat[top_genes, , drop = FALSE],
    top_log_cpm_mat = log_cpm[top_genes, , drop = FALSE],
    gene_variance = data.frame(
      Gene = rownames(log_cpm),
      Variance = gene_var,
      stringsAsFactors = FALSE
    )[order(gene_var, decreasing = TRUE), ]
  ))
}

res_H  <- select_top_genes(expr_H,  top_n = 5000)
res_TI <- select_top_genes(expr_TI, top_n = 5000)
res_TM <- select_top_genes(expr_TM, top_n = 5000)

# ---------- 5. 转置矩阵 ----------
# 原始矩阵转置后：行=sample, 列=gene
transpose_expr <- function(mat) {
  tmat <- as.data.frame(t(mat), check.names = FALSE)
  tmat <- cbind(Sample = rownames(tmat), tmat)
  rownames(tmat) <- NULL
  return(tmat)
}

# 原始转置
expr_H_t  <- transpose_expr(expr_H)
expr_TI_t <- transpose_expr(expr_TI)
expr_TM_t <- transpose_expr(expr_TM)

# top5000 转置（这里用 top_count_mat；如果你更想用 logCPM，也可以换成 top_log_cpm_mat）
top_H_t  <- transpose_expr(res_H$top_count_mat)
top_TI_t <- transpose_expr(res_TI$top_count_mat)
top_TM_t <- transpose_expr(res_TM$top_count_mat)

# ---------- 6. 读取先验调控对 ----------
# 你的文件虽然叫 csv，但示例看起来像制表符分隔，fread 会自动识别
prior_df <- fread(prior_file, data.table = FALSE)

# 统一列名
colnames(prior_df)[1:2] <- c("Gene1", "Gene2")

# 去掉缺失
prior_df <- prior_df[!is.na(prior_df$Gene1) & !is.na(prior_df$Gene2), ]

# 去重
prior_df <- unique(prior_df)

# ---------- 7. 构建三个数据专用先验调控对 ----------
build_prior_for_dataset <- function(prior_df, gene_set) {
  prior_sub <- prior_df[
    prior_df$Gene1 %in% gene_set & prior_df$Gene2 %in% gene_set,
    ,
    drop = FALSE
  ]
  prior_sub <- unique(prior_sub)
  return(prior_sub)
}

prior_H  <- build_prior_for_dataset(prior_df, res_H$top_genes)
prior_TI <- build_prior_for_dataset(prior_df, res_TI$top_genes)
prior_TM <- build_prior_for_dataset(prior_df, res_TM$top_genes)

# ---------- 8. 保存结果 ----------
# 8.1 原始转置矩阵
write.csv(expr_H_t,  file.path(data_dir, "GSE117221_H_transposed.csv"),  row.names = FALSE)
write.csv(expr_TI_t, file.path(data_dir, "GSE117221_TI_transposed.csv"), row.names = FALSE)
write.csv(expr_TM_t, file.path(data_dir, "GSE117221_TM_transposed.csv"), row.names = FALSE)

# 8.2 top5000 转置矩阵
write.csv(top_H_t,  file.path(data_dir, "GSE117221_H_top5000_transposed.csv"),  row.names = FALSE)
write.csv(top_TI_t, file.path(data_dir, "GSE117221_TI_top5000_transposed.csv"), row.names = FALSE)
write.csv(top_TM_t, file.path(data_dir, "GSE117221_TM_top5000_transposed.csv"), row.names = FALSE)

# 8.3 top5000基因列表
write.csv(data.frame(Gene = res_H$top_genes),  file.path(data_dir, "GSE117221_H_top5000_genes.csv"),  row.names = FALSE)
write.csv(data.frame(Gene = res_TI$top_genes), file.path(data_dir, "GSE117221_TI_top5000_genes.csv"), row.names = FALSE)
write.csv(data.frame(Gene = res_TM$top_genes), file.path(data_dir, "GSE117221_TM_top5000_genes.csv"), row.names = FALSE)

# 8.4 方差表
write.csv(res_H$gene_variance,  file.path(data_dir, "GSE117221_H_gene_variance.csv"),  row.names = FALSE)
write.csv(res_TI$gene_variance, file.path(data_dir, "GSE117221_TI_gene_variance.csv"), row.names = FALSE)
write.csv(res_TM$gene_variance, file.path(data_dir, "GSE117221_TM_gene_variance.csv"), row.names = FALSE)

# 8.5 三个数据专用先验调控对
write.csv(prior_H,  file.path(data_dir, "GSE117221_H_top5000_prior_pairs.csv"),  row.names = FALSE)
write.csv(prior_TI, file.path(data_dir, "GSE117221_TI_top5000_prior_pairs.csv"), row.names = FALSE)
write.csv(prior_TM, file.path(data_dir, "GSE117221_TM_top5000_prior_pairs.csv"), row.names = FALSE)

# ---------- 9. 输出结果概况 ----------
cat("处理完成！\n\n")

cat("H组：\n")
cat("原始基因数 =", nrow(expr_H), "\n")
cat("低表达过滤后基因数 =", nrow(res_H$filtered_count_mat), "\n")
cat("Top基因数 =", length(res_H$top_genes), "\n")
cat("专用先验调控对数 =", nrow(prior_H), "\n\n")

cat("TI组：\n")
cat("原始基因数 =", nrow(expr_TI), "\n")
cat("低表达过滤后基因数 =", nrow(res_TI$filtered_count_mat), "\n")
cat("Top基因数 =", length(res_TI$top_genes), "\n")
cat("专用先验调控对数 =", nrow(prior_TI), "\n\n")

cat("TM组：\n")
cat("原始基因数 =", nrow(expr_TM), "\n")
cat("低表达过滤后基因数 =", nrow(res_TM$filtered_count_mat), "\n")
cat("Top基因数 =", length(res_TM$top_genes), "\n")
cat("专用先验调控对数 =", nrow(prior_TM), "\n\n")

cat("所有文件已保存到：\n", data_dir, "\n")














# ================================
# 目的：
# 1. 从三个 top5000_transposed.csv 中恢复不转置矩阵
# 2. 为 H / TI / TM 三个 top5000 数据分别构建专用 TF 列表
# 3. TF 列表保存为 TXT
# ================================

# ---------- 1. 加载包 ----------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)

# ---------- 2. 文件路径 ----------
data_dir <- "E:/SCD/其他数据/GSE117221"

file_H_t  <- file.path(data_dir, "GSE117221_H_top5000_transposed.csv")
file_TI_t <- file.path(data_dir, "GSE117221_TI_top5000_transposed.csv")
file_TM_t <- file.path(data_dir, "GSE117221_TM_top5000_transposed.csv")

tf_file <- "E:/SCD/数据/构建网络的数据/human_tf_list.txt"

# ---------- 3. 读取转置后的表达矩阵 ----------
read_transposed_expr <- function(file_path) {
  df <- fread(file_path, data.table = FALSE, check.names = FALSE)
  colnames(df)[1] <- "Sample"
  return(df)
}

expr_H_t  <- read_transposed_expr(file_H_t)
expr_TI_t <- read_transposed_expr(file_TI_t)
expr_TM_t <- read_transposed_expr(file_TM_t)

# ---------- 4. 恢复为“不转置”格式 ----------
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

expr_H_notrans  <- restore_non_transposed(expr_H_t)
expr_TI_notrans <- restore_non_transposed(expr_TI_t)
expr_TM_notrans <- restore_non_transposed(expr_TM_t)

# ---------- 5. 保存恢复后的不转置矩阵 ----------
write.csv(expr_H_notrans,
          file.path(data_dir, "GSE117221_H_top5000_notransposed.csv"),
          row.names = FALSE)

write.csv(expr_TI_notrans,
          file.path(data_dir, "GSE117221_TI_top5000_notransposed.csv"),
          row.names = FALSE)

write.csv(expr_TM_notrans,
          file.path(data_dir, "GSE117221_TM_top5000_notransposed.csv"),
          row.names = FALSE)

# ---------- 6. 读取 TF 列表 ----------
tf_df <- fread(tf_file, data.table = FALSE)
colnames(tf_df)[1] <- "TF"

tf_df <- tf_df[!is.na(tf_df$TF) & tf_df$TF != "", , drop = FALSE]
tf_df <- unique(tf_df)

tf_list <- tf_df$TF

# ---------- 7. 提取三个数据各自的 top5000 基因 ----------
genes_H  <- expr_H_notrans$Gene
genes_TI <- expr_TI_notrans$Gene
genes_TM <- expr_TM_notrans$Gene

# ---------- 8. 构建三个数据专用 TF 列表 ----------
tf_H  <- data.frame(TF = intersect(tf_list, genes_H),  stringsAsFactors = FALSE)
tf_TI <- data.frame(TF = intersect(tf_list, genes_TI), stringsAsFactors = FALSE)
tf_TM <- data.frame(TF = intersect(tf_list, genes_TM), stringsAsFactors = FALSE)

tf_H  <- tf_H[order(tf_H$TF), , drop = FALSE]
tf_TI <- tf_TI[order(tf_TI$TF), , drop = FALSE]
tf_TM <- tf_TM[order(tf_TM$TF), , drop = FALSE]

# ---------- 9. 保存三个专用 TF 列表（TXT） ----------
write.table(tf_H,
            file.path(data_dir, "GSE117221_H_top5000_TF_list.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

write.table(tf_TI,
            file.path(data_dir, "GSE117221_TI_top5000_TF_list.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

write.table(tf_TM,
            file.path(data_dir, "GSE117221_TM_top5000_TF_list.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

# ---------- 10. 输出结果统计 ----------
cat("处理完成！\n\n")

cat("恢复后的不转置表达矩阵：\n")
cat("H  :", dim(expr_H_notrans)[1],  "基因 x", dim(expr_H_notrans)[2] - 1,  "样本\n")
cat("TI :", dim(expr_TI_notrans)[1], "基因 x", dim(expr_TI_notrans)[2] - 1, "样本\n")
cat("TM :", dim(expr_TM_notrans)[1], "基因 x", dim(expr_TM_notrans)[2] - 1, "样本\n\n")

cat("专用 TF 数量：\n")
cat("H  :", nrow(tf_H),  "\n")
cat("TI :", nrow(tf_TI), "\n")
cat("TM :", nrow(tf_TM), "\n\n")

cat("输出文件位置：\n", data_dir, "\n")









# ================================
# 目的：
# 1. 基于三个不转置表达矩阵构建独立的基因ID文件
# 2. 基于新的三列先验文件构建三个数据专用先验数据集（CSV）
# 3. 将专用先验数据集稳定地转换为基因ID对（TSV）
# ================================

# ---------- 1. 加载包 ----------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)

# ---------- 2. 文件路径 ----------
data_dir <- "E:/SCD/其他数据/GSE117221"

# 三个不转置表达矩阵（第一列为Gene，后面为样本）
file_H  <- file.path(data_dir, "GSE117221_H_top5000_notransposed.csv")
file_TI <- file.path(data_dir, "GSE117221_TI_top5000_notransposed.csv")
file_TM <- file.path(data_dir, "GSE117221_TM_top5000_notransposed.csv")

# 新先验信息文件（三列：TF, Target, Type）
prior_file <- "C:/Users/Administrator/Desktop/new_GRN_Lung_GEN_counts_genename_c1.csv"

# ---------- 3. 读取表达矩阵 ----------
read_expr_file <- function(file_path) {
  df <- fread(file_path, data.table = FALSE, check.names = FALSE)
  colnames(df)[1] <- "Gene"
  
  # 去掉缺失和重复基因
  df <- df[!is.na(df$Gene) & df$Gene != "", , drop = FALSE]
  df <- df[!duplicated(df$Gene), , drop = FALSE]
  
  return(df)
}

expr_H  <- read_expr_file(file_H)
expr_TI <- read_expr_file(file_TI)
expr_TM <- read_expr_file(file_TM)

# ---------- 4. 构建基因ID文件 ----------
# 按表达矩阵中的基因顺序，从0开始编号
build_gene_id_table <- function(expr_df) {
  gene_df <- data.frame(
    name = expr_df$Gene,
    ids = 0:(nrow(expr_df) - 1),
    stringsAsFactors = FALSE
  )
  return(gene_df)
}

gene_id_H  <- build_gene_id_table(expr_H)
gene_id_TI <- build_gene_id_table(expr_TI)
gene_id_TM <- build_gene_id_table(expr_TM)

# 保存基因ID文件（TXT）
write.table(gene_id_H,
            file.path(data_dir, "GSE117221_H_gene_ids.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

write.table(gene_id_TI,
            file.path(data_dir, "GSE117221_TI_gene_ids.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

write.table(gene_id_TM,
            file.path(data_dir, "GSE117221_TM_gene_ids.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

# ---------- 5. 读取新的先验文件 ----------
prior_df <- fread(prior_file, data.table = FALSE, header = FALSE)

if (ncol(prior_df) < 3) {
  stop("先验文件列数少于3列，请检查文件格式。")
}

prior_df <- prior_df[, 1:3, drop = FALSE]
colnames(prior_df) <- c("Gene1", "Gene2", "Type")

# 去掉缺失
prior_df <- prior_df[
  !is.na(prior_df$Gene1) & prior_df$Gene1 != "" &
    !is.na(prior_df$Gene2) & prior_df$Gene2 != "",
  ,
  drop = FALSE
]

# 去重
prior_df <- unique(prior_df)

# ---------- 6. 构建三个数据专用先验数据集 ----------
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

prior_H  <- build_dataset_prior(prior_df, expr_H$Gene)
prior_TI <- build_dataset_prior(prior_df, expr_TI$Gene)
prior_TM <- build_dataset_prior(prior_df, expr_TM$Gene)

# 保存专用先验数据集（CSV）
write.csv(prior_H,
          file.path(data_dir, "GSE117221_H_prior_dataset.csv"),
          row.names = FALSE)

write.csv(prior_TI,
          file.path(data_dir, "GSE117221_TI_prior_dataset.csv"),
          row.names = FALSE)

write.csv(prior_TM,
          file.path(data_dir, "GSE117221_TM_prior_dataset.csv"),
          row.names = FALSE)

# ---------- 7. 将专用先验数据集映射为基因ID对 ----------
# 改用命名向量映射，避免 merge() 排序导致顺序问题
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

prior_H_ids  <- convert_prior_to_ids(prior_H,  gene_id_H)
prior_TI_ids <- convert_prior_to_ids(prior_TI, gene_id_TI)
prior_TM_ids <- convert_prior_to_ids(prior_TM, gene_id_TM)

# 保存两列ID文件（TSV，无表头）
write.table(prior_H_ids,
            file.path(data_dir, "GSE117221_H_prior_ids.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

write.table(prior_TI_ids,
            file.path(data_dir, "GSE117221_TI_prior_ids.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

write.table(prior_TM_ids,
            file.path(data_dir, "GSE117221_TM_prior_ids.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# ---------- 8. 生成带Type的ID版先验文件 ----------
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

prior_H_ids_type  <- convert_prior_to_ids_with_type(prior_H,  gene_id_H)
prior_TI_ids_type <- convert_prior_to_ids_with_type(prior_TI, gene_id_TI)
prior_TM_ids_type <- convert_prior_to_ids_with_type(prior_TM, gene_id_TM)

# 保存带Type的ID版先验文件（TSV，有表头）
write.table(prior_H_ids_type,
            file.path(data_dir, "GSE117221_H_prior_ids_with_type.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

write.table(prior_TI_ids_type,
            file.path(data_dir, "GSE117221_TI_prior_ids_with_type.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

write.table(prior_TM_ids_type,
            file.path(data_dir, "GSE117221_TM_prior_ids_with_type.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

# ---------- 9. 可选：反查验证文件 ----------
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

check_H  <- build_check_table(prior_H,  gene_id_H)
check_TI <- build_check_table(prior_TI, gene_id_TI)
check_TM <- build_check_table(prior_TM, gene_id_TM)

# 保存核对文件（TSV）
write.table(check_H,
            file.path(data_dir, "GSE117221_H_prior_check_table.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

write.table(check_TI,
            file.path(data_dir, "GSE117221_TI_prior_check_table.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

write.table(check_TM,
            file.path(data_dir, "GSE117221_TM_prior_check_table.tsv"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

# ---------- 10. 输出结果统计 ----------
cat("处理完成！\n\n")

cat("三个基因ID文件（TXT）：\n")
cat("H  :", nrow(gene_id_H),  "个基因\n")
cat("TI :", nrow(gene_id_TI), "个基因\n")
cat("TM :", nrow(gene_id_TM), "个基因\n\n")

cat("三个专用先验数据集（CSV）：\n")
cat("H  :", nrow(prior_H),  "条\n")
cat("TI :", nrow(prior_TI), "条\n")
cat("TM :", nrow(prior_TM), "条\n\n")

cat("三个ID版先验文件（TSV）：\n")
cat("H  :", nrow(prior_H_ids),  "条\n")
cat("TI :", nrow(prior_TI_ids), "条\n")
cat("TM :", nrow(prior_TM_ids), "条\n\n")

cat("所有结果已保存到：\n", data_dir, "\n")







