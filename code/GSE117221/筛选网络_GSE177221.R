# ================================
# GSE117221 - 9个GRN统一处理
# 1. 3DCEMA：按每组特有TF列表筛选
# 2. DeepFGRN：ID转基因名 + 合成EdgeWeight + 按每组TF列表筛选
# 3. DeepSEM：按总人类TF列表筛选
# 4. 对9个筛选后的网络保留前10%高置信边
# ================================

# ---------- 1. 加载包 ----------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)

# ---------- 2. 路径设置 ----------
base_dir <- "E:/SCD/其他数据/GSE117221"

# 3DCEMA
file_3dcema_H  <- file.path(base_dir, "3dcema", "3dcema_H.csv")
file_3dcema_TI <- file.path(base_dir, "3dcema", "3DCEMA_TI.csv")
file_3dcema_TM <- file.path(base_dir, "3dcema", "3dcema_TM.csv")

# DeepFGRN
file_deepfgrn_H  <- file.path(base_dir, "deepfgrn", "deepfgrn_H.tsv")
file_deepfgrn_TI <- file.path(base_dir, "deepfgrn", "deepfgrn_TI.tsv")
file_deepfgrn_TM <- file.path(base_dir, "deepfgrn", "deepfgrn_TM.tsv")

# DeepSEM
file_deepsem_H  <- file.path(base_dir, "deepsem", "deepsem_H.tsv")
file_deepsem_TI <- file.path(base_dir, "deepsem", "deepsem_TI .tsv")   # 如果真实文件名没有空格，请手动改成 deepsem_TI.tsv
file_deepsem_TM <- file.path(base_dir, "deepsem", "deepsem_TM.tsv")

# 每组专用TF列表
tf_H_file  <- file.path(base_dir, "GSE117221_H_top10000_TF_list.txt")
tf_TI_file <- file.path(base_dir, "GSE117221_TI_top10000_TF_list.txt")
tf_TM_file <- file.path(base_dir, "GSE117221_TM_top10000_TF_list.txt")

# 每组基因ID文件
gene_id_H_file  <- file.path(base_dir, "GSE117221_H_gene_ids.txt")
gene_id_TI_file <- file.path(base_dir, "GSE117221_TI_gene_ids.txt")
gene_id_TM_file <- file.path(base_dir, "GSE117221_TM_gene_ids.txt")

# 总人类TF列表
human_tf_file <- "E:/SCD/数据/构建网络的数据/human_tf_list.txt"

# 输出目录：统一改到 grn
out_dir <- file.path(base_dir, "grn")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- 3. 读取辅助文件 ----------
read_tf_list <- function(file_path, col_name = "TF") {
  df <- fread(file_path, data.table = FALSE)
  colnames(df)[1] <- col_name
  tf_vec <- trimws(as.character(df[[col_name]]))
  tf_vec <- unique(tf_vec[!is.na(tf_vec) & tf_vec != ""])
  return(tf_vec)
}

read_gene_id <- function(file_path) {
  df <- fread(file_path, data.table = FALSE)
  colnames(df)[1:2] <- c("name", "ids")
  df$name <- trimws(as.character(df$name))
  df$ids  <- as.integer(df$ids)
  df <- df[!is.na(df$name) & df$name != "", , drop = FALSE]
  df <- df[!duplicated(df$ids), , drop = FALSE]
  return(df)
}

tf_H  <- read_tf_list(tf_H_file)
tf_TI <- read_tf_list(tf_TI_file)
tf_TM <- read_tf_list(tf_TM_file)

human_tf <- read_tf_list(human_tf_file)

gene_id_H  <- read_gene_id(gene_id_H_file)
gene_id_TI <- read_gene_id(gene_id_TI_file)
gene_id_TM <- read_gene_id(gene_id_TM_file)

# ---------- 4. 通用函数 ----------
save_tsv <- function(df, file_path) {
  fwrite(df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# 取前10%高置信边
# mode = "raw_desc"  按 EdgeWeight 从大到小
# mode = "abs_desc"  按 abs(EdgeWeight) 从大到小
filter_top10 <- function(df, mode = "raw_desc") {
  if (nrow(df) == 0) return(df)
  
  if (mode == "raw_desc") {
    score <- df$EdgeWeight
  } else if (mode == "abs_desc") {
    score <- abs(df$EdgeWeight)
  } else {
    stop("mode 只能是 'raw_desc' 或 'abs_desc'")
  }
  
  top_n <- ceiling(nrow(df) * 0.10)
  top_n <- max(top_n, 1)
  
  ord <- order(score, decreasing = TRUE)
  df2 <- df[ord, , drop = FALSE]
  df2 <- df2[1:top_n, , drop = FALSE]
  
  rownames(df2) <- NULL
  return(df2)
}

# ---------- 5. 处理 3DCEMA ----------
process_3dcema <- function(file_path, tf_list) {
  df <- fread(file_path, data.table = FALSE)
  
  # 统一列名
  colnames(df)[1:3] <- c("TF", "Target", "EdgeWeight")
  
  df$TF <- trimws(as.character(df$TF))
  df$Target <- trimws(as.character(df$Target))
  df$EdgeWeight <- as.numeric(df$EdgeWeight)
  
  # 只保留本组特有TF
  df <- df[df$TF %in% tf_list, c("TF", "Target", "EdgeWeight"), drop = FALSE]
  
  # 去掉缺失
  df <- df[!is.na(df$TF) & !is.na(df$Target) & !is.na(df$EdgeWeight), , drop = FALSE]
  
  rownames(df) <- NULL
  return(df)
}

net_3dcema_H  <- process_3dcema(file_3dcema_H,  tf_H)
net_3dcema_TI <- process_3dcema(file_3dcema_TI, tf_TI)
net_3dcema_TM <- process_3dcema(file_3dcema_TM, tf_TM)

save_tsv(net_3dcema_H,  file.path(out_dir, "3dcema_H_filtered.tsv"))
save_tsv(net_3dcema_TI, file.path(out_dir, "3dcema_TI_filtered.tsv"))
save_tsv(net_3dcema_TM, file.path(out_dir, "3dcema_TM_filtered.tsv"))

# ---------- 6. 处理 DeepFGRN ----------
# 新规则：
# pred_label == 0 -> EdgeWeight = 1 - prob_no
# pred_label == 1 -> EdgeWeight = prob_act
# pred_label == 2 -> EdgeWeight = -prob_rep
process_deepfgrn <- function(file_path, gene_id_df, tf_list) {
  df <- fread(file_path, data.table = FALSE)
  
  colnames(df)[1:6] <- c("TF_id", "Target_id", "prob_no", "prob_act", "prob_rep", "pred_label")
  
  df$TF_id <- as.integer(df$TF_id)
  df$Target_id <- as.integer(df$Target_id)
  df$prob_no <- as.numeric(df$prob_no)
  df$prob_act <- as.numeric(df$prob_act)
  df$prob_rep <- as.numeric(df$prob_rep)
  df$pred_label <- as.integer(df$pred_label)
  
  # gene_id 映射
  id_to_gene <- setNames(gene_id_df$name, gene_id_df$ids)
  
  df$TF <- unname(id_to_gene[as.character(df$TF_id)])
  df$Target <- unname(id_to_gene[as.character(df$Target_id)])
  
  # 合成 EdgeWeight
  df$EdgeWeight <- ifelse(
    df$pred_label == 0, 1 - df$prob_no,
    ifelse(df$pred_label == 1, df$prob_act,
           ifelse(df$pred_label == 2, -df$prob_rep, NA))
  )
  
  # 只保留本组TF列表中的TF
  df <- df[df$TF %in% tf_list, c("TF", "Target", "EdgeWeight"), drop = FALSE]
  
  # 去掉缺失
  df <- df[!is.na(df$TF) & !is.na(df$Target) & !is.na(df$EdgeWeight), , drop = FALSE]
  
  rownames(df) <- NULL
  return(df)
}

net_deepfgrn_H  <- process_deepfgrn(file_deepfgrn_H,  gene_id_H,  tf_H)
net_deepfgrn_TI <- process_deepfgrn(file_deepfgrn_TI, gene_id_TI, tf_TI)
net_deepfgrn_TM <- process_deepfgrn(file_deepfgrn_TM, gene_id_TM, tf_TM)

save_tsv(net_deepfgrn_H,  file.path(out_dir, "deepfgrn_H_filtered.tsv"))
save_tsv(net_deepfgrn_TI, file.path(out_dir, "deepfgrn_TI_filtered.tsv"))
save_tsv(net_deepfgrn_TM, file.path(out_dir, "deepfgrn_TM_filtered.tsv"))

# ---------- 7. 处理 DeepSEM ----------
process_deepsem <- function(file_path, human_tf) {
  df <- fread(file_path, data.table = FALSE)
  
  colnames(df)[1:3] <- c("TF", "Target", "EdgeWeight")
  
  df$TF <- trimws(as.character(df$TF))
  df$Target <- trimws(as.character(df$Target))
  df$EdgeWeight <- as.numeric(df$EdgeWeight)
  
  # 用总人类TF列表筛选
  df <- df[df$TF %in% human_tf, c("TF", "Target", "EdgeWeight"), drop = FALSE]
  
  # 去掉缺失
  df <- df[!is.na(df$TF) & !is.na(df$Target) & !is.na(df$EdgeWeight), , drop = FALSE]
  
  rownames(df) <- NULL
  return(df)
}

net_deepsem_H  <- process_deepsem(file_deepsem_H,  human_tf)
net_deepsem_TI <- process_deepsem(file_deepsem_TI, human_tf)
net_deepsem_TM <- process_deepsem(file_deepsem_TM, human_tf)

save_tsv(net_deepsem_H,  file.path(out_dir, "deepsem_H_filtered.tsv"))
save_tsv(net_deepsem_TI, file.path(out_dir, "deepsem_TI_filtered.tsv"))
save_tsv(net_deepsem_TM, file.path(out_dir, "deepsem_TM_filtered.tsv"))

# ---------- 8. 9个网络按前10%高置信边筛选 ----------
# 3DCEMA：按原始 EdgeWeight 从大到小
net_3dcema_H_top10  <- filter_top10(net_3dcema_H,  mode = "raw_desc")
net_3dcema_TI_top10 <- filter_top10(net_3dcema_TI, mode = "raw_desc")
net_3dcema_TM_top10 <- filter_top10(net_3dcema_TM, mode = "raw_desc")

# DeepFGRN：按 abs(EdgeWeight) 从大到小
net_deepfgrn_H_top10  <- filter_top10(net_deepfgrn_H,  mode = "abs_desc")
net_deepfgrn_TI_top10 <- filter_top10(net_deepfgrn_TI, mode = "abs_desc")
net_deepfgrn_TM_top10 <- filter_top10(net_deepfgrn_TM, mode = "abs_desc")

# DeepSEM：按原始 EdgeWeight 从大到小
net_deepsem_H_top10  <- filter_top10(net_deepsem_H,  mode = "raw_desc")
net_deepsem_TI_top10 <- filter_top10(net_deepsem_TI, mode = "raw_desc")
net_deepsem_TM_top10 <- filter_top10(net_deepsem_TM, mode = "raw_desc")

# 保存 top10% 高置信网络
save_tsv(net_3dcema_H_top10,  file.path(out_dir, "3dcema_H_top10pct.tsv"))
save_tsv(net_3dcema_TI_top10, file.path(out_dir, "3dcema_TI_top10pct.tsv"))
save_tsv(net_3dcema_TM_top10, file.path(out_dir, "3dcema_TM_top10pct.tsv"))

save_tsv(net_deepfgrn_H_top10,  file.path(out_dir, "deepfgrn_H_top10pct.tsv"))
save_tsv(net_deepfgrn_TI_top10, file.path(out_dir, "deepfgrn_TI_top10pct.tsv"))
save_tsv(net_deepfgrn_TM_top10, file.path(out_dir, "deepfgrn_TM_top10pct.tsv"))

save_tsv(net_deepsem_H_top10,  file.path(out_dir, "deepsem_H_top10pct.tsv"))
save_tsv(net_deepsem_TI_top10, file.path(out_dir, "deepsem_TI_top10pct.tsv"))
save_tsv(net_deepsem_TM_top10, file.path(out_dir, "deepsem_TM_top10pct.tsv"))

# ---------- 9. 输出统计 ----------
cat("处理完成！\n\n")

cat("==== 3DCEMA 筛选后 ====\n")
cat("H  :", nrow(net_3dcema_H),  "条边\n")
cat("TI :", nrow(net_3dcema_TI), "条边\n")
cat("TM :", nrow(net_3dcema_TM), "条边\n\n")

cat("==== DeepFGRN 筛选后 ====\n")
cat("H  :", nrow(net_deepfgrn_H),  "条边\n")
cat("TI :", nrow(net_deepfgrn_TI), "条边\n")
cat("TM :", nrow(net_deepfgrn_TM), "条边\n\n")

cat("==== DeepSEM 筛选后 ====\n")
cat("H  :", nrow(net_deepsem_H),  "条边\n")
cat("TI :", nrow(net_deepsem_TI), "条边\n")
cat("TM :", nrow(net_deepsem_TM), "条边\n\n")

cat("==== Top10% 高置信边 ====\n")
cat("3DCEMA H  :", nrow(net_3dcema_H_top10), "\n")
cat("3DCEMA TI :", nrow(net_3dcema_TI_top10), "\n")
cat("3DCEMA TM :", nrow(net_3dcema_TM_top10), "\n")
cat("DeepFGRN H  :", nrow(net_deepfgrn_H_top10), "\n")
cat("DeepFGRN TI :", nrow(net_deepfgrn_TI_top10), "\n")
cat("DeepFGRN TM :", nrow(net_deepfgrn_TM_top10), "\n")
cat("DeepSEM H  :", nrow(net_deepsem_H_top10), "\n")
cat("DeepSEM TI :", nrow(net_deepsem_TI_top10), "\n")
cat("DeepSEM TM :", nrow(net_deepsem_TM_top10), "\n\n")

cat("所有输出文件保存在：\n", out_dir, "\n")