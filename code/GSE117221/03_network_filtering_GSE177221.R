# ==============================================================================
# Script: Unified Processing of 6 GRNs from GSE117221 (H and TI groups only)
# 
# Description:
# This script processes gene regulatory networks (GRNs) inferred by three 
# different methods (3DCEMA, DeepFGRN, DeepSEM) for two sample groups 
# (H and TI) from the GSE117221 dataset.
#
# Processing steps for each method:
#   1. 3DCEMA: Filter edges by group-specific TF lists
#   2. DeepFGRN: Convert gene IDs to gene symbols, synthesize EdgeWeight values,
#                then filter by group-specific TF lists
#   3. DeepSEM: Filter edges by the general human TF list
#   4. For all 6 filtered networks, retain the top 10% highest-confidence edges
#
# Input data:
#   - Network files:
#     * 3DCEMA: 3dcema_H.csv, 3DCEMA_TI.csv
#     * DeepFGRN: deepfgrn_H.tsv, deepfgrn_TI.tsv
#     * DeepSEM: deepsem_H.tsv, deepsem_TI.tsv
#   - Group-specific TF lists: GSE117221_H_top10000_TF_list.txt, 
#                              GSE117221_TI_top10000_TF_list.txt
#   - Gene ID mapping files: GSE117221_H_gene_ids.txt, GSE117221_TI_gene_ids.txt
#   - General human TF list: human_tf_list.txt
#
# Output data (saved in the 'grn' subdirectory):
#   - Filtered networks (all edges after TF filtering):
#     * 3dcema_H_filtered.tsv
#     * 3dcema_TI_filtered.tsv
#     * deepfgrn_H_filtered.tsv
#     * deepfgrn_TI_filtered.tsv
#     * deepsem_H_filtered.tsv
#     * deepsem_TI_filtered.tsv
#   - Top 10% high-confidence networks:
#     * 3dcema_H_top10pct.tsv
#     * 3dcema_TI_top10pct.tsv
#     * deepfgrn_H_top10pct.tsv
#     * deepfgrn_TI_top10pct.tsv
#     * deepsem_H_top10pct.tsv
#     * deepsem_TI_top10pct.tsv
# ==============================================================================

# ---------- 1. Load required packages ----------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)

# ---------- 2. Set file paths ----------
base_dir <- "E:/SCD/其他数据/GSE117221"

# 3DCEMA files
file_3dcema_H  <- file.path(base_dir, "3dcema", "3dcema_H.csv")
file_3dcema_TI <- file.path(base_dir, "3dcema", "3DCEMA_TI.csv")

# DeepFGRN files
file_deepfgrn_H  <- file.path(base_dir, "deepfgrn", "deepfgrn_H.tsv")
file_deepfgrn_TI <- file.path(base_dir, "deepfgrn", "deepfgrn_TI.tsv")

# DeepSEM files
file_deepsem_H  <- file.path(base_dir, "deepsem", "deepsem_H.tsv")
file_deepsem_TI <- file.path(base_dir, "deepsem", "deepsem_TI .tsv")   # Please remove the space if the actual filename has no space

# Group-specific TF list files
tf_H_file  <- file.path(base_dir, "GSE117221_H_top10000_TF_list.txt")
tf_TI_file <- file.path(base_dir, "GSE117221_TI_top10000_TF_list.txt")

# Gene ID mapping files
gene_id_H_file  <- file.path(base_dir, "GSE117221_H_gene_ids.txt")
gene_id_TI_file <- file.path(base_dir, "GSE117221_TI_gene_ids.txt")

# General human TF list file
human_tf_file <- "E:/SCD/数据/构建网络的数据/human_tf_list.txt"

# Output directory
out_dir <- file.path(base_dir, "grn")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- 3. Read auxiliary files ----------
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

human_tf <- read_tf_list(human_tf_file)

gene_id_H  <- read_gene_id(gene_id_H_file)
gene_id_TI <- read_gene_id(gene_id_TI_file)

# ---------- 4. Utility functions ----------
save_tsv <- function(df, file_path) {
  fwrite(df, file = file_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Filter to keep top 10% edges based on confidence scores
# mode = "raw_desc": sort by EdgeWeight in descending order
# mode = "abs_desc": sort by abs(EdgeWeight) in descending order
filter_top10 <- function(df, mode = "raw_desc") {
  if (nrow(df) == 0) return(df)
  
  if (mode == "raw_desc") {
    score <- df$EdgeWeight
  } else if (mode == "abs_desc") {
    score <- abs(df$EdgeWeight)
  } else {
    stop("mode must be either 'raw_desc' or 'abs_desc'")
  }
  
  top_n <- ceiling(nrow(df) * 0.10)
  top_n <- max(top_n, 1)
  
  ord <- order(score, decreasing = TRUE)
  df2 <- df[ord, , drop = FALSE]
  df2 <- df2[1:top_n, , drop = FALSE]
  
  rownames(df2) <- NULL
  return(df2)
}

# ---------- 5. Process 3DCEMA ----------
process_3dcema <- function(file_path, tf_list) {
  df <- fread(file_path, data.table = FALSE)
  
  # Standardize column names
  colnames(df)[1:3] <- c("TF", "Target", "EdgeWeight")
  
  df$TF <- trimws(as.character(df$TF))
  df$Target <- trimws(as.character(df$Target))
  df$EdgeWeight <- as.numeric(df$EdgeWeight)
  
  # Keep only edges with group-specific TFs
  df <- df[df$TF %in% tf_list, c("TF", "Target", "EdgeWeight"), drop = FALSE]
  
  # Remove rows with missing values
  df <- df[!is.na(df$TF) & !is.na(df$Target) & !is.na(df$EdgeWeight), , drop = FALSE]
  
  rownames(df) <- NULL
  return(df)
}

net_3dcema_H  <- process_3dcema(file_3dcema_H,  tf_H)
net_3dcema_TI <- process_3dcema(file_3dcema_TI, tf_TI)

save_tsv(net_3dcema_H,  file.path(out_dir, "3dcema_H_filtered.tsv"))
save_tsv(net_3dcema_TI, file.path(out_dir, "3dcema_TI_filtered.tsv"))

# ---------- 6. Process DeepFGRN ----------
# Rules for synthesizing EdgeWeight:
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
  
  # Map gene IDs to gene symbols
  id_to_gene <- setNames(gene_id_df$name, gene_id_df$ids)
  
  df$TF <- unname(id_to_gene[as.character(df$TF_id)])
  df$Target <- unname(id_to_gene[as.character(df$Target_id)])
  
  # Synthesize EdgeWeight
  df$EdgeWeight <- ifelse(
    df$pred_label == 0, 1 - df$prob_no,
    ifelse(df$pred_label == 1, df$prob_act,
           ifelse(df$pred_label == 2, -df$prob_rep, NA))
  )
  
  # Keep only edges with group-specific TFs
  df <- df[df$TF %in% tf_list, c("TF", "Target", "EdgeWeight"), drop = FALSE]
  
  # Remove rows with missing values
  df <- df[!is.na(df$TF) & !is.na(df$Target) & !is.na(df$EdgeWeight), , drop = FALSE]
  
  rownames(df) <- NULL
  return(df)
}

net_deepfgrn_H  <- process_deepfgrn(file_deepfgrn_H,  gene_id_H,  tf_H)
net_deepfgrn_TI <- process_deepfgrn(file_deepfgrn_TI, gene_id_TI, tf_TI)

save_tsv(net_deepfgrn_H,  file.path(out_dir, "deepfgrn_H_filtered.tsv"))
save_tsv(net_deepfgrn_TI, file.path(out_dir, "deepfgrn_TI_filtered.tsv"))

# ---------- 7. Process DeepSEM ----------
process_deepsem <- function(file_path, human_tf) {
  df <- fread(file_path, data.table = FALSE)
  
  colnames(df)[1:3] <- c("TF", "Target", "EdgeWeight")
  
  df$TF <- trimws(as.character(df$TF))
  df$Target <- trimws(as.character(df$Target))
  df$EdgeWeight <- as.numeric(df$EdgeWeight)
  
  # Keep only edges with TFs from the general human TF list
  df <- df[df$TF %in% human_tf, c("TF", "Target", "EdgeWeight"), drop = FALSE]
  
  # Remove rows with missing values
  df <- df[!is.na(df$TF) & !is.na(df$Target) & !is.na(df$EdgeWeight), , drop = FALSE]
  
  rownames(df) <- NULL
  return(df)
}

net_deepsem_H  <- process_deepsem(file_deepsem_H,  human_tf)
net_deepsem_TI <- process_deepsem(file_deepsem_TI, human_tf)

save_tsv(net_deepsem_H,  file.path(out_dir, "deepsem_H_filtered.tsv"))
save_tsv(net_deepsem_TI, file.path(out_dir, "deepsem_TI_filtered.tsv"))

# ---------- 8. Filter top 10% high-confidence edges for all 6 networks ----------
# 3DCEMA: sort by raw EdgeWeight in descending order
net_3dcema_H_top10  <- filter_top10(net_3dcema_H,  mode = "raw_desc")
net_3dcema_TI_top10 <- filter_top10(net_3dcema_TI, mode = "raw_desc")

# DeepFGRN: sort by abs(EdgeWeight) in descending order
net_deepfgrn_H_top10  <- filter_top10(net_deepfgrn_H,  mode = "abs_desc")
net_deepfgrn_TI_top10 <- filter_top10(net_deepfgrn_TI, mode = "abs_desc")

# DeepSEM: sort by raw EdgeWeight in descending order
net_deepsem_H_top10  <- filter_top10(net_deepsem_H,  mode = "raw_desc")
net_deepsem_TI_top10 <- filter_top10(net_deepsem_TI, mode = "raw_desc")

# Save top 10% high-confidence networks
save_tsv(net_3dcema_H_top10,  file.path(out_dir, "3dcema_H_top10pct.tsv"))
save_tsv(net_3dcema_TI_top10, file.path(out_dir, "3dcema_TI_top10pct.tsv"))

save_tsv(net_deepfgrn_H_top10,  file.path(out_dir, "deepfgrn_H_top10pct.tsv"))
save_tsv(net_deepfgrn_TI_top10, file.path(out_dir, "deepfgrn_TI_top10pct.tsv"))

save_tsv(net_deepsem_H_top10,  file.path(out_dir, "deepsem_H_top10pct.tsv"))
save_tsv(net_deepsem_TI_top10, file.path(out_dir, "deepsem_TI_top10pct.tsv"))

# ---------- 9. Output statistics ----------
cat("Processing complete!\n\n")

cat("==== 3DCEMA after filtering ====\n")
cat("H  :", nrow(net_3dcema_H),  "edges\n")
cat("TI :", nrow(net_3dcema_TI), "edges\n\n")

cat("==== DeepFGRN after filtering ====\n")
cat("H  :", nrow(net_deepfgrn_H),  "edges\n")
cat("TI :", nrow(net_deepfgrn_TI), "edges\n\n")

cat("==== DeepSEM after filtering ====\n")
cat("H  :", nrow(net_deepsem_H),  "edges\n")
cat("TI :", nrow(net_deepsem_TI), "edges\n\n")

cat("==== Top 10% high-confidence edges ====\n")
cat("3DCEMA H  :", nrow(net_3dcema_H_top10), "\n")
cat("3DCEMA TI :", nrow(net_3dcema_TI_top10), "\n")
cat("DeepFGRN H  :", nrow(net_deepfgrn_H_top10), "\n")
cat("DeepFGRN TI :", nrow(net_deepfgrn_TI_top10), "\n")
cat("DeepSEM H  :", nrow(net_deepsem_H_top10), "\n")
cat("DeepSEM TI :", nrow(net_deepsem_TI_top10), "\n\n")

cat("All output files saved in:\n", out_dir, "\n")