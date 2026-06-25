# ==============================================================================
# Script: GSE117221 Data Processing Pipeline (H and TI groups only)
# 
# Description:
# This script processes expression count data from GSE117221 for two sample 
# groups (H and TI) through a comprehensive pipeline:
#
# Part 1: Merge individual count files into group-level expression matrices
# Part 2: Filter top 5000 most variable genes and transpose matrices
# Part 3: Generate group-specific TF lists from top 5000 genes
# Part 4: Build gene ID mapping files and convert prior regulatory pairs
#
# Input data:
#   - Individual count files: *_counts.txt files for each sample
#   - General human TF list: human_tf_list.txt
#   - Prior regulatory network: label_c1.csv / new_GRN_Lung_GEN_counts_genename_c1.csv
#
# Output data (saved in the same directory):
#   Part 1 - Merged count matrices:
#     * GSE117221_H_counts_merged.txt
#     * GSE117221_TI_counts_merged.txt
#   Part 2 - Transposed expression matrices:
#     * GSE117221_H_transposed.csv, GSE117221_TI_transposed.csv
#     * GSE117221_H_top5000_transposed.csv, GSE117221_TI_top5000_transposed.csv
#     * GSE117221_H_top5000_notransposed.csv, GSE117221_TI_top5000_notransposed.csv
#   Part 3 - Gene lists and TF lists:
#     * GSE117221_H_top5000_genes.csv, GSE117221_TI_top5000_genes.csv
#     * GSE117221_H_top5000_TF_list.txt, GSE117221_TI_top5000_TF_list.txt
#     * GSE117221_H_gene_variance.csv, GSE117221_TI_gene_variance.csv
#   Part 4 - Gene ID mappings and prior pairs:
#     * GSE117221_H_gene_ids.txt, GSE117221_TI_gene_ids.txt
#     * GSE117221_H_prior_dataset.csv, GSE117221_TI_prior_dataset.csv
#     * GSE117221_H_prior_ids.tsv, GSE117221_TI_prior_ids.tsv
#     * GSE117221_H_prior_ids_with_type.tsv, GSE117221_TI_prior_ids_with_type.tsv
#     * GSE117221_H_prior_check_table.tsv, GSE117221_TI_prior_check_table.tsv
# ==============================================================================

# ==============================================================================
# PART 1: Merge individual count files into group-level matrices
# ==============================================================================

# ---------- 1. Set working directory ----------
data_dir <- "E:/SCD/其他数据/GSE117221"

# ---------- 2. Get all count files ----------
files <- list.files(
  path = data_dir,
  pattern = "_counts\\.txt$",
  full.names = TRUE
)

cat("Total files detected:", length(files), "\n")

# ---------- 3. Define function to read individual count files ----------
# Assume each file has two columns: gene and count, no header
read_count_file <- function(file_path) {
  df <- read.table(
    file_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  
  # If not tab-separated, try reading with default separator
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
    stop(paste("File format issue, cannot read two columns:", file_path))
  }
  
  df <- df[, 1:2]
  colnames(df) <- c("Gene", "Count")
  
  # Extract sample name: remove path and suffix
  sample_name <- basename(file_path)
  sample_name <- sub("_counts\\.txt$", "", sample_name)
  
  colnames(df)[2] <- sample_name
  return(df)
}

# ---------- 4. Read all files ----------
count_list <- lapply(files, read_count_file)

# ---------- 5. Extract sample names ----------
sample_names <- sapply(count_list, function(x) colnames(x)[2])

# ---------- 6. Group samples by filename pattern ----------
# Note: Check TM first, then TI/TITD, then H to avoid mismatches
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

group_info <- data.frame(
  sample = sample_names,
  group = group_labels,
  stringsAsFactors = FALSE
)

cat("Sample grouping:\n")
print(group_info)

if (any(is.na(group_labels))) {
  cat("\nThe following samples were not successfully grouped, please check filenames:\n")
  print(group_info[is.na(group_labels), ])
}

# ---------- 7. Define merge function ----------
merge_counts <- function(df_list) {
  merged_df <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), df_list)
  
  # Fill NAs with 0
  merged_df[is.na(merged_df)] <- 0
  
  # Ensure count columns are numeric
  for (i in 2:ncol(merged_df)) {
    merged_df[[i]] <- as.numeric(merged_df[[i]])
  }
  
  return(merged_df)
}

# ---------- 8. Extract H and TI groups only ----------
H_list  <- count_list[group_labels == "H"]
TI_list <- count_list[group_labels == "TI"]

# ---------- 9. Merge ----------
H_merged  <- merge_counts(H_list)
TI_merged <- merge_counts(TI_list)

# ---------- 10. Sort by gene name ----------
H_merged  <- H_merged[order(H_merged$Gene), ]
TI_merged <- TI_merged[order(TI_merged$Gene), ]

# ---------- 11. Save merged matrices ----------
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

cat("\nMerge complete! Output files:\n")
cat(file.path(data_dir, "GSE117221_H_counts_merged.txt"), "\n")
cat(file.path(data_dir, "GSE117221_TI_counts_merged.txt"), "\n")

# ---------- 12. Check dimensions ----------
cat("\nMatrix dimensions:\n")
cat("H group :", dim(H_merged)[1], "genes x", dim(H_merged)[2] - 1, "samples\n")
cat("TI group:", dim(TI_merged)[1], "genes x", dim(TI_merged)[2] - 1, "samples\n")


# ==============================================================================
# PART 2: Filter top 5000 genes and transpose matrices
# ==============================================================================

# ---------- 1. Load required packages ----------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("edgeR", quietly = TRUE)) install.packages("edgeR")
if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")

library(data.table)
library(edgeR)
library(matrixStats)

# ---------- 2. Set file paths ----------
file_H  <- file.path(data_dir, "GSE117221_H_counts_merged.txt")
file_TI <- file.path(data_dir, "GSE117221_TI_counts_merged.txt")

# Prior regulatory pairs file
prior_file <- "E:/SCD/数据/构建网络的数据/label_c1.csv"

# ---------- 3. Read expression matrices ----------
read_expr_matrix <- function(file_path) {
  df <- fread(file_path, data.table = FALSE)
  
  # Rename first column to Gene
  colnames(df)[1] <- "Gene"
  
  # Remove duplicates: keep first occurrence
  df <- df[!duplicated(df$Gene), ]
  
  # Set gene names as rownames
  rownames(df) <- df$Gene
  df$Gene <- NULL
  
  # Convert to numeric matrix
  mat <- as.matrix(df)
  mode(mat) <- "numeric"
  
  return(mat)   # rows=genes, columns=samples
}

expr_H  <- read_expr_matrix(file_H)
expr_TI <- read_expr_matrix(file_TI)

# ---------- 4. Define gene selection function ----------
# Strategy:
# (1) Filter low expression genes by CPM
# (2) Select top 5000 genes by log2(CPM+1) variance
select_top_genes <- function(count_mat, top_n = 5000, cpm_cutoff = 1, min_prop = 0.2) {
  # count_mat: rows=genes, columns=samples
  
  # Calculate minimum number of samples with expression
  min_samples <- ceiling(ncol(count_mat) * min_prop)
  
  # CPM
  cpm_mat <- edgeR::cpm(count_mat)
  
  # Low expression filtering
  keep <- rowSums(cpm_mat > cpm_cutoff) >= min_samples
  filtered_mat <- count_mat[keep, , drop = FALSE]
  
  # Calculate log2(CPM+1)
  log_cpm <- log2(edgeR::cpm(filtered_mat) + 1)
  
  # Calculate variance
  gene_var <- matrixStats::rowVars(log_cpm)
  
  # If fewer than top_n genes remain, keep all
  top_n_use <- min(top_n, length(gene_var))
  
  # Select top genes by variance
  top_idx <- order(gene_var, decreasing = TRUE)[1:top_n_use]
  top_genes <- rownames(log_cpm)[top_idx]
  
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

# ---------- 5. Transpose matrices ----------
# After transposition: rows=samples, columns=genes
transpose_expr <- function(mat) {
  tmat <- as.data.frame(t(mat), check.names = FALSE)
  tmat <- cbind(Sample = rownames(tmat), tmat)
  rownames(tmat) <- NULL
  return(tmat)
}

# Full matrices transposed
expr_H_t  <- transpose_expr(expr_H)
expr_TI_t <- transpose_expr(expr_TI)

# Top 5000 matrices transposed
top_H_t  <- transpose_expr(res_H$top_count_mat)
top_TI_t <- transpose_expr(res_TI$top_count_mat)

# ---------- 6. Restore non-transposed top5000 matrices from transposed versions ----------
read_transposed_expr <- function(file_path) {
  df <- fread(file_path, data.table = FALSE, check.names = FALSE)
  colnames(df)[1] <- "Sample"
  return(df)
}

# Save top5000 transposed first
write.csv(top_H_t,  file.path(data_dir, "GSE117221_H_top5000_transposed.csv"),  row.names = FALSE)
write.csv(top_TI_t, file.path(data_dir, "GSE117221_TI_top5000_transposed.csv"), row.names = FALSE)

# Then restore non-transposed version
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

expr_H_notrans  <- restore_non_transposed(top_H_t)
expr_TI_notrans <- restore_non_transposed(top_TI_t)

# Save restored non-transposed matrices
write.csv(expr_H_notrans,
          file.path(data_dir, "GSE117221_H_top5000_notransposed.csv"),
          row.names = FALSE)

write.csv(expr_TI_notrans,
          file.path(data_dir, "GSE117221_TI_top5000_notransposed.csv"),
          row.names = FALSE)

# ---------- 7. Save other outputs ----------
# Full transposed matrices
write.csv(expr_H_t,  file.path(data_dir, "GSE117221_H_transposed.csv"),  row.names = FALSE)
write.csv(expr_TI_t, file.path(data_dir, "GSE117221_TI_transposed.csv"), row.names = FALSE)

# Top 5000 gene lists
write.csv(data.frame(Gene = res_H$top_genes),  file.path(data_dir, "GSE117221_H_top5000_genes.csv"),  row.names = FALSE)
write.csv(data.frame(Gene = res_TI$top_genes), file.path(data_dir, "GSE117221_TI_top5000_genes.csv"), row.names = FALSE)

# Gene variance tables
write.csv(res_H$gene_variance,  file.path(data_dir, "GSE117221_H_gene_variance.csv"),  row.names = FALSE)
write.csv(res_TI$gene_variance, file.path(data_dir, "GSE117221_TI_gene_variance.csv"), row.names = FALSE)

# ---------- 8. Build group-specific prior pairs ----------
prior_df <- fread(prior_file, data.table = FALSE)
colnames(prior_df)[1:2] <- c("Gene1", "Gene2")
prior_df <- prior_df[!is.na(prior_df$Gene1) & !is.na(prior_df$Gene2), ]
prior_df <- unique(prior_df)

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

write.csv(prior_H,  file.path(data_dir, "GSE117221_H_top5000_prior_pairs.csv"),  row.names = FALSE)
write.csv(prior_TI, file.path(data_dir, "GSE117221_TI_top5000_prior_pairs.csv"), row.names = FALSE)

# ---------- 9. Output summary ----------
cat("\nPart 2 complete!\n\n")
cat("H group:\n")
cat("Original genes =", nrow(expr_H), "\n")
cat("After low-expression filtering =", nrow(res_H$filtered_count_mat), "\n")
cat("Top genes =", length(res_H$top_genes), "\n")
cat("Prior pairs =", nrow(prior_H), "\n\n")

cat("TI group:\n")
cat("Original genes =", nrow(expr_TI), "\n")
cat("After low-expression filtering =", nrow(res_TI$filtered_count_mat), "\n")
cat("Top genes =", length(res_TI$top_genes), "\n")
cat("Prior pairs =", nrow(prior_TI), "\n\n")

cat("Files saved to:\n", data_dir, "\n")


# ==============================================================================
# PART 3: Generate group-specific TF lists from top 5000 genes
# ==============================================================================

# ---------- 1. Load required packages ----------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)

# ---------- 2. File paths ----------
file_H_notrans  <- file.path(data_dir, "GSE117221_H_top5000_notransposed.csv")
file_TI_notrans <- file.path(data_dir, "GSE117221_TI_top5000_notransposed.csv")

tf_file <- "E:/SCD/数据/构建网络的数据/human_tf_list.txt"

# ---------- 3. Read non-transposed top5000 matrices ----------
read_expr_file <- function(file_path) {
  df <- fread(file_path, data.table = FALSE, check.names = FALSE)
  colnames(df)[1] <- "Gene"
  df <- df[!is.na(df$Gene) & df$Gene != "", , drop = FALSE]
  df <- df[!duplicated(df$Gene), , drop = FALSE]
  return(df)
}

expr_H  <- read_expr_file(file_H_notrans)
expr_TI <- read_expr_file(file_TI_notrans)

# ---------- 4. Read TF list ----------
tf_df <- fread(tf_file, data.table = FALSE)
colnames(tf_df)[1] <- "TF"
tf_df <- tf_df[!is.na(tf_df$TF) & tf_df$TF != "", , drop = FALSE]
tf_df <- unique(tf_df)
tf_list <- tf_df$TF

# ---------- 5. Extract top5000 genes ----------
genes_H  <- expr_H$Gene
genes_TI <- expr_TI$Gene

# ---------- 6. Build group-specific TF lists ----------
tf_H  <- data.frame(TF = intersect(tf_list, genes_H),  stringsAsFactors = FALSE)
tf_TI <- data.frame(TF = intersect(tf_list, genes_TI), stringsAsFactors = FALSE)

tf_H  <- tf_H[order(tf_H$TF), , drop = FALSE]
tf_TI <- tf_TI[order(tf_TI$TF), , drop = FALSE]

# ---------- 7. Save TF lists ----------
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

# ---------- 8. Output summary ----------
cat("\nPart 3 complete!\n\n")
cat("Group-specific TF counts:\n")
cat("H  :", nrow(tf_H),  "\n")
cat("TI :", nrow(tf_TI), "\n\n")
cat("Files saved to:\n", data_dir, "\n")


# ==============================================================================
# PART 4: Build gene ID mappings and convert prior pairs to IDs
# ==============================================================================

# ---------- 1. Load required packages ----------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
library(data.table)

# ---------- 2. File paths ----------
file_H  <- file.path(data_dir, "GSE117221_H_top5000_notransposed.csv")
file_TI <- file.path(data_dir, "GSE117221_TI_top5000_notransposed.csv")

# New prior file with three columns (TF, Target, Type)
prior_file_new <- "C:/Users/Administrator/Desktop/new_GRN_Lung_GEN_counts_genename_c1.csv"

# ---------- 3. Read expression matrices ----------
expr_H  <- read_expr_file(file_H)
expr_TI <- read_expr_file(file_TI)

# ---------- 4. Build gene ID tables ----------
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

# Save gene ID files
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

# ---------- 5. Read new prior file ----------
prior_df_new <- fread(prior_file_new, data.table = FALSE, header = FALSE)

if (ncol(prior_df_new) < 3) {
  stop("Prior file has less than 3 columns. Please check file format.")
}

prior_df_new <- prior_df_new[, 1:3, drop = FALSE]
colnames(prior_df_new) <- c("Gene1", "Gene2", "Type")

# Remove missing values
prior_df_new <- prior_df_new[
  !is.na(prior_df_new$Gene1) & prior_df_new$Gene1 != "" &
    !is.na(prior_df_new$Gene2) & prior_df_new$Gene2 != "",
  ,
  drop = FALSE
]

# Remove duplicates
prior_df_new <- unique(prior_df_new)

# ---------- 6. Build group-specific prior datasets ----------
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

prior_H_new  <- build_dataset_prior(prior_df_new, expr_H$Gene)
prior_TI_new <- build_dataset_prior(prior_df_new, expr_TI$Gene)

# Save prior datasets
write.csv(prior_H_new,
          file.path(data_dir, "GSE117221_H_prior_dataset.csv"),
          row.names = FALSE)

write.csv(prior_TI_new,
          file.path(data_dir, "GSE117221_TI_prior_dataset.csv"),
          row.names = FALSE)

# ---------- 7. Convert prior pairs to gene IDs ----------
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

prior_H_ids  <- convert_prior_to_ids(prior_H_new,  gene_id_H)
prior_TI_ids <- convert_prior_to_ids(prior_TI_new, gene_id_TI)

# Save ID-only prior files (TSV, no header)
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

# ---------- 8. Convert prior pairs to IDs with Type column ----------
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

prior_H_ids_type  <- convert_prior_to_ids_with_type(prior_H_new,  gene_id_H)
prior_TI_ids_type <- convert_prior_to_ids_with_type(prior_TI_new, gene_id_TI)

# Save ID prior files with Type (TSV, with header)
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

# ---------- 9. Build verification tables ----------
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

check_H  <- build_check_table(prior_H_new,  gene_id_H)
check_TI <- build_check_table(prior_TI_new, gene_id_TI)

# Save check tables
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

# ---------- 10. Output summary ----------
cat("\nPart 4 complete!\n\n")

cat("Gene ID files:\n")
cat("H  :", nrow(gene_id_H),  "genes\n")
cat("TI :", nrow(gene_id_TI), "genes\n\n")

cat("Prior datasets:\n")
cat("H  :", nrow(prior_H_new),  "pairs\n")
cat("TI :", nrow(prior_TI_new), "pairs\n\n")

cat("ID version prior files:\n")
cat("H  :", nrow(prior_H_ids),  "pairs\n")
cat("TI :", nrow(prior_TI_ids), "pairs\n\n")

cat("All results saved to:\n", data_dir, "\n")