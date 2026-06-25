# ==============================================================================
# Script: GSE117221 Complete Data Processing Pipeline (H and TI groups only)
# 
# Description:
# This script processes expression count data from GSE117221 for two sample 
# groups (H and TI) through a comprehensive pipeline:
#
# Part 1: Filter top 10000 genes by mean expression and generate outputs
#   - Read merged count matrices with robust removal of summary rows
#   - Filter low-expression genes by CPM
#   - Select top 10000 genes by log2(CPM+1) mean expression
#   - Save transposed matrices, gene lists, mean expression tables, and prior pairs
#
# Part 2: Restore non-transposed matrices and build group-specific TF lists
#   - Reconstruct non-transposed matrices from transposed top10000 files
#   - Generate group-specific TF lists using human TF reference
#
# Part 3: Build gene ID mappings and convert prior pairs to IDs
#   - Create gene ID files (0-based indexing)
#   - Filter prior regulatory pairs for each dataset
#   - Convert gene names to IDs with and without Type information
#   - Generate verification tables
#
# Input data:
#   - Merged count matrices: GSE117221_H_counts_merged.txt, GSE117221_TI_counts_merged.txt
#   - Prior regulatory network 1 (two columns): label_c1.csv
#   - Prior regulatory network 2 (three columns): new_GRN_Lung_GEN_counts_genename_c1.csv
#   - Human TF list: human_tf_list.txt
#
# Output data (saved in the same directory):
#   Part 1 outputs:
#     * GSE117221_H_transposed.csv, GSE117221_TI_transposed.csv
#     * GSE117221_H_top10000_transposed.csv, GSE117221_TI_top10000_transposed.csv
#     * GSE117221_H_top10000_genes.csv, GSE117221_TI_top10000_genes.csv
#     * GSE117221_H_gene_mean_expression.csv, GSE117221_TI_gene_mean_expression.csv
#     * GSE117221_H_top10000_prior_pairs.csv, GSE117221_TI_top10000_prior_pairs.csv
#   Part 2 outputs:
#     * GSE117221_H_top10000_notransposed.csv, GSE117221_TI_top10000_notransposed.csv
#     * GSE117221_H_top10000_TF_list.txt, GSE117221_TI_top10000_TF_list.txt
#   Part 3 outputs:
#     * GSE117221_H_gene_ids.txt, GSE117221_TI_gene_ids.txt
#     * GSE117221_H_prior_dataset.csv, GSE117221_TI_prior_dataset.csv
#     * GSE117221_H_prior_ids.tsv, GSE117221_TI_prior_ids.tsv
#     * GSE117221_H_prior_ids_with_type.tsv, GSE117221_TI_prior_ids_with_type.tsv
#     * GSE117221_H_prior_check_table.tsv, GSE117221_TI_prior_check_table.tsv
# ==============================================================================

# ---------- 1. Load required packages ----------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("edgeR", quietly = TRUE)) install.packages("edgeR")

library(data.table)
library(edgeR)

# ---------- 2. Set file paths ----------
data_dir <- "E:/SCD/其他数据/GSE117221"

# Three merged count matrices (only H and TI)
file_H  <- file.path(data_dir, "GSE117221_H_counts_merged.txt")
file_TI <- file.path(data_dir, "GSE117221_TI_counts_merged.txt")

# First prior regulatory pairs file (two columns)
prior_file_1 <- "E:/SCD/数据/构建网络的数据/label_c1.csv"

# TF list file
tf_file <- "E:/SCD/数据/构建网络的数据/human_tf_list.txt"

# Second prior information file (three columns: Gene1, Gene2, Type)
prior_file_2 <- "C:/Users/Administrator/Desktop/new_GRN_Lung_GEN_counts_genename_c1.csv"

# ---------- 3. Define summary row identification function ----------
# More robust than exact matching: supports spaces, case variations, underscores, etc.
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

# ---------- 4. Expression matrix cleaning function ----------
clean_expr_df <- function(df, dataset_name = "dataset") {
  colnames(df)[1] <- "Gene"
  
  # Convert to character and trim whitespace
  df$Gene <- trimws(as.character(df$Gene))
  
  # Remove empty values
  df <- df[!is.na(df$Gene) & df$Gene != "", , drop = FALSE]
  
  # Identify summary rows
  bad_idx <- is_bad_feature(df$Gene)
  removed_genes <- unique(df$Gene[bad_idx])
  
  if (length(removed_genes) > 0) {
    cat("\n", dataset_name, " - removed non-gene rows:\n", sep = "")
    print(removed_genes)
  } else {
    cat("\n", dataset_name, " - no summary rows detected.\n", sep = "")
  }
  
  # Remove summary rows
  df <- df[!bad_idx, , drop = FALSE]
  
  # Remove duplicates
  dup_genes <- unique(df$Gene[duplicated(df$Gene)])
  if (length(dup_genes) > 0) {
    cat(dataset_name, " - duplicate genes detected, keeping first occurrence. Examples:\n")
    print(head(dup_genes, 20))
  }
  df <- df[!duplicated(df$Gene), , drop = FALSE]
  
  rownames(df) <- NULL
  return(df)
}

# =========================================================
# Part 1: Filter top 10000 genes by mean expression
# =========================================================

# ---------- 5. Read expression matrices (cleaned version) ----------
read_expr_matrix <- function(file_path, dataset_name = "dataset") {
  df <- fread(file_path, data.table = FALSE)
  df <- clean_expr_df(df, dataset_name = dataset_name)
  
  rownames(df) <- df$Gene
  df$Gene <- NULL
  
  mat <- as.matrix(df)
  mode(mat) <- "numeric"
  
  return(mat)   # rows=genes, columns=samples
}

expr_H  <- read_expr_matrix(file_H,  dataset_name = "H merged counts")
expr_TI <- read_expr_matrix(file_TI, dataset_name = "TI merged counts")

# ---------- 6. Define gene selection function ----------
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

# ---------- 7. Transpose matrices ----------
transpose_expr <- function(mat) {
  tmat <- as.data.frame(t(mat), check.names = FALSE)
  tmat <- cbind(Sample = rownames(tmat), tmat)
  rownames(tmat) <- NULL
  return(tmat)
}

expr_H_t  <- transpose_expr(expr_H)
expr_TI_t <- transpose_expr(expr_TI)

top_H_t  <- transpose_expr(res_H$top_count_mat)
top_TI_t <- transpose_expr(res_TI$top_count_mat)

# ---------- 8. Read first prior regulatory pairs ----------
prior_df_1 <- fread(prior_file_1, data.table = FALSE)
colnames(prior_df_1)[1:2] <- c("Gene1", "Gene2")
prior_df_1$Gene1 <- trimws(as.character(prior_df_1$Gene1))
prior_df_1$Gene2 <- trimws(as.character(prior_df_1$Gene2))
prior_df_1 <- prior_df_1[!is.na(prior_df_1$Gene1) & !is.na(prior_df_1$Gene2), , drop = FALSE]
prior_df_1 <- unique(prior_df_1)

# ---------- 9. Build dataset-specific prior pairs ----------
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

# ---------- 10. Save Part 1 results ----------
write.csv(expr_H_t,  file.path(data_dir, "GSE117221_H_transposed.csv"),  row.names = FALSE)
write.csv(expr_TI_t, file.path(data_dir, "GSE117221_TI_transposed.csv"), row.names = FALSE)

write.csv(top_H_t,  file.path(data_dir, "GSE117221_H_top10000_transposed.csv"),  row.names = FALSE)
write.csv(top_TI_t, file.path(data_dir, "GSE117221_TI_top10000_transposed.csv"), row.names = FALSE)

write.csv(data.frame(Gene = res_H$top_genes),  file.path(data_dir, "GSE117221_H_top10000_genes.csv"),  row.names = FALSE)
write.csv(data.frame(Gene = res_TI$top_genes), file.path(data_dir, "GSE117221_TI_top10000_genes.csv"), row.names = FALSE)

write.csv(res_H$gene_mean_table,  file.path(data_dir, "GSE117221_H_gene_mean_expression.csv"),  row.names = FALSE)
write.csv(res_TI$gene_mean_table, file.path(data_dir, "GSE117221_TI_gene_mean_expression.csv"), row.names = FALSE)

write.csv(prior_H_1,  file.path(data_dir, "GSE117221_H_top10000_prior_pairs.csv"),  row.names = FALSE)
write.csv(prior_TI_1, file.path(data_dir, "GSE117221_TI_top10000_prior_pairs.csv"), row.names = FALSE)

# ---------- 11. Output Part 1 summary ----------
cat("\nPart 1 processing complete!\n\n")

cat("H group:\n")
cat("Genes after summary row removal =", nrow(expr_H), "\n")
cat("Genes after low-expression filtering =", nrow(res_H$filtered_count_mat), "\n")
cat("Top genes =", length(res_H$top_genes), "\n")
cat("Prior pairs =", nrow(prior_H_1), "\n\n")

cat("TI group:\n")
cat("Genes after summary row removal =", nrow(expr_TI), "\n")
cat("Genes after low-expression filtering =", nrow(res_TI$filtered_count_mat), "\n")
cat("Top genes =", length(res_TI$top_genes), "\n")
cat("Prior pairs =", nrow(prior_TI_1), "\n\n")

cat("Part 1 files saved to:\n", data_dir, "\n\n")


# =========================================================
# Part 2: Restore non-transposed matrices and build TF lists
# =========================================================

# ---------- 12. File paths ----------
file_H_t  <- file.path(data_dir, "GSE117221_H_top10000_transposed.csv")
file_TI_t <- file.path(data_dir, "GSE117221_TI_top10000_transposed.csv")

# ---------- 13. Read transposed expression matrices ----------
read_transposed_expr <- function(file_path) {
  df <- fread(file_path, data.table = FALSE, check.names = FALSE)
  colnames(df)[1] <- "Sample"
  return(df)
}

expr_H_t2  <- read_transposed_expr(file_H_t)
expr_TI_t2 <- read_transposed_expr(file_TI_t)

# ---------- 14. Restore to non-transposed format ----------
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

# Clean again to ensure output matrices are absolutely clean
expr_H_notrans  <- clean_expr_df(expr_H_notrans,  dataset_name = "H top10000 notransposed")
expr_TI_notrans <- clean_expr_df(expr_TI_notrans, dataset_name = "TI top10000 notransposed")

# ---------- 15. Save restored non-transposed matrices ----------
write.csv(expr_H_notrans,
          file.path(data_dir, "GSE117221_H_top10000_notransposed.csv"),
          row.names = FALSE)

write.csv(expr_TI_notrans,
          file.path(data_dir, "GSE117221_TI_top10000_notransposed.csv"),
          row.names = FALSE)

# ---------- 16. Read TF list ----------
tf_df <- fread(tf_file, data.table = FALSE)
colnames(tf_df)[1] <- "TF"
tf_df$TF <- trimws(as.character(tf_df$TF))
tf_df <- tf_df[!is.na(tf_df$TF) & tf_df$TF != "", , drop = FALSE]
tf_df <- unique(tf_df)
tf_list <- tf_df$TF

# ---------- 17. Extract top10000 genes from each dataset ----------
genes_H  <- expr_H_notrans$Gene
genes_TI <- expr_TI_notrans$Gene

# ---------- 18. Build group-specific TF lists ----------
tf_H  <- data.frame(TF = intersect(tf_list, genes_H),  stringsAsFactors = FALSE)
tf_TI <- data.frame(TF = intersect(tf_list, genes_TI), stringsAsFactors = FALSE)

tf_H  <- tf_H[order(tf_H$TF), , drop = FALSE]
tf_TI <- tf_TI[order(tf_TI$TF), , drop = FALSE]

# ---------- 19. Save TF lists ----------
write.table(tf_H,
            file.path(data_dir, "GSE117221_H_top10000_TF_list.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(tf_TI,
            file.path(data_dir, "GSE117221_TI_top10000_TF_list.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ---------- 20. Output Part 2 summary ----------
cat("Part 2 processing complete!\n\n")

cat("Restored non-transposed expression matrices:\n")
cat("H  :", dim(expr_H_notrans)[1],  "genes x", dim(expr_H_notrans)[2] - 1, "samples\n")
cat("TI :", dim(expr_TI_notrans)[1], "genes x", dim(expr_TI_notrans)[2] - 1, "samples\n\n")

cat("Group-specific TF counts:\n")
cat("H  :", nrow(tf_H), "\n")
cat("TI :", nrow(tf_TI), "\n\n")

cat("Part 2 files saved to:\n", data_dir, "\n\n")


# =========================================================
# Part 3: Build gene ID mappings and convert prior pairs to IDs
# =========================================================

# ---------- 21. File paths ----------
file_H_notrans  <- file.path(data_dir, "GSE117221_H_top10000_notransposed.csv")
file_TI_notrans <- file.path(data_dir, "GSE117221_TI_top10000_notransposed.csv")

# ---------- 22. Read expression matrices (cleaned) ----------
read_expr_file <- function(file_path, dataset_name = "dataset") {
  df <- fread(file_path, data.table = FALSE, check.names = FALSE)
  df <- clean_expr_df(df, dataset_name = dataset_name)
  return(df)
}

expr_H3  <- read_expr_file(file_H_notrans,  dataset_name = "H final notransposed")
expr_TI3 <- read_expr_file(file_TI_notrans, dataset_name = "TI final notransposed")

# ---------- 23. Build gene ID files ----------
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

write.table(gene_id_H,
            file.path(data_dir, "GSE117221_H_gene_ids.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(gene_id_TI,
            file.path(data_dir, "GSE117221_TI_gene_ids.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ---------- 24. Read second prior file ----------
prior_df_2 <- fread(prior_file_2, data.table = FALSE, header = FALSE)

if (ncol(prior_df_2) < 3) {
  stop("Prior file has less than 3 columns. Please check file format.")
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

# ---------- 25. Build dataset-specific prior datasets ----------
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

write.csv(prior_H_2,
          file.path(data_dir, "GSE117221_H_prior_dataset.csv"),
          row.names = FALSE)

write.csv(prior_TI_2,
          file.path(data_dir, "GSE117221_TI_prior_dataset.csv"),
          row.names = FALSE)

# ---------- 26. Convert prior pairs to gene IDs ----------
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

write.table(prior_H_ids,
            file.path(data_dir, "GSE117221_H_prior_ids.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(prior_TI_ids,
            file.path(data_dir, "GSE117221_TI_prior_ids.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ---------- 27. Convert prior pairs to IDs with Type column ----------
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

write.table(prior_H_ids_type,
            file.path(data_dir, "GSE117221_H_prior_ids_with_type.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(prior_TI_ids_type,
            file.path(data_dir, "GSE117221_TI_prior_ids_with_type.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ---------- 28. Build verification tables ----------
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

write.table(check_H,
            file.path(data_dir, "GSE117221_H_prior_check_table.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(check_TI,
            file.path(data_dir, "GSE117221_TI_prior_check_table.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ---------- 29. Output Part 3 summary ----------
cat("Part 3 processing complete!\n\n")

cat("Gene ID files:\n")
cat("H  :", nrow(gene_id_H),  "genes\n")
cat("TI :", nrow(gene_id_TI), "genes\n\n")

cat("Dataset-specific prior datasets:\n")
cat("H  :", nrow(prior_H_2),  "pairs\n")
cat("TI :", nrow(prior_TI_2), "pairs\n\n")

cat("ID version prior files:\n")
cat("H  :", nrow(prior_H_ids),  "pairs\n")
cat("TI :", nrow(prior_TI_ids), "pairs\n\n")

cat("All results saved to:\n", data_dir, "\n")