# ==============================================================================
# Script: GSE254951 pre-HU Data Processing Pipeline
# 
# Description:
# This script processes the GSE254951 raw count data to:
#   1. Filter and extract pre-HU samples (excluding washout samples HU21 and HU35)
#   2. Convert Ensembl IDs to gene symbols using biomaRt
#   3. Select the top 10,000 most highly expressed genes based on total counts
#   4. Generate transposed expression matrix (samples as rows, genes as columns)
#   5. Filter prior regulatory pairs to include only top 10,000 genes
#
# Input data:
#   - Expression file: GSE254951_geo_rawcounts.txt
#     (Rows: Ensembl IDs, Columns: HU samples)
#   - Prior regulatory network: label_c1.csv
#     (Two columns: Gene1, Gene2)
#
# Output data (saved in the same output directory):
#   - GSE254951_preHU_Top10000_counts.csv: Top 10,000 gene expression matrix
#     (Rows: gene symbols, Columns: pre-HU samples)
#   - GSE254951_preHU_Top10000_counts_transposed.csv: Transposed matrix
#     (Rows: samples, Columns: gene symbols)
#   - label_c1_prior_filtered_by_GSE254951_preHU_Top10000.csv: Filtered prior pairs
#   - GSE254951_ensembl_to_symbol_mapping.csv: Ensembl ID to gene symbol mapping
# ==============================================================================

library(data.table)
library(dplyr)
library(stringr)
library(biomaRt)

# ------------------------------------------------------------------------------
# Set file paths
# ------------------------------------------------------------------------------
expr_file  <- "D:/放假/SCD/其他数据/GSE254951/GSE254951_geo_rawcounts.txt"
prior_file <- "D:/放假/SCD/数据/构建网络的数据/label_c1.csv"
out_dir    <- "D:/放假/SCD/其他数据/GSE254951"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_expr_top10k_csv   <- file.path(out_dir, "GSE254951_preHU_Top10000_counts.csv")
out_expr_top10k_T_csv <- file.path(out_dir, "GSE254951_preHU_Top10000_counts_transposed.csv")
out_prior_top10k_csv  <- file.path(out_dir, "label_c1_prior_filtered_by_GSE254951_preHU_Top10000.csv")
out_idmap_csv         <- file.path(out_dir, "GSE254951_ensembl_to_symbol_mapping.csv")

# ------------------------------------------------------------------------------
# Step 1) Read raw counts and filter pre-HU samples
# ------------------------------------------------------------------------------
df <- fread(expr_file, data.table = FALSE, check.names = FALSE)

colnames(df)[1] <- "ensembl_id"

hu_cols <- setdiff(colnames(df), "ensembl_id")
hu_num  <- as.integer(str_remove(hu_cols, "^HU"))

# pre-HU = odd-numbered columns, excluding washout samples HU21 and HU35
prehu_cols <- hu_cols[(hu_num %% 2 == 1) & !(hu_num %in% c(21, 35))]

expr_prehu <- df[, c("ensembl_id", prehu_cols)]

# ------------------------------------------------------------------------------
# Step 2) ID conversion: Ensembl -> Gene Symbol
# ------------------------------------------------------------------------------
expr_prehu$ensembl_id_nover <- str_remove(expr_prehu$ensembl_id, "\\.\\d+$")

mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

map <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = unique(expr_prehu$ensembl_id_nover),
  mart       = mart
)

write.csv(map, out_idmap_csv, row.names = FALSE, quote = TRUE)

# Merge mapping
expr2 <- expr_prehu %>%
  dplyr::select(-ensembl_id) %>%
  dplyr::left_join(map, by = c("ensembl_id_nover" = "ensembl_gene_id")) %>%
  dplyr::rename(GeneSymbol = hgnc_symbol)

# Remove rows without gene symbols
expr2 <- expr2 %>% dplyr::filter(!is.na(GeneSymbol), GeneSymbol != "")

# If multiple Ensembl IDs map to the same symbol, sum their counts
expr2_num <- expr2 %>%
  dplyr::select(GeneSymbol, all_of(prehu_cols))

expr2_num[prehu_cols] <- lapply(expr2_num[prehu_cols], as.numeric)

expr_symbol <- expr2_num %>%
  dplyr::group_by(GeneSymbol) %>%
  dplyr::summarise(across(all_of(prehu_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# ------------------------------------------------------------------------------
# Step 3) Select top 10,000 genes and transpose
# ------------------------------------------------------------------------------
mat <- as.matrix(expr_symbol[, prehu_cols])
mat[is.na(mat)] <- 0

expr_symbol$TotalCount <- rowSums(mat)

expr_top10k <- expr_symbol %>%
  dplyr::arrange(desc(TotalCount)) %>%
  dplyr::slice_head(n = 10000) %>%
  dplyr::select(-TotalCount)

write.csv(expr_top10k, out_expr_top10k_csv, row.names = FALSE, quote = TRUE)

# Transpose: samples x genes
mat10k <- as.matrix(expr_top10k[, prehu_cols])
rownames(mat10k) <- expr_top10k$GeneSymbol
mat10k[is.na(mat10k)] <- 0

mat10k_t <- t(mat10k)
out_t <- data.frame(Sample = rownames(mat10k_t), mat10k_t, check.names = FALSE)
write.csv(out_t, out_expr_top10k_T_csv, row.names = FALSE, quote = TRUE)

gene_set_10k <- expr_top10k$GeneSymbol

# ------------------------------------------------------------------------------
# Step 4) Read prior pairs and filter
# ------------------------------------------------------------------------------
read_prior_pairs <- function(path) {
  pri <- tryCatch(
    read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(pri) || ncol(pri) < 2) {
    pri <- read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  pri <- pri[, 1:2]
  colnames(pri) <- c("Gene1", "Gene2")
  pri <- pri %>%
    dplyr::filter(!is.na(Gene1), !is.na(Gene2), Gene1 != "", Gene2 != "")
  pri
}

prior <- read_prior_pairs(prior_file)

prior_filt <- prior %>%
  dplyr::filter(Gene1 %in% gene_set_10k, Gene2 %in% gene_set_10k) %>%
  dplyr::distinct()

write.csv(prior_filt, out_prior_top10k_csv, row.names = FALSE, quote = TRUE)

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
cat("=== GSE254951 pre-HU processing complete ===\n")
cat("pre-HU sample count:", length(prehu_cols), "\n")
cat("pre-HU sample columns:\n", paste(prehu_cols, collapse = ", "), "\n\n")
cat("Top 10,000 genes:", nrow(expr_top10k), "\n")
cat("Filtered prior edges:", nrow(prior_filt), "\n\n")
cat("Output expression matrix:", out_expr_top10k_csv, "\n")
cat("Output transposed matrix:", out_expr_top10k_T_csv, "\n")
cat("Output prior edges:", out_prior_top10k_csv, "\n")