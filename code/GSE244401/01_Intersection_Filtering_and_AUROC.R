# ==============================================================================
# Script: Extract Edges Predicted by ≥2 Methods and Validate with ChIP-seq
# 
# Description:
# This script processes predictions from three gene regulatory network inference 
# methods (DeepSEM, DeepFGRN, 3DCEMA) and:
#   1. Extracts edges that are predicted by at least two methods (consensus edges)
#   2. Computes the average EdgeWeight for consensus edges
#   3. Validates these consensus edges using ChIP-seq peaks for multiple TFs
#   4. Calculates AUROC and AUPR for each TF-specific validation
# 
# Input data:
#   - Method prediction files (all edges with EdgeWeight):
#     * deepsem_predicted_all_edges_named_c1.txt
#     * deepfgrn_predicted_all_edges_named_c1.txt
#     * 3dcema_predicted_all_edges_named_c1.txt
#   - ChIP-seq BED files for each TF:
#     * GSM970257_GATA1-F_peaks.bed
#     * GSM1816086_NFE2-F5.peak.bed
#     * GSM1816083_TAL1-F5.peak.bed
#     * GSE103445_HUDEP2_BCL11A_optimal_idr_peaks.bed
# 
# Output data:
#   - two_or_more_overlap_c1.txt: Consensus edges predicted by ≥2 methods
#   - Individual TF evaluation directories containing:
#     * ROC_curve.png: ROC curve visualization
#     * PR_curve.png: Precision-Recall curve visualization
#     * summary.txt: Detailed evaluation metrics
#   - auc_ap_summary_all.txt: Combined summary for all TFs
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Clean environment and load packages
# ------------------------------------------------------------------------------
rm(list = ls()); gc()
if (!requireNamespace("data.table", quietly=TRUE))   install.packages("data.table")
if (!requireNamespace("dplyr",      quietly=TRUE))   install.packages("dplyr")
if (!requireNamespace("ChIPseeker", quietly=TRUE))   BiocManager::install("ChIPseeker")
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly=TRUE))
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
if (!requireNamespace("org.Hs.eg.db", quietly=TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("pROC",    quietly=TRUE))       install.packages("pROC")
if (!requireNamespace("PRROC",   quietly=TRUE))       install.packages("PRROC")
if (!requireNamespace("ggplot2", quietly=TRUE))       install.packages("ggplot2")

library(data.table)
library(dplyr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(pROC)
library(PRROC)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Define cleaning function
# ------------------------------------------------------------------------------
clean_name <- function(x) {
  x %>%
    trimws() %>%
    tolower() %>%
    gsub("[^a-z0-9]", "", .)
}

# ------------------------------------------------------------------------------
# 2. Read three method prediction files and build dt_list
# ------------------------------------------------------------------------------
paths <- list(
  deepsem  = "E:/1数据/deepsem/filtered_networks/deepsem_predicted_all_edges_named_c1.txt",
  deepfgrn = "E:/1数据/deepfgrn/filtered_networks/deepfgrn_predicted_all_edges_named_c1.txt",
  dcema    = "E:/1数据/3dcema/filtered_networks/3dcema_predicted_all_edges_named_c1.txt"
)

dt_list <- lapply(names(paths), function(meth) {
  dt <- fread(paths[[meth]], sep="\t", header=TRUE)
  dt[, TF     := clean_name(TF)]
  dt[, Target := clean_name(Target)]
  dt[, pair   := paste0(TF, "_", Target)]
  # Rename EdgeWeight to method-specific column
  ew_col <- paste0("EW_", meth)
  setnames(dt, "EdgeWeight", ew_col)
  # Keep only required columns
  keep_cols <- c("pair", "TF", "Target", ew_col)
  dt <- dt[, ..keep_cols]
  # If duplicate pairs exist within the same method, aggregate by mean
  dt <- dt[, .(
    TF = TF[1L],
    Target = Target[1L],
    tmp = mean(get(ew_col), na.rm = TRUE)
  ), by = pair]
  setnames(dt, "tmp", ew_col)
  dt
})
names(dt_list) <- names(paths)

# ------------------------------------------------------------------------------
# 3. Build complete pairs_map
# ------------------------------------------------------------------------------
pairs_map <- rbindlist(
  lapply(dt_list, function(dt) dt[, .(pair, TF, Target)]),
  use.names = TRUE
)
pairs_map <- unique(pairs_map, by = "pair")

# ------------------------------------------------------------------------------
# 4. Merge weights from all methods
# ------------------------------------------------------------------------------
merged_dt <- copy(pairs_map)
for (meth in names(dt_list)) {
  ew_col <- paste0("EW_", meth)
  merged_dt <- merge(
    merged_dt,
    dt_list[[meth]][, c("pair", ew_col), with = FALSE],
    by   = "pair",
    all.x = TRUE
  )
}

# ------------------------------------------------------------------------------
# 5. Filter edges with ≥2 methods and calculate average weight
# ------------------------------------------------------------------------------
ew_cols <- paste0("EW_", names(paths))
merged_dt[, count_present := rowSums(!is.na(.SD)), .SDcols = ew_cols]
two_or_more <- merged_dt[count_present >= 2]
two_or_more[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = ew_cols]
two_or_more <- two_or_more[, .(TF, Target, EdgeWeight)]

# ------------------------------------------------------------------------------
# 6. Write consensus subset
# ------------------------------------------------------------------------------
out_pred <- "E:/1数据/filtered_networks/two_or_more_overlap_c1.txt"
fwrite(two_or_more, out_pred, sep="\t", quote=FALSE)

# ------------------------------------------------------------------------------
# 7. Define ChIP annotation function
# ------------------------------------------------------------------------------
get_chip_pairs <- function(bed_path, tf_name, tss_region = c(-2000, 500)) {
  lines <- readLines(bed_path)
  lines <- lines[!grepl("^(track|browser|#)", lines)]
  tmp_bed <- tempfile(fileext = ".bed")
  writeLines(lines, tmp_bed)
  pa <- annotatePeak(
    tmp_bed,
    tssRegion = tss_region,
    TxDb      = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb    = "org.Hs.eg.db"
  )
  anno_df <- as.data.frame(pa)
  genes   <- unique(na.omit(anno_df$SYMBOL))
  paste0(clean_name(tf_name), "_", clean_name(genes))
}

# ------------------------------------------------------------------------------
# 8. Define evaluation function: AUROC/AUPR and plotting
# ------------------------------------------------------------------------------
evaluate_network <- function(pred_dt, chip_pairs, out_dir) {
  df <- copy(pred_dt)[
    , pair := paste0(TF, "_", Target)
  ][
    , `:=`(
      label = as.integer(pair %in% chip_pairs),
      score = EdgeWeight
    )
  ]
  if (length(unique(df$label)) < 2) return(NULL)
  roc_obj <- roc(df$label, df$score, quiet=TRUE)
  pr_obj  <- pr.curve(
    scores.class0 = df$score[df$label == 1],
    scores.class1 = df$score[df$label == 0],
    curve = TRUE
  )
  # Output plots
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  png(file.path(out_dir, "ROC_curve.png"), 600,600)
  plot(roc_obj, main=sprintf("ROC (AUROC=%.3f)", auc(roc_obj)))
  dev.off()
  png(file.path(out_dir, "PR_curve.png"), 600,600)
  plot(pr_obj$curve[,1], pr_obj$curve[,2],
       type="l", xlab="Recall", ylab="Precision",
       main=sprintf("PR (AUPR=%.3f)", pr_obj$auc.integral))
  dev.off()
  data.frame(
    Total = nrow(df),
    Pos   = sum(df$label),
    Neg   = sum(1 - df$label),
    AUROC = round(auc(roc_obj),   4),
    AUPR  = round(pr_obj$auc.integral, 4),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 9. Batch ChIP-seq validation
# ------------------------------------------------------------------------------
tf_beds <- list(
  GATA1  = "D:/Users/29321/Desktop/GSM970257_GATA1-F_peaks.bed",
  NFE2   = "D:/Users/29321/Desktop/GSM1816086_NFE2-F5.peak.bed",
  TAL1   = "D:/Users/29321/Desktop/GSM1816083_TAL1-F5.peak.bed",
  BCL11A = "D:/Users/29321/Desktop/GSE103445_HUDEP2_BCL11A_optimal_idr_peaks.bed"
)

results <- list()
for (tf in names(tf_beds)) {
  cat(">> Validating TF:", tf, "\n")
  chip_pairs <- get_chip_pairs(tf_beds[[tf]], tf)
  res <- evaluate_network(two_or_more, chip_pairs,
                          out_dir = file.path("E:/1数据/evaluation", tf))
  if (is.null(res)) {
    cat("  Warning: Single label, cannot compute ROC/PR.\n")
  } else {
    results[[tf]] <- cbind(TF = tf, res)
    print(results[[tf]])
    write.table(results[[tf]],
                file   = file.path("E:/1数据/evaluation", tf, "summary.txt"),
                sep    = "\t", quote = FALSE, row.names = FALSE)
  }
}

# ------------------------------------------------------------------------------
# 10. Summarize all results
# ------------------------------------------------------------------------------
summary_all <- rbindlist(lapply(results, as.data.frame), fill=TRUE)
fwrite(summary_all,
       "E:/1数据/evaluation/auc_ap_summary_all.txt",
       sep="\t", quote=FALSE)
cat("\nAll processing complete! Results saved in E:/1数据/evaluation.\n")























# -----------------------------------------------------------------------------
# Batch ChIP-seq Validation for Three Transcription Factors (GATA1, TAL1, BCL11A)
# -----------------------------------------------------------------------------
# Purpose:
# This script performs batch validation of predicted regulatory networks for three
# transcription factors (GATA1, TAL1, BCL11A) against ChIP-seq ground truth data.
# It evaluates four overlapping edge subsets: three pairwise overlaps and one 
# triple overlap across three prediction methods (DeepSEM, DeepFGRN, 3D-CEMA).
#
# Input:
#   1. Three prediction files (tab-separated) from DeepSEM, DeepFGRN, and 3D-CEMA
#      Each file must contain columns: TF, Target, EdgeWeight (or Score)
#   2. Three ChIP-seq BED files (one per transcription factor)
#      - GATA1:  GSM970257_GATA1-F_peaks.bed
#      - TAL1:   GSM1816083_TAL1-F5.peak.bed
#      - BCL11A: GSE103445_HUDEP2_BCL11A_optimal_idr_peaks.bed
#
# Output:
#   1. Four overlap subset files (tab-separated):
#      - sem_fgrn_overlap_s2.txt (DeepSEM ∩ DeepFGRN)
#      - sem_dcema_overlap_s2.txt (DeepSEM ∩ 3D-CEMA)
#      - fgrn_dcema_overlap_s2.txt (DeepFGRN ∩ 3D-CEMA)
#      - three_overlap_s2.txt (DeepSEM ∩ DeepFGRN ∩ 3D-CEMA)
#   2. For each TF and each subset:
#      - AUROC and AUPR curves (PNG format)
#      - Summary statistics file (tab-separated)
#   3. Console output showing evaluation metrics for each combination
#
# Note: All gene/peak names are cleaned by converting to lowercase and removing
#       non-alphanumeric characters to ensure consistent matching.
# -----------------------------------------------------------------------------

# 0. Set working directory and define input file paths for s2 condition
deepsem_file  <- "E:/1数据/deepsem/filtered_networks/deepsem_predicted_all_edges_named_s2.txt"
deepfgrn_file <- "E:/1数据/deepfgrn/filtered_networks/deepfgrn_predicted_all_edges_named_s2.txt"
dcema_file    <- "E:/1数据/3dcema/filtered_networks/3dcema_predicted_all_edges_named_s2.txt"

# 1. Define name cleaning function (convert TF/Target to lowercase, remove non-alphanumeric)
clean_name <- function(x) {
  x <- trimws(x)
  x <- tolower(x)
  gsub("[^a-z0-9]", "", x)
}

# 2. Read three prediction files and generate "pair" column
library(data.table)
df_deepsem  <- fread(deepsem_file,  header = TRUE, sep = "\t", data.table = FALSE)
df_deepfgrn <- fread(deepfgrn_file, header = TRUE, sep = "\t", data.table = FALSE)
df_dcema    <- fread(dcema_file,    header = TRUE, sep = "\t", data.table = FALSE)

df_deepsem$TF     <- clean_name(df_deepsem$TF)
df_deepsem$Target <- clean_name(df_deepsem$Target)
df_deepsem$pair   <- paste0(df_deepsem$TF, "_", df_deepsem$Target)

df_deepfgrn$TF     <- clean_name(df_deepfgrn$TF)
df_deepfgrn$Target <- clean_name(df_deepfgrn$Target)
df_deepfgrn$pair   <- paste0(df_deepfgrn$TF, "_", df_deepfgrn$Target)

df_dcema$TF     <- clean_name(df_dcema$TF)
df_dcema$Target <- clean_name(df_dcema$Target)
df_dcema$pair   <- paste0(df_dcema$TF, "_", df_dcema$Target)

pairs_deepsem  <- df_deepsem$pair
pairs_deepfgrn <- df_deepfgrn$pair
pairs_dcema    <- df_dcema$pair

# 3. Calculate triple overlap and three pairwise overlaps
three_overlap_pairs <- Reduce(intersect, list(pairs_deepsem, pairs_deepfgrn, pairs_dcema))
sem_fgrn_pairs     <- intersect(pairs_deepsem, pairs_deepfgrn)
sem_dcema_pairs    <- intersect(pairs_deepsem, pairs_dcema)
fgrn_dcema_pairs   <- intersect(pairs_deepfgrn, pairs_dcema)

# 4. Filter four edge subsets from df_deepsem, keep TF/Target/EdgeWeight
two_sem_fgrn_df   <- df_deepsem[df_deepsem$pair %in% sem_fgrn_pairs,  c("TF","Target","EdgeWeight"), drop = FALSE]
two_sem_dcema_df  <- df_deepsem[df_deepsem$pair %in% sem_dcema_pairs, c("TF","Target","EdgeWeight"), drop = FALSE]
two_fgrn_dcema_df <- df_deepsem[df_deepsem$pair %in% fgrn_dcema_pairs, c("TF","Target","EdgeWeight"), drop = FALSE]
three_overlap_df  <- df_deepsem[df_deepsem$pair %in% three_overlap_pairs, c("TF","Target","EdgeWeight"), drop = FALSE]

# 5. Write the four subset files
out_sem_fgrn   <- "sem_fgrn_overlap_s2.txt"
out_sem_dcema  <- "sem_dcema_overlap_s2.txt"
out_fgrn_dcema <- "fgrn_dcema_overlap_s2.txt"
out_three      <- "three_overlap_s2.txt"

fwrite(two_sem_fgrn_df,   file = out_sem_fgrn,   sep = "\t", quote = FALSE, col.names = TRUE)
fwrite(two_sem_dcema_df,  file = out_sem_dcema,  sep = "\t", quote = FALSE, col.names = TRUE)
fwrite(two_fgrn_dcema_df, file = out_fgrn_dcema, sep = "\t", quote = FALSE, col.names = TRUE)
fwrite(three_overlap_df,  file = out_three,      sep = "\t", quote = FALSE, col.names = TRUE)

# -----------------------------------------------------------------------------
# 6. Batch validation for three transcription factors (GATA1, TAL1, BCL11A)
# -----------------------------------------------------------------------------
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(pROC)
library(PRROC)
library(ggplot2)
library(dplyr)

# 6.1 Define get_chip_pairs(): Annotate BED file to "TF_Target" pair set
get_chip_pairs <- function(bed_path, tf_name, tss_region = c(-2000, 2000)) {
  bed_raw <- fread(bed_path, header = FALSE)[, 1:3]
  colnames(bed_raw) <- c("chr","start","end")
  
  bed_clean <- bed_raw %>%
    filter(!is.na(chr), !is.na(start), !is.na(end)) %>%
    filter(suppressWarnings(is.numeric(start)), suppressWarnings(is.numeric(end))) %>%
    filter(start <= end)
  
  if (nrow(bed_clean) == 0) {
    stop("❌ Cleaned BED file is empty: ", bed_path)
  }
  
  tmp_bed <- tempfile(fileext = ".bed")
  fwrite(bed_clean, tmp_bed, sep = "\t", col.names = FALSE, quote = FALSE)
  
  peakAnno <- annotatePeak(
    tmp_bed,
    tssRegion = tss_region,
    TxDb      = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb    = "org.Hs.eg.db"
  )
  
  genes <- as.data.frame(peakAnno@anno)$SYMBOL
  chip_pairs <- paste0(clean_name(tf_name), "_", clean_name(genes))
  unique(chip_pairs)
}

# 6.2 Define evaluate_network(): Evaluate single subset file, save plots & return stats
evaluate_network <- function(pred_file, chip_pairs_set, out_dir) {
  df <- fread(pred_file, header = TRUE, sep = "\t")
  if (!all(c("TF","Target") %in% colnames(df))) {
    warning("Skipping, missing TF/Target columns: ", pred_file)
    return(NULL)
  }
  score_col <- if      ("EdgeWeight" %in% colnames(df)) "EdgeWeight" 
  else if ("Score"      %in% colnames(df)) "Score"
  else {
    warning("Skipping, missing Score/EdgeWeight column: ", pred_file)
    return(NULL)
  }
  
  df$TF     <- clean_name(df$TF)
  df$Target <- clean_name(df$Target)
  pairs    <- paste0(df$TF,"_",df$Target)
  
  labels <- as.integer(pairs %in% chip_pairs_set)
  scores <- df[[score_col]]
  
  total <- length(labels)
  pos   <- sum(labels)
  neg   <- total - pos
  
  if (length(unique(labels)) < 2) {
    auc_roc <- NA
    auc_pr  <- NA
  } else {
    roc_obj <- roc(labels, scores, quiet = TRUE)
    auc_roc <- as.numeric(auc(roc_obj))
    
    pr_obj <- pr.curve(
      scores.class0 = scores[labels == 1],
      scores.class1 = scores[labels == 0],
      curve = TRUE
    )
    auc_pr <- pr_obj$auc.integral
    
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    png(file.path(out_dir, paste0("ROC_", basename(pred_file), ".png")), 600, 600)
    plot(roc_obj, main = paste("ROC", basename(pred_file)), lwd = 2)
    dev.off()
    
    png(file.path(out_dir, paste0("PR_", basename(pred_file), ".png")), 600, 600)
    plot(pr_obj$curve[,1], pr_obj$curve[,2],
         type = "l", xlab = "Recall", ylab = "Precision",
         main = paste("PR", basename(pred_file)), lwd = 2)
    dev.off()
  }
  
  data.frame(
    File           = basename(pred_file),
    Total_Edges    = total,
    Positive_Edges = pos,
    Negative_Edges = neg,
    Positive_Rate  = round(pos/total, 4),
    AUROC          = round(auc_roc, 4),
    AUPR           = round(auc_pr, 4),
    stringsAsFactors = FALSE
  )
}

# 6.3 Define BED paths for three transcription factors
tf_list <- list(
  GATA1  = "C:/Users/29321/Desktop/GSM970257_GATA1-F_peaks.bed",
  TAL1   = "C:/Users/29321/Desktop/GSM1816083_TAL1-F5.peak.bed",
  BCL11A = "C:/Users/29321/Desktop/GSE103445_HUDEP2_BCL11A_optimal_idr_peaks.bed"
)

# 6.4 Vector of four subset filenames
subset_files <- c(
  Sem_Fgrn   = out_sem_fgrn,
  Sem_Dcema  = out_sem_dcema,
  Fgrn_Dcema = out_fgrn_dcema,
  Three      = out_three
)

# 7. Validate each transcription factor and print results
for (tf_name in names(tf_list)) {
  bed_path <- tf_list[[tf_name]]
  cat("\n===== Validating transcription factor:", tf_name, " =====\n")
  
  # 7.1 Generate ChIP ground truth pair set for this TF
  chip_pairs <- get_chip_pairs(bed_path = bed_path, tf_name = tf_name)
  
  # 7.2 Evaluate each of the four subsets
  for (subset_name in names(subset_files)) {
    pred_file <- subset_files[[subset_name]]
    output_dir <- file.path("evaluation_results", paste0(subset_name, "_", tf_name))
    
    res <- evaluate_network(
      pred_file      = pred_file,
      chip_pairs_set = chip_pairs,
      out_dir        = output_dir
    )
    
    if (!is.null(res)) {
      if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      # Save summary table
      fwrite(res,
             file      = file.path(output_dir, paste0("summary_", subset_name, "_", tf_name, ".txt")),
             sep       = "\t",
             quote     = FALSE,
             col.names = TRUE)
      
      # Print to console
      cat("\n---", subset_name, "subset validation results (", tf_name, ")---\n")
      print(res)
    } else {
      cat("\n---", subset_name, "subset validation failed (", tf_name, "), possibly single positive/negative class ---\n")
    }
  }
}

# -----------------------------------------------------------------------------
# 8. Clean environment and free memory
# -----------------------------------------------------------------------------
rm(list = ls())
graphics.off()
gc()









# -----------------------------------------------------------------------------
# Batch ChIP-seq Validation for Transcription Factors
# -----------------------------------------------------------------------------
# Purpose:
# This script performs batch validation of predicted regulatory networks for 
# multiple transcription factors against ChIP-seq ground truth data. It evaluates
# prediction performance using AUROC Area Under ROC Curve and AUPR 
# Area Under Precision-Recall Curve metrics.
#
# Input:
#   1. Prediction files tab-separated containing predicted TF-target interactions
#      Each file must contain columns TF Target and EdgeWeight or Score
#      Currently configured for deepsem_predicted_all_edges_named_c1/c2/s1/s2.txt
#   2. ChIP-seq BED files one per transcription factor as ground truth
#      - GATA1  GSM970257_GATA1-F_peaks.bed
#      - TAL1   GSM1816083_TAL1-F5.peak.bed
#      - BCL11A GSE103445_HUDEP2_BCL11A_optimal_idr_peaks.bed
#
# Output:
#   1. For each transcription factor
#      - ROC curves PNG format for each prediction file
#      - PR curves PNG format for each prediction file
#      - Summary statistics file tab-separated with AUROC/AUPR values
#   2. Console output showing evaluation progress and completion status
#
# Note: All gene/peak names are cleaned by converting to lowercase and removing
#       non-alphanumeric characters to ensure consistent matching
# -----------------------------------------------------------------------------

# 1. Load required packages
library(data.table)
library(dplyr)
library(pROC)
library(PRROC)
library(ggplot2)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# 2. Clean and standardize names convert to lowercase remove non-alphanumeric
clean_name <- function(x) {
  x <- trimws(x)
  x <- tolower(x)
  gsub("[^a-z0-9]", "", x)
}

# 3. Annotate BED file to generate TF-Target pairs
#    Returns a character vector of TF_Target format for matching
get_chip_pairs <- function(bed_path, tf_name, tss_region = c(-2000, 2000)) {
  bed_raw <- fread(bed_path, header = FALSE)[, 1:3]
  colnames(bed_raw) <- c("chr", "start", "end")
  
  bed_clean <- bed_raw %>%
    filter(!is.na(chr), !is.na(start), !is.na(end)) %>%
    filter(suppressWarnings(is.numeric(start)), suppressWarnings(is.numeric(end))) %>%
    filter(start <= end)
  
  if (nrow(bed_clean) == 0) {
    stop("Cleaned BED file is empty ", bed_path)
  }
  
  tmp_bed <- tempfile(fileext = ".bed")
  fwrite(bed_clean, tmp_bed, sep = "\t", col.names = FALSE, quote = FALSE)
  
  peakAnno <- annotatePeak(
    tmp_bed,
    tssRegion = tss_region,
    TxDb      = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb    = "org.Hs.eg.db"
  )
  
  genes <- as.data.frame(peakAnno@anno)$SYMBOL
  chip_pairs <- paste0(clean_name(tf_name), "_", clean_name(genes))
  unique(chip_pairs)
}

# 4. Evaluate a single prediction file for AUROC/AUPR save plots and return stats
evaluate_network <- function(pred_file, chip_pairs_set, out_dir) {
  df <- fread(pred_file, header = TRUE, sep = "\t")
  # Must contain TF/Target columns
  if (!all(c("TF", "Target") %in% colnames(df))) {
    warning("Skipping missing TF/Target columns ", pred_file)
    return(NULL)
  }
  
  # Identify score column
  score_col <- if      ("EdgeWeight" %in% colnames(df)) "EdgeWeight"
  else if ("Score"      %in% colnames(df)) "Score"
  else {
    warning("Skipping missing Score/EdgeWeight column ", pred_file)
    return(NULL)
  }
  
  df$TF     <- clean_name(df$TF)
  df$Target <- clean_name(df$Target)
  pairs    <- paste0(df$TF, "_", df$Target)
  
  labels <- as.integer(pairs %in% chip_pairs_set)
  scores <- df[[score_col]]
  
  # Statistics
  total <- length(labels)
  pos   <- sum(labels)
  neg   <- total - pos
  
  # If only one class exists cannot compute curves
  if (length(unique(labels)) < 2) {
    auc_roc  <- NA
    auc_pr   <- NA
  } else {
    roc_obj  <- roc(labels, scores, quiet = TRUE)
    auc_roc  <- as.numeric(auc(roc_obj))
    pr_obj   <- pr.curve(
      scores.class0 = scores[labels == 1],
      scores.class1 = scores[labels == 0],
      curve = TRUE
    )
    auc_pr   <- pr_obj$auc.integral
    
    # Save ROC plot
    png(file.path(out_dir, paste0("ROC_", basename(pred_file), ".png")), 600, 600)
    plot(roc_obj, main = paste("ROC", basename(pred_file)), lwd = 2)
    dev.off()
    
    # Save PR plot
    png(file.path(out_dir, paste0("PR_", basename(pred_file), ".png")), 600, 600)
    plot(pr_obj$curve[,1], pr_obj$curve[,2],
         type = "l", xlab = "Recall", ylab = "Precision",
         main = paste("PR", basename(pred_file)), lwd = 2)
    dev.off()
  }
  
  data.frame(
    File           = basename(pred_file),
    Total_Edges    = total,
    Positive_Edges = pos,
    Negative_Edges = neg,
    Positive_Rate  = round(pos/total, 4),
    AUROC          = round(auc_roc, 4),
    AUPR           = round(auc_pr, 4),
    stringsAsFactors = FALSE
  )
}

# 5. Main workflow Define TFs and BED files prediction file list then execute
#    Modify the lists below as needed
tf_list <- list(
  list(name = "GATA1", bed = "C:/Users/29321/Desktop/GSM970257_GATA1-F_peaks.bed"),
  list(name = "TAL1",  bed = "C:/Users/29321/Desktop/GSM1816083_TAL1-F5.peak.bed"),
  list(name = "BCL11A",bed = "C:/Users/29321/Desktop/GSE103445_HUDEP2_BCL11A_optimal_idr_peaks.bed")
)

pred_files <- c(
  "deepsem_predicted_all_edges_named_c1.txt",
  "deepsem_predicted_all_edges_named_c2.txt",
  "deepsem_predicted_all_edges_named_s1.txt",
  "deepsem_predicted_all_edges_named_s2.txt"
)

# Evaluate each transcription factor
for (tf in tf_list) {
  cat("\n### Evaluating transcription factor", tf$name, " ###\n")
  chip_pairs <- get_chip_pairs(tf$bed, tf$name)
  
  out_dir <- file.path("evaluation_results", tf$name)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Summarize results from all networks
  results <- lapply(pred_files, evaluate_network, chip_pairs_set = chip_pairs, out_dir = out_dir)
  results <- do.call(rbind, results)
  
  # Save summary table
  fwrite(results,
         file = file.path(out_dir, paste0("summary_", tf$name, ".txt")),
         sep  = "\t")
  
  cat("Completed", tf$name, "evaluation results in", out_dir, "\n")
}

# Clear all objects from global environment
rm(list = ls())

# Close all graphics devices
graphics.off()

# Force garbage collection to free memory
gc()













