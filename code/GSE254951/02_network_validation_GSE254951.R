# ==============================================================================
# Script: Comprehensive AUROC/AUPR Evaluation for GSE254951 Dataset
# 
# Description:
# This script evaluates the performance of gene regulatory network predictions
# for the GSE254951 dataset. It constructs and evaluates 8 different network 
# combinations:
#   1. Individual methods (3): deepfgrn, 3dcema, deepsem
#   2. Pairwise intersections (3): all combinations of two methods
#   3. Three-way intersection (1): common edges from all three methods
#   4. Union of two or more methods (1): edges supported by at least two methods
# 
# Total: 8 networks evaluated
#
# Input data:
#   - Network files (top 10% confidence edges):
#     * GSE254951_deepfgrn_top10pct.tsv
#     * GSE254951_3dcema_top10pct.tsv
#     * GSE254951_deepsem_top10pct.tsv
#   - Gold standard ChIP data: Integrated_ChIP_Top200.txt
#
# Output data (saved in AUROC directory):
#   - ALL_NETWORKS_SUMMARY.tsv: Summary table with AUROC and AUPR for all networks
#   - Individual network subdirectories with:
#     * network_edges.tsv: Edge list for that network
#     * evaluation_results.tsv: Detailed metrics
#     * ROC_curve.png: ROC curve visualization
#     * PR_curve.png: Precision-Recall curve visualization
#   - Visualizations:
#     * performance_comparison.png: Bar chart comparing all networks
#     * aurocs_by_type.png: AUROC distribution by network type
#     * aurocs_vs_network_size.png: Relationship between network size and performance
#     * method_overlap_matrix.tsv: Jaccard overlap between methods
#   - README.txt: Summary report
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Clean environment and load packages
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# 1. Define cleaning function
# ------------------------------------------------------------------------------
clean_name <- function(x) {
  x %>% trimws() %>% tolower() %>% gsub("[^a-z0-9]", "", .)
}

# ------------------------------------------------------------------------------
# 2. Input: network files (3 top10pct.tsv files)
# ------------------------------------------------------------------------------
net_paths <- list(
  deepfgrn = "E:/SCD/其他数据/GSE254951/GSE254951_deepfgrn_top10pct.tsv",
  `3dcema` = "E:/SCD/其他数据/GSE254951/GSE254951_3dcema_top10pct.tsv",
  deepsem  = "E:/SCD/其他数据/GSE254951/GSE254951_deepsem_top10pct.tsv"
)

# Output root directory
out_root <- "E:/SCD/其他数据/GSE254951/AUROC"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 3. Read three method prediction files and build dt_list
# ------------------------------------------------------------------------------
dt_list <- lapply(names(net_paths), function(meth) {
  dt <- fread(net_paths[[meth]], sep = "\t", header = TRUE)
  
  # Handle different weight column names
  if (!("EdgeWeight" %in% names(dt))) {
    cand <- intersect(names(dt), c("weight", "Weight", "score", "Score", "edgeweight", "Edge_Weight"))
    if (length(cand) >= 1) setnames(dt, cand[1], "EdgeWeight")
  }
  
  dt[, TF     := clean_name(TF)]
  dt[, Target := clean_name(Target)]
  dt[, pair   := paste0(TF, "_", Target)]
  
  # If duplicate pairs exist within the same method, aggregate by mean
  dt <- dt[, .(
    EdgeWeight = mean(EdgeWeight, na.rm = TRUE)
  ), by = .(TF, Target, pair)]
  
  # Add method name column
  dt[, Method := meth]
  dt
})
names(dt_list) <- names(net_paths)

# ------------------------------------------------------------------------------
# 4. Read gold standard data
# ------------------------------------------------------------------------------
chip_file <- "E:/SCD/其他数据/gse133181_单细胞/GSE133181_RAW/chip/Integrated_ChIP_Top200.txt"

cat("Reading gold standard file:", chip_file, "\n")
chip_dt <- fread(chip_file, sep = "\t", header = TRUE)

# Check column names
if (!all(c("TF", "Target") %in% names(chip_dt))) {
  stop("ChIP file must contain 'TF' and 'Target' columns")
}

# Clean and generate pairs
chip_dt[, TF_clean := clean_name(TF)]
chip_dt[, Target_clean := clean_name(Target)]
chip_dt[, pair := paste0(TF_clean, "_", Target_clean)]

# Get all gold standard pairs
all_chip_pairs <- unique(chip_dt$pair)

cat("Total gold standard edges:", nrow(chip_dt), "\n")
cat("Unique gold standard edges:", length(all_chip_pairs), "\n")
cat("Number of TFs:", length(unique(chip_dt$TF_clean)), "\n\n")

# ------------------------------------------------------------------------------
# 5. Define evaluation function
# ------------------------------------------------------------------------------
evaluate_network <- function(network_dt, network_name, gold_pairs, out_dir) {
  cat("\nEvaluating network:", network_name, "\n")
  cat("  Number of edges:", nrow(network_dt), "\n")
  
  # Prepare data
  df <- copy(network_dt)[
    , `:=`(
      label = as.integer(pair %in% gold_pairs),
      score = EdgeWeight
    )
  ]
  
  # Statistics
  pos_count <- sum(df$label)
  neg_count <- sum(1 - df$label)
  cat("  Positive samples:", pos_count, "(", round(pos_count/nrow(df)*100, 2), "%)\n")
  cat("  Negative samples:", neg_count, "(", round(neg_count/nrow(df)*100, 2), "%)\n")
  
  # Check if evaluation is possible
  if (length(unique(df$label)) < 2) {
    cat("  Warning: Single label, cannot compute ROC/PR\n")
    return(NULL)
  }
  
  # Compute AUROC
  roc_obj <- roc(df$label, df$score, quiet = TRUE, direction = "<")
  auroc <- auc(roc_obj)
  
  # Compute AUPR
  pr_obj <- tryCatch({
    pr.curve(
      scores.class0 = df$score[df$label == 1],
      scores.class1 = df$score[df$label == 0],
      curve = TRUE
    )
  }, error = function(e) {
    cat("  Warning: PR calculation error, using alternative method\n")
    return(NULL)
  })
  
  if (!is.null(pr_obj) && is.list(pr_obj) && !is.null(pr_obj$auc.integral)) {
    aupr <- pr_obj$auc.integral
  } else {
    # Approximate AUPR using pROC
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
  
  # Create network-specific directory
  net_dir <- file.path(out_dir, gsub("[^a-zA-Z0-9]", "_", network_name))
  dir.create(net_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Plot ROC curve
  png(file.path(net_dir, "ROC_curve.png"), 800, 600)
  plot(roc_obj, 
       main = sprintf("%s\nAUROC = %.4f", network_name, auroc),
       col = "blue", lwd = 2, legacy.axes = TRUE)
  abline(a = 0, b = 1, lty = 2, col = "gray")
  legend("bottomright", legend = c("ROC curve", "Random classifier"), 
         col = c("blue", "gray"), lty = c(1, 2), lwd = c(2, 1))
  dev.off()
  
  # Plot PR curve if successful
  if (!is.null(pr_obj) && is.list(pr_obj) && !is.null(pr_obj$curve)) {
    png(file.path(net_dir, "PR_curve.png"), 800, 600)
    plot(pr_obj$curve[, 1], pr_obj$curve[, 2],
         type = "l", lwd = 2, col = "red",
         xlab = "Recall", ylab = "Precision",
         main = sprintf("%s\nAUPR = %.4f", network_name, aupr))
    grid()
    dev.off()
  }
  
  # Save results
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
  
  # Save network edges
  fwrite(network_dt, file.path(net_dir, "network_edges.tsv"), sep = "\t", quote = FALSE)
  
  # Save results
  fwrite(result, file.path(net_dir, "evaluation_results.tsv"), sep = "\t", quote = FALSE)
  
  cat("  AUROC =", round(auroc, 4), ", AUPR =", round(aupr, 4), "\n")
  
  return(result)
}

# ------------------------------------------------------------------------------
# 6. Build all networks
# ------------------------------------------------------------------------------

cat("\n", rep("=", 80), "\n", sep = "")
cat("Building all network combinations\n")
cat(rep("=", 80), "\n", sep = "")

all_networks <- list()
method_names <- names(net_paths)

# 6.1 Individual method networks
cat("\n1. Building individual method networks...\n")
for (meth in method_names) {
  network_dt <- dt_list[[meth]][, .(TF, Target, pair, EdgeWeight)]
  all_networks[[meth]] <- network_dt
  cat("  ", meth, ":", nrow(network_dt), "edges\n")
}

# 6.2 Pairwise intersections
cat("\n2. Building pairwise intersection networks...\n")
combinations <- combn(method_names, 2, simplify = FALSE)
for (comb in combinations) {
  comb_name <- paste(comb, collapse = "_vs_")
  
  # Get pair sets from both methods
  pairs1 <- dt_list[[comb[1]]]$pair
  pairs2 <- dt_list[[comb[2]]]$pair
  
  # Calculate intersection
  intersect_pairs <- intersect(pairs1, pairs2)
  
  if (length(intersect_pairs) > 0) {
    # Extract intersection edges, weight = mean of both methods
    dt1_sub <- dt_list[[comb[1]]][pair %in% intersect_pairs, .(pair, TF, Target, Weight1 = EdgeWeight)]
    dt2_sub <- dt_list[[comb[2]]][pair %in% intersect_pairs, .(pair, Weight2 = EdgeWeight)]
    
    network_dt <- merge(dt1_sub, dt2_sub, by = "pair")
    network_dt[, EdgeWeight := (Weight1 + Weight2) / 2]
    network_dt <- network_dt[, .(TF, Target, pair, EdgeWeight)]
    
    all_networks[[comb_name]] <- network_dt
    cat("  ", comb_name, ":", nrow(network_dt), "edges\n")
  } else {
    cat("  ", comb_name, ": empty intersection\n")
  }
}

# 6.3 Three-way intersection
cat("\n3. Building three-way intersection network...\n")
pairs_all <- lapply(dt_list, function(x) x$pair)
common_pairs <- Reduce(intersect, pairs_all)

if (length(common_pairs) > 0) {
  # Extract weights from all three methods and take mean
  network_dt <- data.table(pair = common_pairs)
  for (meth in method_names) {
    dt_sub <- dt_list[[meth]][pair %in% common_pairs, .(pair, Weight = EdgeWeight)]
    setnames(dt_sub, "Weight", paste0("Weight_", meth))
    network_dt <- merge(network_dt, dt_sub, by = "pair", all.x = TRUE)
  }
  
  # Calculate mean weight
  weight_cols <- paste0("Weight_", method_names)
  network_dt[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = weight_cols]
  
  # Add TF and Target information
  network_dt <- merge(network_dt, dt_list[[1]][, .(pair, TF, Target)], by = "pair")
  network_dt <- network_dt[, .(TF, Target, pair, EdgeWeight)]
  
  all_networks[["AllThree_Intersection"]] <- network_dt
  cat("  AllThree_Intersection:", nrow(network_dt), "edges\n")
} else {
  cat("  AllThree_Intersection: empty intersection\n")
}

# 6.4 Union of two or more methods (edges supported by at least two methods)
cat("\n4. Building union of two or more methods network...\n")

# Merge all method weights
merged_dt <- data.table(pair = unique(unlist(lapply(dt_list, function(x) x$pair))))

for (meth in method_names) {
  dt_sub <- dt_list[[meth]][, .(pair, Weight = EdgeWeight)]
  setnames(dt_sub, "Weight", paste0("Weight_", meth))
  merged_dt <- merge(merged_dt, dt_sub, by = "pair", all.x = TRUE)
}

# Count how many methods support each pair
weight_cols <- paste0("Weight_", method_names)
merged_dt[, count_present := rowSums(!is.na(.SD)), .SDcols = weight_cols]

# Filter edges supported by >= 2 methods
two_or_more <- merged_dt[count_present >= 2]

if (nrow(two_or_more) > 0) {
  # Calculate mean weight
  two_or_more[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = weight_cols]
  
  # Add TF and Target information
  two_or_more <- merge(two_or_more, dt_list[[1]][, .(pair, TF, Target)], by = "pair")
  network_dt <- two_or_more[, .(TF, Target, pair, EdgeWeight)]
  
  all_networks[["TwoOrMore_Union"]] <- network_dt
  cat("  TwoOrMore_Union:", nrow(network_dt), "edges\n")
} else {
  cat("  TwoOrMore_Union: no edges\n")
}

# ------------------------------------------------------------------------------
# 7. Evaluate all networks
# ------------------------------------------------------------------------------
cat("\n", rep("=", 80), "\n", sep = "")
cat("Evaluating all networks\n")
cat(rep("=", 80), "\n", sep = "")

all_results <- list()

for (net_name in names(all_networks)) {
  result <- evaluate_network(all_networks[[net_name]], net_name, all_chip_pairs, out_root)
  if (!is.null(result)) {
    all_results[[net_name]] <- result
  }
}

# ------------------------------------------------------------------------------
# 8. Summarize all results
# ------------------------------------------------------------------------------
cat("\n", rep("=", 80), "\n", sep = "")
cat("Summary results\n")
cat(rep("=", 80), "\n", sep = "")

if (length(all_results) > 0) {
  summary_all <- rbindlist(all_results)
  
  # Sort by AUROC
  setorder(summary_all, -AUROC)
  
  # Save summary
  fwrite(summary_all, file.path(out_root, "ALL_NETWORKS_SUMMARY.tsv"), 
         sep = "\t", quote = FALSE)
  
  # Print summary table
  cat("\nNetwork performance summary:\n")
  print(summary_all)
  
  # ----------------------------------------------------------------------------
  # 9. Visualization
  # ----------------------------------------------------------------------------
  
  # 9.1 AUROC and AUPR bar chart
  summary_plot <- summary_all %>%
    pivot_longer(cols = c(AUROC, AUPR), names_to = "Metric", values_to = "Value")
  
  p1 <- ggplot(summary_plot, aes(x = reorder(Network, Value), y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "All Networks Performance Comparison",
         x = "Network", y = "Score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "performance_comparison.png"), p1, width = 12, height = 7)
  
  # 9.2 AUROC by type comparison
  p2 <- ggplot(summary_all, aes(x = Type, y = AUROC, fill = Type)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    theme_minimal() +
    labs(title = "AUROC by Network Type",
         x = "Network Type", y = "AUROC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "aurocs_by_type.png"), p2, width = 8, height = 6)
  
  # 9.3 AUROC vs network size
  p3 <- ggplot(summary_all, aes(x = N_Edges, y = AUROC, color = Type, label = Network)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
    theme_minimal() +
    labs(title = "AUROC vs Network Size",
         x = "Number of Edges", y = "AUROC") +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(file.path(out_root, "aurocs_vs_network_size.png"), p3, width = 10, height = 7)
  
  # 9.4 Method overlap matrix
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
  
  fwrite(as.data.table(overlap_matrix, keep.rownames = "Method"), 
         file.path(out_root, "method_overlap_matrix.tsv"), sep = "\t")
  
  # 9.5 Best network information
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Best network (by AUROC):\n")
  best_network <- summary_all[which.max(AUROC)]
  print(best_network)
  
  cat("\nBest network (by AUPR):\n")
  best_pr_network <- summary_all[which.max(AUPR)]
  print(best_pr_network)
  
  cat("\nNetwork size vs performance:\n")
  cat("  Smallest network:", summary_all[which.min(N_Edges), .(Network, N_Edges, AUROC)] %>% 
        mutate(across(where(is.numeric), ~round(., 4))), "\n")
  cat("  Largest network:", summary_all[which.max(N_Edges), .(Network, N_Edges, AUROC)] %>% 
        mutate(across(where(is.numeric), ~round(., 4))), "\n")
  
  cat("\nAll evaluations complete!\n")
  cat("Results saved in:", out_root, "\n")
  cat("  - Summary file: ALL_NETWORKS_SUMMARY.tsv\n")
  cat("  - Detailed results in network subdirectories\n")
  cat("  - Visualizations: performance_comparison.png, aurocs_by_type.png, aurocs_vs_network_size.png\n")
  
} else {
  cat("\nWarning: No networks were successfully evaluated.\n")
}

# ------------------------------------------------------------------------------
# 10. Generate report
# ------------------------------------------------------------------------------
report_file <- file.path(out_root, "README.txt")
writeLines(c(
  "========================================",
  "Network Performance Evaluation Report",
  "========================================",
  paste("Generated:", Sys.time()),
  "",
  "Evaluated Network Types:",
  "1. Individual method networks (3): deepfgrn, 3dcema, deepsem",
  "2. Pairwise intersections (3): deepfgrn_vs_3dcema, deepfgrn_vs_deepsem, 3dcema_vs_deepsem",
  "3. Three-way intersection (1): AllThree_Intersection",
  "4. Union of two or more methods (1): TwoOrMore_Union",
  "",
  "Output Files:",
  "- ALL_NETWORKS_SUMMARY.tsv: Summary of all networks",
  "- performance_comparison.png: Performance bar chart",
  "- aurocs_by_type.png: AUROC by network type",
  "- aurocs_vs_network_size.png: AUROC vs network size relationship",
  "- method_overlap_matrix.tsv: Method overlap matrix",
  "",
  "Each network has a dedicated subdirectory with:",
  "  - ROC and PR curves",
  "  - Network edge list",
  "  - Detailed evaluation results",
  "",
  "Gold Standard Source: Integrated_ChIP_Top200.txt",
  paste("Gold standard edges:", length(all_chip_pairs)),
  paste("Total predicted networks:", length(all_networks)),
  paste("Successfully evaluated:", length(all_results))
), con = report_file)

cat("\nReport saved to:", report_file, "\n")