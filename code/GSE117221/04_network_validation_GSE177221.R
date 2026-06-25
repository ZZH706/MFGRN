# ==============================================================================
# Script: Comprehensive AUROC/AUPR Evaluation for GSE117221 Dataset (H and TI)
# 
# Description:
# This script evaluates the performance of gene regulatory network predictions
# for both H and TI groups from the GSE117221 dataset. It constructs and evaluates
# 8 different network combinations for each group:
#   1. Individual methods (3): deepfgrn, 3dcema, deepsem
#   2. Pairwise intersections (3): all combinations of two methods
#   3. Three-way intersection (1): common edges from all three methods
#   4. Union of two or more methods (1): edges supported by at least two methods
# 
# Total: 8 networks per group, 16 networks across both groups
#
# Input data:
#   - Network files (top 10% confidence edges):
#     * deepfgrn_H_top10pct.tsv, deepfgrn_TI_top10pct.tsv
#     * 3dcema_H_top10pct.tsv, 3dcema_TI_top10pct.tsv
#     * deepsem_H_top10pct.tsv, deepsem_TI_top10pct.tsv
#   - Gold standard ChIP data: Integrated_ChIP_Top200.txt
#
# Output data (saved in separate directories for H and TI):
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
# Load required packages
# ------------------------------------------------------------------------------
rm(list = ls()); gc()

if (!requireNamespace("data.table", quietly=TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
if (!requireNamespace("pROC", quietly=TRUE)) install.packages("pROC")
if (!requireNamespace("PRROC", quietly=TRUE)) install.packages("PRROC")
if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly=TRUE)) install.packages("tidyr")

library(data.table)
library(dplyr)
library(pROC)
library(PRROC)
library(ggplot2)
library(tidyr)

# ------------------------------------------------------------------------------
# Define helper functions
# ------------------------------------------------------------------------------

# Clean gene names: trim whitespace, convert to lowercase, remove special characters
clean_name <- function(x) {
  x %>% trimws() %>% tolower() %>% gsub("[^a-z0-9]", "", .)
}

# Evaluate a single network: compute AUROC, AUPR, and generate plots
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

# Process a single group (H or TI)
process_group <- function(group_name, net_paths, gold_pairs, out_root) {
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Processing", group_name, "group\n")
  cat(rep("=", 80), "\n", sep = "")
  
  # Create output directory for this group
  out_dir <- file.path(out_root, paste0("AUROC_", group_name))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read network files
  dt_list <- lapply(names(net_paths), function(meth) {
    dt <- fread(net_paths[[meth]], sep = "\t", header = TRUE)
    
    # Handle different weight column names
    if (!("EdgeWeight" %in% names(dt))) {
      cand <- intersect(names(dt), c("weight", "Weight", "score", "Score", "edgeweight", "Edge_Weight"))
      if (length(cand) >= 1) setnames(dt, cand[1], "EdgeWeight")
    }
    
    dt[, TF := clean_name(TF)]
    dt[, Target := clean_name(Target)]
    dt[, pair := paste0(TF, "_", Target)]
    
    # Aggregate duplicate pairs by taking mean weight
    dt <- dt[, .(
      EdgeWeight = mean(EdgeWeight, na.rm = TRUE)
    ), by = .(TF, Target, pair)]
    
    dt[, Method := meth]
    dt
  })
  names(dt_list) <- names(net_paths)
  
  # Build all networks
  cat("\nBuilding all network combinations...\n")
  all_networks <- list()
  method_names <- names(net_paths)
  
  # 1. Individual methods
  cat("\n1. Building individual method networks...\n")
  for (meth in method_names) {
    network_dt <- dt_list[[meth]][, .(TF, Target, pair, EdgeWeight)]
    all_networks[[meth]] <- network_dt
    cat("  ", meth, ":", nrow(network_dt), "edges\n")
  }
  
  # 2. Pairwise intersections
  cat("\n2. Building pairwise intersection networks...\n")
  combinations <- combn(method_names, 2, simplify = FALSE)
  for (comb in combinations) {
    comb_name <- paste(comb, collapse = "_vs_")
    
    pairs1 <- dt_list[[comb[1]]]$pair
    pairs2 <- dt_list[[comb[2]]]$pair
    intersect_pairs <- intersect(pairs1, pairs2)
    
    if (length(intersect_pairs) > 0) {
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
  
  # 3. Three-way intersection
  cat("\n3. Building three-way intersection network...\n")
  pairs_all <- lapply(dt_list, function(x) x$pair)
  common_pairs <- Reduce(intersect, pairs_all)
  
  if (length(common_pairs) > 0) {
    network_dt <- data.table(pair = common_pairs)
    for (meth in method_names) {
      dt_sub <- dt_list[[meth]][pair %in% common_pairs, .(pair, Weight = EdgeWeight)]
      setnames(dt_sub, "Weight", paste0("Weight_", meth))
      network_dt <- merge(network_dt, dt_sub, by = "pair", all.x = TRUE)
    }
    
    weight_cols <- paste0("Weight_", method_names)
    network_dt[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = weight_cols]
    
    network_dt <- merge(network_dt, dt_list[[1]][, .(pair, TF, Target)], by = "pair")
    network_dt <- network_dt[, .(TF, Target, pair, EdgeWeight)]
    
    all_networks[["AllThree_Intersection"]] <- network_dt
    cat("  AllThree_Intersection:", nrow(network_dt), "edges\n")
  } else {
    cat("  AllThree_Intersection: empty intersection\n")
  }
  
  # 4. Union of two or more methods
  cat("\n4. Building union of two or more methods network...\n")
  
  merged_dt <- data.table(pair = unique(unlist(lapply(dt_list, function(x) x$pair))))
  
  for (meth in method_names) {
    dt_sub <- dt_list[[meth]][, .(pair, Weight = EdgeWeight)]
    setnames(dt_sub, "Weight", paste0("Weight_", meth))
    merged_dt <- merge(merged_dt, dt_sub, by = "pair", all.x = TRUE)
  }
  
  weight_cols <- paste0("Weight_", method_names)
  merged_dt[, count_present := rowSums(!is.na(.SD)), .SDcols = weight_cols]
  
  two_or_more <- merged_dt[count_present >= 2]
  
  if (nrow(two_or_more) > 0) {
    two_or_more[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = weight_cols]
    two_or_more <- merge(two_or_more, dt_list[[1]][, .(pair, TF, Target)], by = "pair")
    network_dt <- two_or_more[, .(TF, Target, pair, EdgeWeight)]
    
    all_networks[["TwoOrMore_Union"]] <- network_dt
    cat("  TwoOrMore_Union:", nrow(network_dt), "edges\n")
  } else {
    cat("  TwoOrMore_Union: no edges\n")
  }
  
  # Evaluate all networks
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Evaluating all networks\n")
  cat(rep("=", 80), "\n", sep = "")
  
  all_results <- list()
  
  for (net_name in names(all_networks)) {
    result <- evaluate_network(all_networks[[net_name]], net_name, gold_pairs, out_dir)
    if (!is.null(result)) {
      all_results[[net_name]] <- result
    }
  }
  
  # Summarize results
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Summary results\n")
  cat(rep("=", 80), "\n", sep = "")
  
  if (length(all_results) > 0) {
    summary_all <- rbindlist(all_results)
    setorder(summary_all, -AUROC)
    
    fwrite(summary_all, file.path(out_dir, "ALL_NETWORKS_SUMMARY.tsv"), 
           sep = "\t", quote = FALSE)
    
    cat("\nNetwork performance summary:\n")
    print(summary_all)
    
    # Generate visualizations
    # Bar chart comparing all networks
    summary_plot <- summary_all %>%
      pivot_longer(cols = c(AUROC, AUPR), names_to = "Metric", values_to = "Value")
    
    p1 <- ggplot(summary_plot, aes(x = reorder(Network, Value), y = Value, fill = Metric)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      labs(title = paste("GSE117221", group_name, "Dataset: All Networks Performance Comparison"),
           x = "Network", y = "Score") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      coord_cartesian(ylim = c(0, 1))
    
    ggsave(file.path(out_dir, "performance_comparison.png"), p1, width = 12, height = 7)
    
    # AUROC by network type
    p2 <- ggplot(summary_all, aes(x = Type, y = AUROC, fill = Type)) +
      geom_boxplot() +
      geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
      theme_minimal() +
      labs(title = paste("GSE117221", group_name, "Dataset: AUROC by Network Type"),
           x = "Network Type", y = "AUROC") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      coord_cartesian(ylim = c(0, 1))
    
    ggsave(file.path(out_dir, "aurocs_by_type.png"), p2, width = 8, height = 6)
    
    # AUROC vs network size
    p3 <- ggplot(summary_all, aes(x = N_Edges, y = AUROC, color = Type, label = Network)) +
      geom_point(size = 3) +
      geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
      theme_minimal() +
      labs(title = paste("GSE117221", group_name, "Dataset: AUROC vs Network Size"),
           x = "Number of Edges", y = "AUROC") +
      scale_x_log10() +
      coord_cartesian(ylim = c(0, 1))
    
    ggsave(file.path(out_dir, "aurocs_vs_network_size.png"), p3, width = 10, height = 7)
    
    # Method overlap matrix
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
           file.path(out_dir, "method_overlap_matrix.tsv"), sep = "\t")
    
    # Print best network information
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("Best network (by AUROC):\n")
    best_network <- summary_all[which.max(AUROC)]
    print(best_network)
    
    cat("\nBest network (by AUPR):\n")
    best_pr_network <- summary_all[which.max(AUPR)]
    print(best_pr_network)
    
    cat("\nNetwork size vs performance:\n")
    min_net <- summary_all[which.min(N_Edges), .(Network, N_Edges, AUROC)]
    max_net <- summary_all[which.max(N_Edges), .(Network, N_Edges, AUROC)]
    cat("  Smallest network:", 
        paste(round(min_net[, .(N_Edges, AUROC)], 4) %>% 
                mutate(Network = min_net$Network), collapse = " "), "\n")
    cat("  Largest network:", 
        paste(round(max_net[, .(N_Edges, AUROC)], 4) %>% 
                mutate(Network = max_net$Network), collapse = " "), "\n")
    
    cat("\nAll evaluations complete!\n")
    cat("Results saved in:", out_dir, "\n")
    cat("  - Summary file: ALL_NETWORKS_SUMMARY.tsv\n")
    cat("  - Detailed results in network subdirectories\n")
    cat("  - Visualizations: performance_comparison.png, aurocs_by_type.png, aurocs_vs_network_size.png\n")
    
    # Generate report
    report_file <- file.path(out_dir, "README.txt")
    writeLines(c(
      "========================================",
      paste("Network Performance Evaluation Report - GSE117221", group_name, "Dataset"),
      "========================================",
      paste("Generated:", Sys.time()),
      "",
      "Dataset Information:",
      paste("  - GSE117221 (", group_name, " series samples)", sep = ""),
      "  - Gene regulatory network predictions",
      "",
      "Evaluated Network Types:",
      "  1. Individual method networks (3): deepfgrn, 3dcema, deepsem",
      "  2. Pairwise intersections (3): all combinations",
      "  3. Three-way intersection (1): AllThree_Intersection",
      "  4. Union of two or more methods (1): TwoOrMore_Union",
      "",
      "Input Files:",
      paste("  - deepfgrn:", net_paths$deepfgrn),
      paste("  - 3dcema:", net_paths$`3dcema`),
      paste("  - deepsem:", net_paths$deepsem),
      "",
      "Gold Standard: Integrated_ChIP_Top200.txt",
      paste("Gold standard edges:", length(gold_pairs)),
      "",
      "Output Files:",
      "  - ALL_NETWORKS_SUMMARY.tsv: Summary of all networks",
      "  - performance_comparison.png: Performance bar chart",
      "  - aurocs_by_type.png: AUROC by network type",
      "  - aurocs_vs_network_size.png: AUROC vs network size relationship",
      "  - method_overlap_matrix.tsv: Method overlap matrix",
      "",
      "Each network has a dedicated subdirectory with:",
      "  - ROC and PR curves",
      "  - Network edge list",
      "  - Detailed evaluation results",
      "",
      paste("Total predicted networks:", length(all_networks)),
      paste("Successfully evaluated:", length(all_results))
    ), con = report_file)
    
    cat("\nReport saved to:", report_file, "\n")
    
  } else {
    cat("\nWarning: No networks were successfully evaluated.\n")
  }
  
  return(all_results)
}

# ------------------------------------------------------------------------------
# Main execution
# ------------------------------------------------------------------------------

# Set paths for H group
net_paths_H <- list(
  deepfgrn = "E:/SCD/其他数据/GSE117221/grn/deepfgrn_H_top10pct.tsv",
  `3dcema` = "E:/SCD/其他数据/GSE117221/grn/3dcema_H_top10pct.tsv",
  deepsem  = "E:/SCD/其他数据/GSE117221/grn/deepsem_H_top10pct.tsv"
)

# Set paths for TI group
net_paths_TI <- list(
  deepfgrn = "E:/SCD/其他数据/GSE117221/grn/deepfgrn_TI_top10pct.tsv",
  `3dcema` = "E:/SCD/其他数据/GSE117221/grn/3dcema_TI_top10pct.tsv",
  deepsem  = "E:/SCD/其他数据/GSE117221/grn/deepsem_TI_top10pct.tsv"
)

# Gold standard file
chip_file <- "E:/SCD/其他数据/gse133181_单细胞/GSE133181_RAW/chip/Integrated_ChIP_Top200.txt"

# Check if gold standard file exists
if (!file.exists(chip_file)) {
  stop(paste("Gold standard file not found:", chip_file, 
             "\nPlease update the correct file path!"))
}

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

# Output root directory
out_root <- "E:/SCD/其他数据/GSE117221"

# Process H group
cat("\n", rep("#", 80), "\n", sep = "")
cat("# PROCESSING H GROUP\n")
cat(rep("#", 80), "\n", sep = "")

results_H <- process_group("H", net_paths_H, all_chip_pairs, out_root)

# Process TI group
cat("\n", rep("#", 80), "\n", sep = "")
cat("# PROCESSING TI GROUP\n")
cat(rep("#", 80), "\n", sep = "")

results_TI <- process_group("TI", net_paths_TI, all_chip_pairs, out_root)

# Final summary
cat("\n", rep("=", 80), "\n", sep = "")
cat("ALL PROCESSING COMPLETE!\n")
cat(rep("=", 80), "\n", sep = "")

cat("\nSummary:\n")
cat("  - H group: Processed", length(results_H), "networks\n")
cat("  - TI group: Processed", length(results_TI), "networks\n")
cat("\nResults saved in:\n")
cat("  -", file.path(out_root, "AUROC_H"), "\n")
cat("  -", file.path(out_root, "AUROC_TI"), "\n")
cat("\n")