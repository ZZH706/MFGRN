# ============================================================================
# Poisson-based Key Transcription Factor Identification
# ============================================================================
# Purpose:
# This script identifies key transcription factors in regulatory networks using
# Poisson distribution modeling of out-degree centrality. It processes four
# network edge lists, filters TFs using a whitelist, calculates Poisson
# p-values and FDR-adjusted q-values, and outputs significant TFs for each network.
#
# Input:
#   1. Four edge list files (tab-separated) containing TF-Target interactions
#      - c1: Control pre-exercise (two_or_more_overlap_c1_TF_TARGET.txt)
#      - c2: Control post-exercise (two_or_more_overlap_c2_TF_TARGET.txt)
#      - s1: SCD pre-exercise (two_or_more_overlap_s1_TF_TARGET.txt)
#      - s2: SCD post-exercise (two_or_more_overlap_s2_TF_TARGET.txt)
#   2. TF whitelist file: Homo_TF_clean.txt containing known TF gene symbols
#
# Output:
#   1. Full results for each network (all TFs with degrees, p-values, q-values)
#      - *_poisson_all.txt
#   2. Significant TFs only for each network (FDR < 0.05 and degree >= threshold)
#      - *_poisson_TF.txt
#   3. Console summary showing lambda, degree threshold, and key TF count
#
# Notes:
#   - Poisson distribution parameter lambda = edges / nodes
#   - Degree threshold set at 99.9th percentile of Poisson distribution
#   - FDR correction using Benjamini-Hochberg method
#   - All gene names are matched against TF whitelist
# ============================================================================

library(data.table)
library(igraph)

# 1. Input file list with descriptive tags
edge_files <- c(
  c1 = "E:/1数据/≥2/evaluation_results-c1/two_or_more_overlap_c1_TF_TARGET.txt",
  c2 = "E:/1数据/≥2/evaluation_results-c2/two_or_more_overlap_c2_TF_TARGET.txt",
  s1 = "E:/1数据/≥2/evaluation_results-s1/two_or_more_overlap_s1_TF_TARGET.txt",
  s2 = "E:/1数据/≥2/evaluation_results-s2/two_or_more_overlap_s2_TF_TARGET.txt"
)

# 2. Descriptive labels for console output
label_map <- c(
  c1 = "Control Pre-exercise",
  c2 = "Control Post-exercise",
  s1 = "SCD Pre-exercise",
  s2 = "SCD Post-exercise"
)

# 3. Output directory
out_dir <- "E:/1数据/泊松分布筛选关键转录因子"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 4. Parameters
alpha     <- 0.05      # FDR threshold
tail_prob <- 0.999     # Degree threshold quantile

# 5. Read TF whitelist
tf_list <- fread("E:/1数据/Homo_TF_clean.txt",
                 header = TRUE, sep = "\t", colClasses = "character")[, TF]

# 6. Main processing loop
for (tag in names(edge_files)) {
  
  ef <- edge_files[[tag]]
  
  # 6.1 Read edge list
  edges <- fread(ef, header = TRUE, sep = "\t", colClasses = "character")
  setnames(edges, 1:2, c("TF", "Target"))
  
  # 6.2 Build graph
  g <- graph_from_data_frame(edges, directed = TRUE)
  
  # 6.3 Basic statistics
  n <- vcount(g)               # Number of nodes
  m <- ecount(g)               # Number of edges
  lambda <- m / n              # Expected out-degree (approx (n-1)*p)
  
  # 6.4 Calculate out-degree only for TF whitelist nodes
  deg_out <- degree(g, mode = "out")
  deg_tf  <- deg_out[names(deg_out) %in% tf_list]
  
  # 6.5 Poisson right-tail p-value and FDR
  p_raw <- 1 - ppois(deg_tf - 1, lambda = lambda)
  q_val <- p.adjust(p_raw, method = "BH")
  
  # 6.6 Determine significance
  deg_cut <- qpois(tail_prob, lambda = lambda)
  sig_idx <- which(q_val < alpha & deg_tf >= deg_cut)
  
  # 6.7 Save results
  net_tag <- sub("\\.txt$", "", basename(ef))
  
  res_all <- data.frame(
    Gene    = names(deg_tf),
    Degree  = as.integer(deg_tf),
    p_value = signif(p_raw,  3),
    q_value = signif(q_val,  3),
    stringsAsFactors = FALSE
  )
  res_sig <- res_all[sig_idx, ]
  
  fwrite(res_all,
         file = file.path(out_dir, paste0(net_tag, "_poisson_all.txt")),
         sep  = "\t", quote = FALSE)
  fwrite(res_sig,
         file = file.path(out_dir, paste0(net_tag, "_poisson_TF.txt")),
         sep  = "\t", quote = FALSE)
  
  # 6.8 Console summary
  cat(sprintf(
    "\n[%s] lambda = %.2f, degree threshold = %d, key TFs = %d\n",
    label_map[tag], lambda, deg_cut, nrow(res_sig)
  ))
}

cat("\nAll networks processed. Results saved to:", out_dir, "\n")