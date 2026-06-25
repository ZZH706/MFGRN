# ============================================================================
# TF Degree Centrality Analysis for Regulatory Networks
# ============================================================================
# Purpose:
# This script calculates degree centrality for transcription factors only across
# multiple regulatory networks. It filters TFs using a whitelist, computes
# degree centrality for each network, generates distribution plots, and creates
# heatmaps for visualization across networks.
#
# Input:
#   1. TF whitelist file: Homo_TF_clean.txt containing TF gene symbols
#   2. Network edge files: Tab-separated files with TF, Target, Weight columns
#      - c1_TF_TARGET_cleaned.txt (Control_T1)
#      - c2_TF_TARGET_cleaned.txt (Control_T2)
#      - s1_TF_TARGET_cleaned.txt (SCD_T1)
#      - s2_TF_TARGET_cleaned.txt (SCD_T2)
#
# Output:
#   1. Per-network degree centrality tables tab-separated
#      - degree_centrality_Control_T1.txt
#      - degree_centrality_Control_T2.txt
#      - degree_centrality_SCD_T1.txt
#      - degree_centrality_SCD_T2.txt
#   2. Combined degree data in R environment all_degree_data
#   3. Distribution plots: degree_centrality_TF_only.pdf
#   4. Heatmap with row names: TF_degree_heatmap_log10.pdf
#   5. Heatmap without row names: TF_degree_heatmap_log10_no_rownames.pdf
#
# Note: All gene names are trimmed and matched against TF whitelist
# ============================================================================

library(igraph)
library(dplyr)
library(tidyr)
library(pheatmap)

# 1. TF whitelist
tf_list <- read.table("E:/1数据/Homo_TF_clean.txt",
                      header = FALSE,
                      stringsAsFactors = FALSE)
all_tfs <- unique(trimws(tf_list[[1]]))

# 2. Network files
network_files <- c(
  "E:/1数据/≥2/c1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/c2_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s2_TF_TARGET_cleaned.txt"
)
network_names <- c("Control_T1", "Control_T2", "SCD_T1", "SCD_T2")

degree_centrality_list <- list()

# 3. Main calculation loop
for (i in seq_along(network_files)) {
  
  # Read edge table
  edge_data <- read.table(
    network_files[i],
    header = TRUE,
    sep    = "\t",
    stringsAsFactors = FALSE
  )
  
  # Build igraph object
  g <- graph_from_data_frame(edge_data[, 1:2], directed = TRUE)
  
  # Filter TFs using whitelist
  tf_in_graph <- intersect(all_tfs, V(g)$name)
  
  degree_values <- degree(g, mode = "all")[tf_in_graph]
  degree_values <- degree_values[!is.na(degree_values)]
  
  degree_df <- data.frame(
    Gene    = names(degree_values),
    Degree  = as.numeric(degree_values),
    Network = network_names[i],
    stringsAsFactors = FALSE
  )
  
  degree_centrality_list[[i]] <- degree_df
  
  # Save per-network results sorted by degree
  write.table(
    degree_df[order(-degree_df$Degree), ],
    file = paste0("degree_centrality_", network_names[i], ".txt"),
    sep  = "\t", quote = FALSE, row.names = FALSE
  )
}

# 4. Combine all network results
all_degree_data <- bind_rows(degree_centrality_list)

# 5. Generate distribution plots
pdf("degree_centrality_TF_only.pdf", width = 10, height = 8)
par(mfrow = c(2, 2))

for (nw in network_names) {
  
  sub <- all_degree_data %>% 
    filter(Network == nw) %>% 
    arrange(desc(Degree))
  
  plot(
    sub$Degree,
    type = "h", col = "black",
    main = paste("TF Degree -", nw),
    xlab = "TFs (sorted)",
    ylab = "Degree",
    lwd  = 1
  )
}
dev.off()

cat("\nTF degree calculation complete\n",
    "- Per-network results: degree_centrality_<Network>.txt\n",
    "- Combined data: all_degree_data in R environment\n",
    "- Plot file: degree_centrality_TF_only.pdf\n")

# 6. Convert to wide format matrix
degree_wide <- pivot_wider(all_degree_data,
                           names_from = Network,
                           values_from = Degree,
                           values_fill = 0)

degree_matrix <- as.matrix(degree_wide[,-1])
rownames(degree_matrix) <- degree_wide$Gene

# 7. Log10 transformation for heatmap
degree_log_matrix <- log10(degree_matrix + 1)

# 8. Generate heatmap with row names
pdf("TF_degree_heatmap_log10.pdf", width = 8, height = 10)
pheatmap(degree_log_matrix,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 6)
dev.off()

# 9. Generate heatmap without row names
pdf("TF_degree_heatmap_log10_no_rownames.pdf", width = 8, height = 10)
pheatmap(degree_log_matrix,
         scale = "none",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 6,
         show_rownames = FALSE)
dev.off()