###############################################################
#  Small-world Network Detection using Humphries sigma Index and Network Density
# ------------------------------------------------------------
#  Input  : TF Target EdgeWeight tab-separated .tsv file
#  Output : network_topology_sigma_<filename>.txt
#  sigma = (C_real / <C_rand>) / (L_real / <L_rand>), sigma > 1 indicates small-world
###############################################################

# -------------------- Required packages --------------------
library(igraph)
library(readr)
library(dplyr)
library(tools)

# -------------------- Parameters ----------------------
edge_file <- "filtered_networks/edges_s2_deepsem.tsv"   # Change this to use different file
score_thr <- 0          # EdgeWeight threshold
N_rand    <- 1000         # Number of random networks
set.seed(123)           # Reproducible results

# -------------------- Read data and build graph ----------------
cat("Reading:", edge_file, "\n")

edges <- read_tsv(edge_file,
                  col_types = cols(
                    TF         = col_character(),
                    Target     = col_character(),
                    EdgeWeight = col_double()
                  )) |>
  dplyr::filter(EdgeWeight > score_thr)

stopifnot(nrow(edges) > 0)

g_dir <- graph_from_data_frame(edges[, c("TF", "Target")],
                               directed = TRUE) |>
  igraph::simplify(remove.multiple = TRUE,
                   remove.loops     = TRUE)

g <- as.undirected(g_dir, mode = "collapse")

cat("Nodes:", vcount(g), "   Edges:", ecount(g), "\n")

# -------------------- Real network metrics ---------------
C_real <- transitivity(g, "average")
L_real <- average.path.length(g, unconnected = TRUE)
D_real <- edge_density(g, loops = FALSE)               # Network density

cat(sprintf("C = %.6f   L = %.6f   Density = %.6f\n",
            C_real, L_real, D_real))

# -------------------- Random network sampling ---------------
get_stats <- function(n_vert, n_edge) {
  gr <- sample_gnm(n_vert, n_edge, directed = FALSE, loops = FALSE)
  c(
    C = transitivity(gr, "average"),
    L = average.path.length(gr, unconnected = TRUE)
  )
}

cat("Generating", N_rand, "random networks...\n")
rand_stats <- replicate(N_rand, get_stats(vcount(g), ecount(g)))

C_rand_mu <- mean(rand_stats["C", ])
L_rand_mu <- mean(rand_stats["L", ])

cat(sprintf("<C_rand> = %.6f   <L_rand> = %.6f\n",
            C_rand_mu, L_rand_mu))

# -------------------- Sigma index and classification ---------------
sigma <- (C_real / C_rand_mu) / (L_real / L_rand_mu)
is_sw <- sigma > 1

cat(sprintf("sigma = %.3f  -> %s small-world network\n",
            sigma, ifelse(is_sw, "Yes", "No")))

# -------------------- Save results -------------------
out <- data.frame(
  Network          = file_path_sans_ext(basename(edge_file)),
  Nodes            = vcount(g),
  Edges            = ecount(g),
  Density          = round(D_real, 6),
  C_real           = round(C_real, 6),
  L_real           = round(L_real, 6),
  C_rand_mean      = round(C_rand_mu, 6),
  L_rand_mean      = round(L_rand_mu, 6),
  sigma            = round(sigma, 3),
  Small_world      = ifelse(is_sw, "Yes", "No"),
  Random_networks  = N_rand
)

out_file <- paste0("network_topology_sigma_",
                   file_path_sans_ext(basename(edge_file)), ".txt")

write.table(out,
            file = out_file,
            sep  = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

cat("\nTask complete, results saved to", out_file, "\n")









# ============================================================================
# Script: Radar Chart Visualization for Network Topological Properties
# ============================================================================
# Description:
#   This script generates radar charts to compare network topological 
#   properties (clustering coefficient, average path length, network density)
#   across different experimental conditions. Two versions are provided:
#   1. Four-panel comparison: predicted vs random networks for each condition
#   2. Three-panel comparison: predicted networks only across groups
#
# Input data:
#   Hard-coded network metrics for four conditions:
#   - ctrl_pre: Control group before exercise
#   - ctrl_post: Control group after exercise
#   - dis_pre: Disease group before exercise
#   - dis_post: Disease group after exercise
#   Metrics: clustering coefficient, average path length, network density
#
# Output data:
#   Radar charts rendered to the active graphics device (no file output)
#   - Version 1: 2x2 layout showing predicted vs random networks
#   - Version 2: 1x3 layout comparing predicted networks across groups
# ============================================================================

## Install and load package (install only once)
# install.packages("fmsb")
library(fmsb)

##-----------------------------
## 1. Enter all data
##-----------------------------
# Unify axis scales by calculating global max/min across all groups
all_C <- c(0.2751, 0.2295, 0.2938, 0.4227,  # Predicted: clustering coefficient
           0.0020, 0.0020, 0.0022, 0.0021) # Random: clustering coefficient
all_L <- c(2.3400, 2.3890, 2.1763, 2.2344, # Predicted: average path length
           2.99154, 2.9242, 2.8975, 2.9017)   # Random: average path length
all_D <- c(0.0020, 0.0020, 0.0022, 0.0021, # Predicted: network density
           0.0020, 0.0020, 0.0022, 0.0021) # Random: network density

# Unified max/min rows (first row = max, second row = min)
maxmin <- data.frame(
  ClusteringCoefficient = c(max(all_C) * 1.1, min(all_C) * 0.9),
  AveragePathLength = c(max(all_L) * 1.1, min(all_L) * 0.9),
  NetworkDensity = c(max(all_D) * 1.1, min(all_D) * 0.9)
)
rownames(maxmin) <- c("max", "min")

##-----------------------------
## 2. Construct data frames for each group
##-----------------------------
# Control group before exercise
dat_ctrl_pre <- rbind(
  maxmin,
  `Predicted_network` = c(0.2751, 2.3400, 0.0020),
  `Random_network` = c(0.0020, 2.99154, 0.0020)
)

# Control group after exercise
dat_ctrl_post <- rbind(
  maxmin,
  `Predicted_network` = c(0.2295, 2.3890, 0.0020),
  `Random_network` = c(0.0020, 2.9242, 0.0020)
)

# Disease group before exercise
dat_dis_pre <- rbind(
  maxmin,
  `Predicted_network` = c(0.2938, 2.1763, 0.0022),
  `Random_network` = c(0.0022, 2.8975, 0.0022)
)

# Disease group after exercise
dat_dis_post <- rbind(
  maxmin,
  `Predicted_network` = c(0.4227, 2.2344, 0.0021),
  `Random_network` = c(0.0021, 2.9017, 0.0021)
)

##-----------------------------
## 3. Draw 4 triangular radar charts (remove percentage scale labels)
##-----------------------------
par(mfrow = c(2, 2), mar = c(1, 2, 3, 2))  # 2x2 layout

radarchart(dat_ctrl_pre,
           axistype = 1,
           caxislabels = rep("", 5),    # Do not display percentage ticks
           pcol = c("red", "blue"),
           pfcol = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)),
           plwd = 2,
           title = "Control group before exercise")
legend("topright", legend = c("Predicted network", "Random network"),
       col = c("red","blue"), pch = 19, bty = "n")

radarchart(dat_ctrl_post,
           axistype = 1,
           caxislabels = rep("", 5),
           pcol = c("red", "blue"),
           pfcol = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)),
           plwd = 2,
           title = "Control group after exercise")
legend("topright", legend = c("Predicted network", "Random network"),
       col = c("red","blue"), pch = 19, bty = "n")

radarchart(dat_dis_pre,
           axistype = 1,
           caxislabels = rep("", 5),
           pcol = c("red", "blue"),
           pfcol = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)),
           plwd = 2,
           title = "Disease group before exercise")
legend("topright", legend = c("Predicted network", "Random network"),
       col = c("red","blue"), pch = 19, bty = "n")

radarchart(dat_dis_post,
           axistype = 1,
           caxislabels = rep("", 5),
           pcol = c("red", "blue"),
           pfcol = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)),
           plwd = 2,
           title = "Disease group after exercise")
legend("topright", legend = c("Predicted network", "Random network"),
       col = c("red","blue"), pch = 19, bty = "n")




## ============================================================================
# Script: Radar Chart Comparison for Predicted Networks Only
# ============================================================================
# Description:
#   This script generates radar charts comparing predicted network topological
#   properties (clustering coefficient, average path length, network density)
#   across three pairwise comparisons: control vs disease before exercise,
#   control vs disease after exercise, and disease before vs after exercise.
#
# Input data:
#   Hard-coded network metrics for predicted networks only:
#   - ctrl_pre: Control group before exercise (0.2751, 2.3400, 0.0020)
#   - ctrl_post: Control group after exercise (0.2295, 2.3890, 0.0020)
#   - dis_pre: Disease group before exercise (0.2938, 2.1763, 0.0022)
#   - dis_post: Disease group after exercise (0.4227, 2.2344, 0.0021)
#   Metrics: clustering coefficient, average path length, network density
#
# Output data:
#   Three radar charts rendered to the active graphics device (no file output)
#   in a 1x3 layout comparing:
#   1. Control before exercise vs Disease before exercise
#   2. Control after exercise vs Disease after exercise
#   3. Disease before exercise vs Disease after exercise
# ============================================================================

## Install and load package (install only once)
# install.packages("fmsb")
library(fmsb)

set.seed(123)

##-----------------------------
## 1. Enter four groups of "predicted network" metrics (no random networks)
##-----------------------------
# Three metrics: clustering coefficient / average path length / network density
ctrl_pre  <- c(0.2751, 2.3400, 0.0020)  # Control group before exercise
ctrl_post <- c(0.2295, 2.3890, 0.0020)  # Control group after exercise
dis_pre   <- c(0.2938, 2.1763, 0.0022)  # Disease group before exercise
dis_post  <- c(0.4227, 2.2344, 0.0021)  # Disease group after exercise

all_mat <- rbind(ctrl_pre, ctrl_post, dis_pre, dis_post)
colnames(all_mat) <- c("ClusteringCoefficient", "AveragePathLength", "NetworkDensity")

# Unified max/min rows (first row = max, second row = min)
maxmin <- data.frame(
  ClusteringCoefficient = c(max(all_mat[, "ClusteringCoefficient"]) * 1.1,     min(all_mat[, "ClusteringCoefficient"]) * 0.9),
  AveragePathLength = c(max(all_mat[, "AveragePathLength"]) * 1.1, min(all_mat[, "AveragePathLength"]) * 0.9),
  NetworkDensity = c(max(all_mat[, "NetworkDensity"]) * 1.1,     min(all_mat[, "NetworkDensity"]) * 0.9)
)
rownames(maxmin) <- c("max", "min")

##-----------------------------
## 2. Construct three comparison data frames (each chart: two lines A vs B)
##-----------------------------
dat_cmp1 <- rbind(
  maxmin,
  `Control_before_exercise` = ctrl_pre,
  `Disease_before_exercise` = dis_pre
)

dat_cmp2 <- rbind(
  maxmin,
  `Control_after_exercise` = ctrl_post,
  `Disease_after_exercise` = dis_post
)

dat_cmp3 <- rbind(
  maxmin,
  `Disease_before_exercise` = dis_pre,
  `Disease_after_exercise` = dis_post
)

##-----------------------------
## 3. Draw 3 radar charts (no scale labels)
##-----------------------------
par(mfrow = c(1, 3), mar = c(1, 2, 3, 2))  # 1x3 layout

# Unified plotting parameters
axis_n <- 5
cols_line <- c("red", "blue")
cols_fill <- c(rgb(1, 0, 0, 0.25), rgb(0, 0, 1, 0.25))

radarchart(dat_cmp1,
           axistype = 1,
           caxislabels = rep("", axis_n),  # Do not display scale labels
           pcol = cols_line,
           pfcol = cols_fill,
           plwd = 2,
           title = "Control before vs Disease before")
legend("topright", legend = c("Control before exercise", "Disease before exercise"),
       col = cols_line, pch = 19, bty = "n", cex = 0.9)

radarchart(dat_cmp2,
           axistype = 1,
           caxislabels = rep("", axis_n),
           pcol = cols_line,
           pfcol = cols_fill,
           plwd = 2,
           title = "Control after vs Disease after")
legend("topright", legend = c("Control after exercise", "Disease after exercise"),
       col = cols_line, pch = 19, bty = "n", cex = 0.9)

radarchart(dat_cmp3,
           axistype = 1,
           caxislabels = rep("", axis_n),
           pcol = cols_line,
           pfcol = cols_fill,
           plwd = 2,
           title = "Disease before vs Disease after")
legend("topright", legend = c("Disease before exercise", "Disease after exercise"),
       col = cols_line, pch = 19, bty = "n", cex = 0.9)
