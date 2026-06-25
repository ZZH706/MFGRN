# ============================================================================
# TF Degree Centrality Comparison Analysis
# ============================================================================
# Purpose:
# This script performs comparative analysis of transcription factor degree
# centrality across different experimental conditions. It identifies TFs with
# significant changes in degree centrality between conditions, generates
# ranked lists of differentially connected TFs, and produces visualization plots.
#
# Input:
#   Tab-separated files containing TF degree centrality results with columns:
#   - Gene: Transcription factor gene symbol
#   - Degree: Degree centrality value
#   - p_value: Statistical significance p-value (optional)
#   - q_value: FDR-adjusted q-value (optional)
#   
#   Specific files analyzed:
#   1. Control_T1 (control before exercise)
#   2. Control_T2 (control after exercise)  
#   3. SCD_T1 (SCD before exercise)
#   4. SCD_T2 (SCD after exercise)
#
# Output:
#   For each comparison (3 comparisons total):
#   1. Full comparison table with all genes (TSV)
#   2. Top 100 increased degree TFs (TSV)
#   3. Top 100 decreased degree TFs (TSV)
#   4. TFs with large degree difference (|diff| >= 5) (TSV)
#   5. Scatter plot visualization (PNG)
#   
#   Comparisons performed:
#   - Control_T1 vs SCD_T1 (effect of disease at baseline)
#   - SCD_T1 vs SCD_T2 (effect of exercise in SCD group)
#   - Control_T2 vs SCD_T2 (effect of disease after exercise)
#
# Note: All gene names are trimmed and matched consistently
# ============================================================================

library(data.table)

# ============================================================================
# Function: Compare TF degree centrality between two groups
# ============================================================================
compare_TF_degree <- function(file_A, file_B, comparison_name) {
  
  #--------------------------------------------
  # 1. Read input files
  #--------------------------------------------
  A <- fread(file_A)
  B <- fread(file_B)
  
  # Display column names for verification
  print(colnames(A))
  print(colnames(B))
  
  #--------------------------------------------
  # 2. Merge data by gene
  #--------------------------------------------
  tf <- merge(
    A, B,
    by = "Gene",
    suffixes = c("_A", "_B"),
    all = TRUE
  )
  
  # Missing Degree set to 0 (only appears in one group)
  tf[is.na(Degree_A), Degree_A := 0]
  tf[is.na(Degree_B), Degree_B := 0]
  
  #--------------------------------------------
  # 3. Calculate difference metrics
  #--------------------------------------------
  tf[, Degree_diff     := Degree_B - Degree_A]
  tf[, Degree_abs_diff := abs(Degree_diff)]
  tf[, Degree_log2FC   := log2((Degree_B + 1) / (Degree_A + 1))]
  
  #--------------------------------------------
  # 4. Extract top 100 differential TFs by direction
  #--------------------------------------------
  Increase_top100 <- tf[Degree_diff > 0][order(-Degree_diff)][1:100]
  Decrease_top100 <- tf[Degree_diff < 0][order(Degree_diff)][1:100]
  
  # Large difference TFs (|diff| >= 5, threshold adjustable)
  tf_bigdiff <- tf[Degree_abs_diff >= 5]
  
  #--------------------------------------------
  # 5. Export results
  #--------------------------------------------
  fwrite(tf,
         paste0(comparison_name, "_compare_all.txt"),
         sep = "\t", quote = FALSE)
  
  fwrite(Increase_top100,
         paste0(comparison_name, "_top100_increase.txt"),
         sep = "\t", quote = FALSE)
  
  fwrite(Decrease_top100,
         paste0(comparison_name, "_top100_decrease.txt"),
         sep = "\t", quote = FALSE)
  
  fwrite(tf_bigdiff,
         paste0(comparison_name, "_bigdiff.txt"),
         sep = "\t", quote = FALSE)
  
  #--------------------------------------------
  # 6. Generate scatter plot
  #--------------------------------------------
  png(paste0(comparison_name, "_scatter.png"),
      width = 1200, height = 1200, res = 150)
  
  plot(tf$Degree_A,
       tf$Degree_B,
       xlab = "Degree (Group A)",
       ylab = "Degree (Group B)",
       main = paste0("TF Degree Comparison: ", comparison_name))
  abline(0, 1, lty = 2)
  
  dev.off()
  
  # Return results as list
  return(list(
    all = tf,
    increase = Increase_top100,
    decrease = Decrease_top100
  ))
}

# ============================================================================
# Execute three specified comparisons
# ============================================================================

# File paths
file_CT1 <- "E:/SCD/下游分析/度中心/degree_centrality_Control_T1.txt"
file_CT2 <- "E:/SCD/下游分析/度中心/degree_centrality_Control_T2.txt"
file_ST1 <- "E:/SCD/下游分析/度中心/degree_centrality_SCD_T1.txt"
file_ST2 <- "E:/SCD/下游分析/度中心/degree_centrality_SCD_T2.txt"

#-----------------------------------------------------------
# Comparison 1: Control_T1 vs SCD_T1 (Disease effect at baseline)
#-----------------------------------------------------------
compare_TF_degree(
  file_A = file_CT1,
  file_B = file_ST1,
  comparison_name = "Control_T1_vs_SCD_T1"
)

#-----------------------------------------------------------
# Comparison 2: SCD_T1 vs SCD_T2 (Exercise effect in SCD group)
#-----------------------------------------------------------
compare_TF_degree(
  file_A = file_ST1,
  file_B = file_ST2,
  comparison_name = "SCD_T1_vs_SCD_T2"
)

#-----------------------------------------------------------
# Comparison 3: Control_T2 vs SCD_T2 (Disease effect after exercise)
#-----------------------------------------------------------
compare_TF_degree(
  file_A = file_CT2,
  file_B = file_ST2,
  comparison_name = "Control_T2_vs_SCD_T2"
)