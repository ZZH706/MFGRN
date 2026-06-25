# ============================================================================
# Script: Infomap Module Detection and Key TF Enrichment Analysis
# ============================================================================
# Description:
#   This script identifies functional modules from gene regulatory networks 
#   using the Infomap algorithm, filters modules by size, and selects the 
#   top 10 modules enriched with key transcription factors (TFs) for each 
#   condition (c1, c2, s1, s2). Results are exported to an Excel workbook.
#
# Input data:
#   1. Four network edge lists (TF-Target interactions):
#      - E:/SCD/数据/合并后的数据/c1_TF_TARGET_cleaned.txt
#      - E:/SCD/数据/合并后的数据/c2_TF_TARGET_cleaned.txt
#      - E:/SCD/数据/合并后的数据/s1_TF_TARGET_cleaned.txt
#      - E:/SCD/数据/合并后的数据/s2_TF_TARGET_cleaned.txt
#      Format: at least two columns (Source, Target), optional third column (Weight)
#   2. Four key transcription factor lists (one per condition):
#      - E:/SCD/下游分析/关键转录因子/two_or_more_overlap_c1_TF_TARGET_poisson_TF_Gene_only.txt
#      - E:/SCD/下游分析/关键转录因子/two_or_more_overlap_c2_TF_TARGET_poisson_TF_Gene_only.txt
#      - E:/SCD/下游分析/关键转录因子/two_or_more_overlap_s1_TF_TARGET_poisson_TF_Gene_only.txt
#      - E:/SCD/下游分析/关键转录因子/two_or_more_overlap_s2_TF_TARGET_poisson_TF_Gene_only.txt
#      Format: one TF gene symbol per line
#
# Output data:
#   1. Rdata file: Infomap_modules_10_300_all.Rdata (all detected modules)
#   2. Excel file: Top10_modules_with_keyTF_infomap_10_300.xlsx
#      - Four condition sheets (c1, c2, s1, s2), each with three worksheets:
#        a) Summary: ModuleID, ModuleSize, KeyTF_Count, KeyTFs
#        b) Genes: ModuleID and all genes in each module
#        c) KeyTFs: ModuleID and key TF genes only
# ============================================================================

library(data.table)
library(openxlsx)
library(miRspongeR)

set.seed(123)

## =========================
## 0. Parameters and output directory
## =========================
min_size <- 10
max_size <- 300

out_dir <- "E:/SCD/下游分析/功能模块/模块"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## =========================
## 1. Read four networks
## =========================
net_files <- list(
  c1 = "E:/SCD/数据/合并后的数据/c1_TF_TARGET_cleaned.txt",
  c2 = "E:/SCD/数据/合并后的数据/c2_TF_TARGET_cleaned.txt",
  s1 = "E:/SCD/数据/合并后的数据/s1_TF_TARGET_cleaned.txt",
  s2 = "E:/SCD/数据/合并后的数据/s2_TF_TARGET_cleaned.txt"
)

# Compatible with different column names/numbers: ensure first two columns are edge ends; third column (if present) as weight
normalize_linklist <- function(dt) {
  if (!is.data.table(dt)) dt <- as.data.table(dt)
  if (ncol(dt) < 2) stop("Network edge file needs at least two columns (TF/Source, Target)")
  
  # Unify column names (without changing content)
  if (ncol(dt) == 2) {
    setnames(dt, names(dt)[1:2], c("Source", "Target"))
  } else {
    setnames(dt, names(dt)[1:3], c("Source", "Target", "Weight"))
  }
  dt
}

linklists <- lapply(net_files, function(fp) {
  if (!file.exists(fp)) stop("Network file does not exist: ", fp)
  normalize_linklist(fread(fp))
})

## =========================
## 2. Use infomap for module detection + filter module size 10~300
## =========================
run_infomap_and_filter <- function(linklist_dt, min_size = 10, max_size = 300) {
  # netModule's modulesize is "lower bound"; upper bound is filtered via post-processing
  clu <- netModule(linklist_dt, method = "infomap", modulesize = min_size)
  
  if (is.null(clu) || length(clu) == 0) return(clu)
  
  # Name modules if not already named
  if (is.null(names(clu)) || any(names(clu) == "")) {
    names(clu) <- paste0("M", seq_along(clu))
  }
  
  # Filter upper bound
  sizes <- sapply(clu, function(x) length(unique(x)))
  clu <- clu[sizes >= min_size & sizes <= max_size]
  
  # Renumber/re-name for clarity
  if (length(clu) > 0) names(clu) <- paste0("M", seq_along(clu))
  clu
}

clusters <- lapply(linklists, run_infomap_and_filter, min_size = min_size, max_size = max_size)

# Save module objects (optional)
save(clusters, file = file.path(out_dir, sprintf("Infomap_modules_%s_%s_all.Rdata", min_size, max_size)))

cat("\nInfomap module counts (after filtering 10~300):\n")
for (nm in names(clusters)) {
  cat(sprintf("%-3s : %d modules\n", nm, length(clusters[[nm]])))
}

## =========================
## 3. Read four key transcription factor lists
## =========================
tf_files <- list(
  c1 = "E:/SCD/下游分析/关键转录因子/two_or_more_overlap_c1_TF_TARGET_poisson_TF_Gene_only.txt",
  c2 = "E:/SCD/下游分析/关键转录因子/two_or_more_overlap_c2_TF_TARGET_poisson_TF_Gene_only.txt",
  s1 = "E:/SCD/下游分析/关键转录因子/two_or_more_overlap_s1_TF_TARGET_poisson_TF_Gene_only.txt",
  s2 = "E:/SCD/下游分析/关键转录因子/two_or_more_overlap_s2_TF_TARGET_poisson_TF_Gene_only.txt"
)

read_tf_list <- function(fp) {
  if (!file.exists(fp)) stop("TF list file does not exist: ", fp)
  x <- readLines(fp, warn = FALSE, encoding = "UTF-8")
  x <- trimws(x)
  x <- x[x != ""]
  unique(x)
}

tf_lists <- lapply(tf_files, read_tf_list)

cat("\nKey TF list sizes:\n")
for (nm in names(tf_lists)) {
  cat(sprintf("%-3s : %d TFs\n", nm, length(tf_lists[[nm]])))
}

## =========================
## 4. Select "top 10 modules containing key TFs"
##    Sorting: KeyTF_Count (descending) -> ModuleSize (descending) -> ModuleID (ascending)
## =========================
get_top10_modules_with_keyTF <- function(cluster_list, key_tfs, top_n = 10) {
  if (is.null(cluster_list) || length(cluster_list) == 0) {
    return(list(summary = data.frame(), module_genes = data.frame(), module_keytfs = data.frame()))
  }
  
  mod_ids <- names(cluster_list)
  if (is.null(mod_ids) || any(mod_ids == "")) {
    mod_ids <- paste0("M", seq_along(cluster_list))
    names(cluster_list) <- mod_ids
  }
  
  overlaps <- lapply(cluster_list, function(genes) intersect(unique(genes), key_tfs))
  keytf_n  <- sapply(overlaps, length)
  mod_size <- sapply(cluster_list, function(genes) length(unique(genes)))
  
  keep <- keytf_n > 0
  if (!any(keep)) {
    return(list(summary = data.frame(), module_genes = data.frame(), module_keytfs = data.frame()))
  }
  
  df <- data.frame(
    ModuleID    = mod_ids[keep],
    ModuleSize  = mod_size[keep],
    KeyTF_Count = keytf_n[keep],
    KeyTFs      = sapply(overlaps[keep], function(v) paste(v, collapse = ";")),
    stringsAsFactors = FALSE
  )
  
  df <- df[order(-df$KeyTF_Count, -df$ModuleSize, df$ModuleID), ]
  df_top <- head(df, top_n)
  
  module_genes <- do.call(rbind, lapply(df_top$ModuleID, function(mid) {
    data.frame(ModuleID = mid, Gene = unique(cluster_list[[mid]]), stringsAsFactors = FALSE)
  }))
  
  module_keytfs <- do.call(rbind, lapply(df_top$ModuleID, function(mid) {
    tfs <- overlaps[[mid]]
    if (length(tfs) == 0) return(NULL)
    data.frame(ModuleID = mid, KeyTF = tfs, stringsAsFactors = FALSE)
  }))
  
  list(summary = df_top, module_genes = module_genes, module_keytfs = module_keytfs)
}

## =========================
## 5. Process c1/c2/s1/s2 individually and export to Excel
## =========================
wb <- createWorkbook()

for (cond in c("c1", "c2", "s1", "s2")) {
  res <- get_top10_modules_with_keyTF(
    cluster_list = clusters[[cond]],
    key_tfs      = tf_lists[[cond]],
    top_n        = 10
  )
  
  sh1 <- paste0(cond, "_summary")
  sh2 <- paste0(cond, "_genes")
  sh3 <- paste0(cond, "_keyTFs")
  
  addWorksheet(wb, sh1)
  addWorksheet(wb, sh2)
  addWorksheet(wb, sh3)
  
  writeData(wb, sh1, res$summary)
  writeData(wb, sh2, res$module_genes)
  writeData(wb, sh3, res$module_keytfs)
  
  freezePane(wb, sh1, firstRow = TRUE)
  freezePane(wb, sh2, firstRow = TRUE)
  freezePane(wb, sh3, firstRow = TRUE)
  
  cat(sprintf("\n[%s] Modules containing key TFs (after filtering): %d; Exported Top10: %d\n",
              cond,
              ifelse(is.null(clusters[[cond]]), 0, length(clusters[[cond]])),
              nrow(res$summary)))
}

out_xlsx <- file.path(out_dir, sprintf("Top10_modules_with_keyTF_infomap_%s_%s.xlsx", min_size, max_size))
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("\nAll done! Results exported to:\n", out_xlsx, "\n")