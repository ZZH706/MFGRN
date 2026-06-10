# -----------------------------------------------------------------------------
# 完整脚本：提取 ≥2 方法预测的边，并用 ChIP-seq 计算 AUROC/AUPR
# -----------------------------------------------------------------------------

# —— 0. 清理环境 & 加载包 —— 
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

# —— 1. 定义清洗函数 —— 
clean_name <- function(x) {
  x %>%
    trimws() %>%
    tolower() %>%
    gsub("[^a-z0-9]", "", .)
}

# —— 2. 读取三个方法预测文件 & 构建 dt_list —— 
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
  # 把 EdgeWeight 改名成方法专属列
  ew_col <- paste0("EW_", meth)
  setnames(dt, "EdgeWeight", ew_col)
  # 只保留需要的列（使用 .. 安全选择动态列名）
  keep_cols <- c("pair", "TF", "Target", ew_col)
  dt <- dt[, ..keep_cols]
  # 若同一方法内同一 pair 多行，取均值聚合（保持一行）
  dt <- dt[, .(
    TF = TF[1L],
    Target = Target[1L],
    tmp = mean(get(ew_col), na.rm = TRUE)
  ), by = pair]
  setnames(dt, "tmp", ew_col)
  dt
})
names(dt_list) <- names(paths)

# —— 3. 构建全集 pairs_map —— 
pairs_map <- rbindlist(
  lapply(dt_list, function(dt) dt[, .(pair, TF, Target)]),
  use.names = TRUE
)
pairs_map <- unique(pairs_map, by = "pair")

# —— 4. 合并各方法权重 —— 
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

# —— 5. 筛选出现 ≥2 方法的边 & 计算平均权重 —— 
ew_cols <- paste0("EW_", names(paths))
merged_dt[, count_present := rowSums(!is.na(.SD)), .SDcols = ew_cols]
two_or_more <- merged_dt[count_present >= 2]
two_or_more[, EdgeWeight := rowMeans(.SD, na.rm = TRUE), .SDcols = ew_cols]
two_or_more <- two_or_more[, .(TF, Target, EdgeWeight)]

# —— 6. 写出子集 —— 
out_pred <- "E:/1数据/filtered_networks/two_or_more_overlap_c1.txt"
fwrite(two_or_more, out_pred, sep="\t", quote=FALSE)

# —— 7. 定义 ChIP 注释函数 —— 
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

# —— 8. 定义评估函数：AUROC/AUPR & 绘图 —— 
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
  # 输出图
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

# —— 9. 批量 ChIP-seq 验证 —— 
tf_beds <- list(
  GATA1  = "D:/Users/29321/Desktop/GSM970257_GATA1-F_peaks.bed",
  NFE2   = "D:/Users/29321/Desktop/GSM1816086_NFE2-F5.peak.bed",
  TAL1   = "D:/Users/29321/Desktop/GSM1816083_TAL1-F5.peak.bed",
  BCL11A = "D:/Users/29321/Desktop/GSE103445_HUDEP2_BCL11A_optimal_idr_peaks.bed"
)

results <- list()
for (tf in names(tf_beds)) {
  cat(">> 验证 TF:", tf, "\n")
  chip_pairs <- get_chip_pairs(tf_beds[[tf]], tf)
  res <- evaluate_network(two_or_more, chip_pairs,
                          out_dir = file.path("E:/1数据/evaluation", tf))
  if (is.null(res)) {
    cat("  ⚠️ 单一标签，无法计算 ROC/PR。\n")
  } else {
    results[[tf]] <- cbind(TF = tf, res)
    print(results[[tf]])
    write.table(results[[tf]],
                file   = file.path("E:/1数据/evaluation", tf, "summary.txt"),
                sep    = "\t", quote = FALSE, row.names = FALSE)
  }
}

# —— 10. 汇总所有结果 —— 
summary_all <- rbindlist(lapply(results, as.data.frame), fill=TRUE)
fwrite(summary_all,
       "E:/1数据/evaluation/auc_ap_summary_all.txt",
       sep="\t", quote=FALSE)
cat("\n✔️ 全部完成，结果已保存在 E:/1数据/evaluation 下。\n")
