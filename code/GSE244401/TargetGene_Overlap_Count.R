# =========================
# 1. 文件路径
# =========================
chip_files <- c(
  BCL11A = "E:/SCD/数据/CHIP-seq数据/BCL11A_ChIP_targets.tsv",
  GATA1  = "E:/SCD/数据/CHIP-seq数据/GATA1_ChIP_targets.tsv",
  TAL1   = "E:/SCD/数据/CHIP-seq数据/TAL1_ChIP_targets.tsv"
)

expr_files <- c(
  c1 = "E:/SCD/数据/构建网络的数据/data_c1_no0_transposed.csv",
  c2 = "E:/SCD/数据/构建网络的数据/data_c2_no0_transposed.csv",
  s1 = "E:/SCD/数据/构建网络的数据/data_s1_no0_transposed.csv",
  s2 = "E:/SCD/数据/构建网络的数据/data_s2_no0_transposed.csv"
)

# 输出文件夹
outdir <- "E:/SCD/数据/CHIP-seq数据/overlap_results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# =========================
# 2. 读取 ChIP-seq 靶基因
# =========================
read_chip_targets <- function(file) {
  df <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # 检查列名
  if (!all(c("TF", "Target") %in% colnames(df))) {
    stop(paste("文件缺少 TF 或 Target 列:", file))
  }
  
  genes <- unique(trimws(df$Target))
  genes <- genes[genes != "" & !is.na(genes)]
  return(genes)
}

chip_gene_list <- lapply(chip_files, read_chip_targets)

# =========================
# 3. 读取表达矩阵第一列基因名
# =========================
read_expr_genes <- function(file) {
  df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  
  # 第一列作为基因名
  genes <- unique(trimws(df[[1]]))
  genes <- genes[genes != "" & !is.na(genes)]
  return(genes)
}

expr_gene_list <- lapply(expr_files, read_expr_genes)

# =========================
# 4. 计算重叠并保存结果
# =========================
summary_list <- list()

for (expr_name in names(expr_gene_list)) {
  expr_genes <- expr_gene_list[[expr_name]]
  
  for (chip_name in names(chip_gene_list)) {
    chip_genes <- chip_gene_list[[chip_name]]
    
    overlap_genes <- intersect(expr_genes, chip_genes)
    overlap_genes <- sort(unique(overlap_genes))
    
    # 保存重叠基因列表
    overlap_df <- data.frame(Gene = overlap_genes, stringsAsFactors = FALSE)
    write.csv(
      overlap_df,
      file = file.path(outdir, paste0(expr_name, "_", chip_name, "_overlap_genes.csv")),
      row.names = FALSE
    )
    
    # 汇总统计
    summary_list[[paste(expr_name, chip_name, sep = "_")]] <- data.frame(
      Expression_File = expr_name,
      ChIP_File = chip_name,
      Expr_Gene_Count = length(expr_genes),
      ChIP_Target_Count = length(chip_genes),
      Overlap_Count = length(overlap_genes),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_list)

# 保存汇总表
write.csv(summary_df,
          file = file.path(outdir, "overlap_summary.csv"),
          row.names = FALSE)

# 在控制台打印
print(summary_df)








# =============================
# 1. 输入数据
# =============================
df <- data.frame(
  Group = c("CON_T1","CON_T1","CON_T1",
            "CON_T2","CON_T2","CON_T2",
            "SCD_T1","SCD_T1","SCD_T1",
            "SCD_T2","SCD_T2","SCD_T2"),
  TF = c("GATA1","TAL1","BCL11A",
         "GATA1","TAL1","BCL11A",
         "GATA1","TAL1","BCL11A",
         "GATA1","TAL1","BCL11A"),
  Left_only = c(605, 470, 661,
                978, 1196, 1868,
                983, 529, 244,
                1967, 311, 885),
  Overlap = c(162, 159, 267,
              249, 498, 896,
              274, 189, 169,
              469, 126, 552),
  Right_only = c(4016, 5627, 5808,
                 3915, 5248, 5190,
                 3773, 5390, 5845,
                 3609, 5478, 5474),
  stringsAsFactors = FALSE
)

# 四组背景数
bg_map <- c(
  CON_T1 = 25341,
  CON_T2 = 25288,
  SCD_T1 = 23587,
  SCD_T2 = 23932
)

# 给每行加背景数
df$Background <- bg_map[df$Group]

# =============================
# 2. Fisher精确检验（单尾，右尾）
# =============================
get_fisher_result <- function(left_only, overlap, right_only, background) {
  
  a <- overlap
  b <- left_only
  c <- right_only
  d <- background - a - b - c
  
  if (d < 0) {
    stop(paste("出现 d < 0，请检查输入数据。当前 d =", d))
  }
  
  # 2×2列联表
  mat <- matrix(c(a, b,
                  c, d),
                nrow = 2,
                byrow = TRUE)
  
  ft <- fisher.test(mat, alternative = "greater")
  
  # 显著性星号
  sig <- if (ft$p.value < 0.001) {
    "***"
  } else if (ft$p.value < 0.01) {
    "**"
  } else if (ft$p.value < 0.05) {
    "*"
  } else {
    "ns"
  }
  
  return(data.frame(
    a = a,
    b = b,
    c = c,
    d = d,
    p_value = ft$p.value,
    odds_ratio = unname(ft$estimate),
    sig = sig,
    stringsAsFactors = FALSE
  ))
}

# 对每一行计算
result_list <- lapply(seq_len(nrow(df)), function(i) {
  res <- get_fisher_result(
    left_only = df$Left_only[i],
    overlap = df$Overlap[i],
    right_only = df$Right_only[i],
    background = df$Background[i]
  )
  cbind(df[i, ], res)
})

result_df <- do.call(rbind, result_list)

# =============================
# 3. 美化输出
# =============================
result_df$p_value_scientific <- format(result_df$p_value, scientific = TRUE, digits = 4)
result_df$odds_ratio_round <- round(result_df$odds_ratio, 4)

# 调整列顺序
result_df <- result_df[, c(
  "Group", "TF", "Background",
  "a", "b", "c", "d",
  "p_value", "p_value_scientific",
  "odds_ratio", "odds_ratio_round",
  "sig"
)]

print(result_df)

# =============================
# 4. 保存结果
# =============================
write.csv(result_df,
          file = "E:/SCD/数据/fisher_one_tailed_results.csv",
          row.names = FALSE)