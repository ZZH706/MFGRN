library(readxl)
library(dplyr)
library(stringr)
library(data.table)

# ---------------------------
# 路径设置
# ---------------------------
expr_xlsx  <- "D:/放假/SCD/其他数据/GSE139912/GSE139912.xlsx"
prior_file <- "D:/放假/SCD/数据/构建网络的数据/label_c1.csv"
out_dir    <- "D:/放假/SCD/其他数据/GSE139912"

out_expr_top20k_csv   <- file.path(out_dir, "GSE139912_SCD_Baseline_Raw_Top20000.csv")
out_prior_top20k_csv  <- file.path(out_dir, "label_c1_prior_filtered_by_Top20000.csv")
out_expr_top20k_T_csv <- file.path(out_dir, "GSE139912_SCD_Baseline_Raw_Top20000_transposed.csv")

# ---------------------------
# 1) 读表达数据：只取 Raw_*_Baseline
# ---------------------------
df <- readxl::read_xlsx(expr_xlsx)

baseline_cols_raw <- names(df)[
  str_detect(names(df), "^Raw_") & str_detect(names(df), "_Baseline$")
]

expr <- df %>%
  select(ID, all_of(baseline_cols_raw)) %>%
  distinct(ID, .keep_all = TRUE)

# 转 numeric
expr[baseline_cols_raw] <- lapply(expr[baseline_cols_raw], function(x) as.numeric(x))

# ---------------------------
# 2) 删除所有 Baseline 样本中表达全为0的基因（NA当0）
# ---------------------------
mat <- as.matrix(expr[, baseline_cols_raw])
mat[is.na(mat)] <- 0

expr_nonzero <- expr[rowSums(mat) > 0, , drop = FALSE]

# ---------------------------
# 3) 选表达量 Top 20000 基因
#    用每个基因在所有 Baseline 样本 raw counts 的“总和”作为筛选指标
# ---------------------------
mat2 <- as.matrix(expr_nonzero[, baseline_cols_raw])
mat2[is.na(mat2)] <- 0
gene_sum <- rowSums(mat2)

expr_nonzero$TotalCount <- gene_sum

expr_top20k <- expr_nonzero %>%
  arrange(desc(TotalCount)) %>%
  slice_head(n = 20000) %>%
  select(-TotalCount)

# 保存 Top20000 表达矩阵（genes x samples）
write.csv(expr_top20k, out_expr_top20k_csv, row.names = FALSE, quote = TRUE)

gene_set_20k <- unique(expr_top20k$ID)

# ---------------------------
# 4) 读先验调控对文件（兼容 csv 或空白分隔）
# ---------------------------
read_prior_pairs <- function(path) {
  pri <- tryCatch(
    read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(pri) || ncol(pri) < 2) {
    pri <- read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  pri <- pri[, 1:2]
  colnames(pri) <- c("Gene1", "Gene2")
  pri <- pri %>%
    filter(!is.na(Gene1), !is.na(Gene2), Gene1 != "", Gene2 != "")
  pri
}

prior <- read_prior_pairs(prior_file)

# ---------------------------
# 5) 用 Top20000 基因集合筛先验调控对
# ---------------------------
prior_top20k <- prior %>%
  filter(Gene1 %in% gene_set_20k, Gene2 %in% gene_set_20k) %>%
  distinct()

write.csv(prior_top20k, out_prior_top20k_csv, row.names = FALSE, quote = TRUE)

# ---------------------------
# 6) 转置 Top20000 表达矩阵：samples x genes
# ---------------------------
# 用 ID 作为行名再转置
mat_top20k <- as.matrix(expr_top20k[, baseline_cols_raw])
rownames(mat_top20k) <- expr_top20k$ID
mat_top20k[is.na(mat_top20k)] <- 0

mat_top20k_t <- t(mat_top20k)

# 输出：第一列为 Sample
out_t <- data.frame(Sample = rownames(mat_top20k_t), mat_top20k_t, check.names = FALSE)
write.csv(out_t, out_expr_top20k_T_csv, row.names = FALSE, quote = TRUE)

# ---------------------------
# 7) 汇总信息
# ---------------------------
cat("Baseline样本数:", length(baseline_cols_raw), "\n")
cat("去掉全0基因后基因数:", nrow(expr_nonzero), "\n")
cat("Top20000 基因数:", nrow(expr_top20k), "\n")
cat("原始先验边数:", nrow(prior), "\n")
cat("Top20000筛选后先验边数:", nrow(prior_top20k), "\n")
cat("输出表达矩阵(Top20000):", out_expr_top20k_csv, "\n")
cat("输出先验边(Top20000):", out_prior_top20k_csv, "\n")
cat("输出转置表达矩阵:", out_expr_top20k_T_csv, "\n")













library(readxl)
library(dplyr)
library(stringr)
library(data.table)

# ---------------------------
# 路径设置
# ---------------------------
expr_xlsx  <- "D:/放假/SCD/其他数据/GSE139912/GSE139912.xlsx"
prior_file <- "D:/放假/SCD/数据/构建网络的数据/label_c1.csv"
tf_file    <- "D:/放假/SCD/数据/构建网络的数据/human_tf_list.txt"
out_dir    <- "D:/放假/SCD/其他数据/GSE139912"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 输出文件路径
out_expr_top10k_csv   <- file.path(out_dir, "GSE139912_SCD_Baseline_Raw_Top10000.csv")
out_prior_top10k_csv  <- file.path(out_dir, "label_c1_prior_filtered_by_Top10000.csv")
out_expr_top10k_T_csv <- file.path(out_dir, "GSE139912_SCD_Baseline_Raw_Top10000_transposed.csv")
out_tf_top10k_txt     <- file.path(out_dir, "GSE139912_TF_list_in_Top10000.txt")

# ---------------------------
# Step 1) 读取表达数据：只取 Raw_*_Baseline
# ---------------------------
df <- readxl::read_xlsx(expr_xlsx)

baseline_cols_raw <- names(df)[
  str_detect(names(df), "^Raw_") & str_detect(names(df), "_Baseline$")
]

expr <- df %>%
  dplyr::select(ID, all_of(baseline_cols_raw)) %>%
  dplyr::distinct(ID, .keep_all = TRUE)

# 转 numeric
expr[baseline_cols_raw] <- lapply(expr[baseline_cols_raw], function(x) as.numeric(x))

# ---------------------------
# Step 2) 删除所有 Baseline 样本中表达全为0的基因
# ---------------------------
mat <- as.matrix(expr[, baseline_cols_raw])
mat[is.na(mat)] <- 0

expr_nonzero <- expr[rowSums(mat) > 0, , drop = FALSE]

# ---------------------------
# Step 3) 选表达量 Top10000 基因（按总表达量）
# ---------------------------
mat2 <- as.matrix(expr_nonzero[, baseline_cols_raw])
mat2[is.na(mat2)] <- 0

expr_nonzero$TotalCount <- rowSums(mat2)

expr_top10k <- expr_nonzero %>%
  dplyr::arrange(desc(TotalCount)) %>%
  dplyr::slice_head(n = 10000) %>%
  dplyr::select(-TotalCount)

# 保存 Top10000 表达矩阵（genes × samples）
write.csv(expr_top10k, out_expr_top10k_csv, row.names = FALSE, quote = TRUE)

gene_set_10k <- unique(expr_top10k$ID)

# ---------------------------
# Step 4) 读取先验调控对文件
# ---------------------------
read_prior_pairs <- function(path) {
  pri <- tryCatch(
    read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(pri) || ncol(pri) < 2) {
    pri <- read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  pri <- pri[, 1:2]
  colnames(pri) <- c("Gene1", "Gene2")
  pri <- pri %>%
    dplyr::filter(!is.na(Gene1), !is.na(Gene2), Gene1 != "", Gene2 != "")
  return(pri)
}

prior <- read_prior_pairs(prior_file)

# ---------------------------
# Step 5) 用 Top10000 基因集合筛选先验调控对
# ---------------------------
prior_top10k <- prior %>%
  dplyr::filter(Gene1 %in% gene_set_10k, Gene2 %in% gene_set_10k) %>%
  dplyr::distinct()

write.csv(prior_top10k, out_prior_top10k_csv, row.names = FALSE, quote = TRUE)

# ---------------------------
# Step 6) 转置 Top10000 表达矩阵（samples × genes）
# ---------------------------
mat_top10k <- as.matrix(expr_top10k[, baseline_cols_raw])
rownames(mat_top10k) <- expr_top10k$ID
mat_top10k[is.na(mat_top10k)] <- 0

mat_top10k_t <- t(mat_top10k)

out_t <- data.frame(Sample = rownames(mat_top10k_t), mat_top10k_t, check.names = FALSE)
write.csv(out_t, out_expr_top10k_T_csv, row.names = FALSE, quote = TRUE)

# ---------------------------
# Step 7) 生成“数据专用 TF 列表”（保存为 TXT）
# ---------------------------
tf_df <- tryCatch(
  read.table(tf_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t"),
  error = function(e) read.table(tf_file, header = TRUE, stringsAsFactors = FALSE)
)

if (!("TF" %in% colnames(tf_df))) {
  colnames(tf_df)[1] <- "TF"
}

tf_list <- unique(tf_df$TF)

tf_in_data <- sort(intersect(tf_list, gene_set_10k))

# 保存为 TXT（单列，无列名）
write.table(
  tf_in_data,
  file = out_tf_top10k_txt,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# ---------------------------
# Step 8) 汇总信息
# ---------------------------
cat("Baseline样本数:", length(baseline_cols_raw), "\n")
cat("去掉全0基因后基因数:", nrow(expr_nonzero), "\n")
cat("Top10000 基因数:", nrow(expr_top10k), "\n")
cat("原始先验边数:", nrow(prior), "\n")
cat("Top10000筛选后先验边数:", nrow(prior_top10k), "\n")
cat("Top10000中TF数量:", length(tf_in_data), "\n\n")

cat("输出表达矩阵:", out_expr_top10k_csv, "\n")
cat("输出转置矩阵:", out_expr_top10k_T_csv, "\n")
cat("输出先验边:", out_prior_top10k_csv, "\n")
cat("输出数据专用TF列表:", out_tf_top10k_txt, "\n")





























library(data.table)
library(dplyr)
library(stringr)
library(biomaRt)

# ---------------------------
# 路径设置
# ---------------------------
expr_file  <- "D:/放假/SCD/其他数据/GSE254951/GSE254951_geo_rawcounts.txt"
prior_file <- "D:/放假/SCD/数据/构建网络的数据/label_c1.csv"
out_dir    <- "D:/放假/SCD/其他数据/GSE254951"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_expr_top10k_csv   <- file.path(out_dir, "GSE254951_preHU_Top10000_counts.csv")
out_expr_top10k_T_csv <- file.path(out_dir, "GSE254951_preHU_Top10000_counts_transposed.csv")
out_prior_top10k_csv  <- file.path(out_dir, "label_c1_prior_filtered_by_GSE254951_preHU_Top10000.csv")
out_idmap_csv         <- file.path(out_dir, "GSE254951_ensembl_to_symbol_mapping.csv")

# ---------------------------
# Step 1) 读入 rawcounts，并筛选 pre-HU 样本
# ---------------------------
df <- fread(expr_file, data.table = FALSE, check.names = FALSE)

colnames(df)[1] <- "ensembl_id"

hu_cols <- setdiff(colnames(df), "ensembl_id")
hu_num  <- as.integer(str_remove(hu_cols, "^HU"))

# pre-HU = 奇数列，排除 washout HU21, HU35
prehu_cols <- hu_cols[(hu_num %% 2 == 1) & !(hu_num %in% c(21, 35))]

expr_prehu <- df[, c("ensembl_id", prehu_cols)]

# ---------------------------
# Step 2) ID 转换：Ensembl -> Gene Symbol
# ---------------------------
expr_prehu$ensembl_id_nover <- str_remove(expr_prehu$ensembl_id, "\\.\\d+$")

mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

map <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = unique(expr_prehu$ensembl_id_nover),
  mart       = mart
)

write.csv(map, out_idmap_csv, row.names = FALSE, quote = TRUE)

# 合并映射（全部用 dplyr:: 防冲突）
expr2 <- expr_prehu %>%
  dplyr::select(-ensembl_id) %>%
  dplyr::left_join(map, by = c("ensembl_id_nover" = "ensembl_gene_id")) %>%
  dplyr::rename(GeneSymbol = hgnc_symbol)

# 去掉没有 symbol 的行
expr2 <- expr2 %>% dplyr::filter(!is.na(GeneSymbol), GeneSymbol != "")

# 如果一个 symbol 对应多个 ensembl，counts 求和
expr2_num <- expr2 %>%
  dplyr::select(GeneSymbol, all_of(prehu_cols))

expr2_num[prehu_cols] <- lapply(expr2_num[prehu_cols], as.numeric)

expr_symbol <- expr2_num %>%
  dplyr::group_by(GeneSymbol) %>%
  dplyr::summarise(across(all_of(prehu_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# ---------------------------
# Step 3) 选 Top10000 基因并转置
# ---------------------------
mat <- as.matrix(expr_symbol[, prehu_cols])
mat[is.na(mat)] <- 0

expr_symbol$TotalCount <- rowSums(mat)

expr_top10k <- expr_symbol %>%
  dplyr::arrange(desc(TotalCount)) %>%
  dplyr::slice_head(n = 10000) %>%
  dplyr::select(-TotalCount)

write.csv(expr_top10k, out_expr_top10k_csv, row.names = FALSE, quote = TRUE)

# 转置：samples x genes
mat10k <- as.matrix(expr_top10k[, prehu_cols])
rownames(mat10k) <- expr_top10k$GeneSymbol
mat10k[is.na(mat10k)] <- 0

mat10k_t <- t(mat10k)
out_t <- data.frame(Sample = rownames(mat10k_t), mat10k_t, check.names = FALSE)
write.csv(out_t, out_expr_top10k_T_csv, row.names = FALSE, quote = TRUE)

gene_set_10k <- expr_top10k$GeneSymbol

# ---------------------------
# Step 4) 读先验调控对并筛选
# ---------------------------
read_prior_pairs <- function(path) {
  pri <- tryCatch(
    read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(pri) || ncol(pri) < 2) {
    pri <- read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }
  pri <- pri[, 1:2]
  colnames(pri) <- c("Gene1", "Gene2")
  pri <- pri %>%
    dplyr::filter(!is.na(Gene1), !is.na(Gene2), Gene1 != "", Gene2 != "")
  pri
}

prior <- read_prior_pairs(prior_file)

prior_filt <- prior %>%
  dplyr::filter(Gene1 %in% gene_set_10k, Gene2 %in% gene_set_10k) %>%
  dplyr::distinct()

write.csv(prior_filt, out_prior_top10k_csv, row.names = FALSE, quote = TRUE)

# ---------------------------
# 汇总
# ---------------------------
cat("=== GSE254951 pre-HU 处理完成 ===\n")
cat("pre-HU 样本数:", length(prehu_cols), "\n")
cat("pre-HU 样本列:\n", paste(prehu_cols, collapse = ", "), "\n\n")
cat("Top10000 基因数:", nrow(expr_top10k), "\n")
cat("筛选后先验边数:", nrow(prior_filt), "\n\n")
cat("输出表达矩阵:", out_expr_top10k_csv, "\n")
cat("输出转置矩阵:", out_expr_top10k_T_csv, "\n")
cat("输出先验边:", out_prior_top10k_csv, "\n")
