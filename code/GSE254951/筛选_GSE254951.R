library(data.table)

# ============================================================
# GSE139912: 基因编号 + 先验对筛选 + 先验对编号 + 保存TXT/CSV/TSV
# 输出目录: D:/放假/SCD/其他数据/GSE139912
# ============================================================

# ---------- 0) 路径 ----------
expr_path  <- "D:/放假/SCD/其他数据/GSE139912/GSE139912_SCD_Baseline_Raw_Top10000.csv"
prior_path <- "D:/Users/29321/Desktop/new_GRN_Lung_GEN_counts_genename.csv"
out_dir    <- "D:/放假/SCD/其他数据/GSE139912"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 输出文件
gene_map_txt     <- file.path(out_dir, "gene_id_map_GSE139912.txt")                 # TXT
prior_filtered_csv <- file.path(out_dir, "new_GRN_Lung_GSE139912_filtered.csv")     # CSV
prior_ids_tsv    <- file.path(out_dir, "new_GRN_Lung_GSE139912_filtered_ids.tsv")   # TSV (no header)

# ---------- 1) 读取表达矩阵并提取基因名 ----------
# 你的表达文件第一列为 ID（基因名），分隔符可能是 tab 或逗号，这里自动识别
expr_dt <- fread(expr_path, check.names = FALSE)

# 第一列当作基因名列
gene_col <- expr_dt[[1]]
genes <- unique(gene_col[!is.na(gene_col) & gene_col != ""])

# 为保证可复现：排序后编号
genes_ord <- sort(genes)

# ---------- 2) 建立基因编号并保存 TXT ----------
gene_map <- data.table(
  name = genes_ord,
  ids  = seq.int(0, length(genes_ord) - 1)
)

# name<tab>ids
fwrite(gene_map, gene_map_txt, sep = "\t", quote = FALSE, col.names = TRUE)

# 构建 name -> id 的查找向量
id_lookup <- gene_map$ids
names(id_lookup) <- gene_map$name

# ---------- 3) 读取先验调控对（兼容有/无表头，tab/逗号） ----------
prior_try <- fread(prior_path, header = TRUE, check.names = FALSE)

if (!all(c("Gene1", "Gene2") %in% colnames(prior_try))) {
  # 无表头
  prior_dt <- fread(prior_path, header = FALSE, check.names = FALSE)
  stopifnot(ncol(prior_dt) >= 2)
  setnames(prior_dt, 1:3, c("Gene1", "Gene2", "Type"), skip_absent = TRUE)
  prior_dt <- prior_dt[, .(Gene1, Gene2, Type)]
} else {
  # 有表头
  prior_dt <- as.data.table(prior_try)
  if (!("Type" %in% colnames(prior_dt))) {
    if (ncol(prior_dt) >= 3) {
      setnames(prior_dt, 1:3, c("Gene1", "Gene2", "Type"))
    } else {
      prior_dt[, Type := NA_character_]
    }
  }
  prior_dt <- prior_dt[, .(Gene1, Gene2, Type)]
}

# 清理空值 + 防止表头被当数据（极端情况）
prior_dt <- prior_dt[!is.na(Gene1) & Gene1 != "" & !is.na(Gene2) & Gene2 != ""]
prior_dt <- prior_dt[!(Gene1 == "Gene1" & Gene2 == "Gene2")]

# ---------- 4) 按表达矩阵基因集合筛选先验对，并保存 CSV ----------
prior_filtered <- prior_dt[Gene1 %in% genes_ord & Gene2 %in% genes_ord]
fwrite(prior_filtered, prior_filtered_csv, sep = ",", quote = FALSE)

# ---------- 5) 先验对按基因编号转为 (id1, id2)，保存 TSV（无表头） ----------
prior_ids <- data.table(
  id1 = unname(id_lookup[prior_filtered$Gene1]),
  id2 = unname(id_lookup[prior_filtered$Gene2])
)

# 去掉无法映射的（理论上筛过后不该有，但保险）
prior_ids <- prior_ids[!is.na(id1) & !is.na(id2)]

# 输出 TSV：两列数字，无表头
fwrite(prior_ids, prior_ids_tsv, sep = "\t", col.names = FALSE, quote = FALSE)

# ---------- 6) 汇总 ----------
cat("✅ 完成！输出目录：", out_dir, "\n\n")
cat("表达矩阵基因数：", length(genes_ord), "\n")
cat("先验对（筛选后）条数：", nrow(prior_filtered), "\n")
cat("先验对（编号后）条数：", nrow(prior_ids), "\n\n")
cat("输出文件：\n")
cat(" - 基因ID编号(TXT): ", gene_map_txt, "\n")
cat(" - 先验对筛选(CSV): ", prior_filtered_csv, "\n")
cat(" - 先验对编号(TSV): ", prior_ids_tsv, "\n")




















library(data.table)

# ============================================================
# GSE254951: 基因编号 + 先验对筛选 + 先验对编号 + 保存TXT/CSV/TSV
# 输出目录: D:/放假/SCD/其他数据/GSE254951
# ============================================================

# ---------- 0) 路径 ----------
expr_path  <- "D:/放假/SCD/其他数据/GSE254951/GSE254951_preHU_Top10000_counts.csv"
prior_path <- "D:/Users/29321/Desktop/new_GRN_Lung_GEN_counts_genename.csv"
out_dir    <- "D:/放假/SCD/其他数据/GSE254951"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 输出文件
gene_map_txt       <- file.path(out_dir, "gene_id_map_GSE254951.txt")                # TXT
prior_filtered_csv <- file.path(out_dir, "new_GRN_Lung_GSE254951_filtered.csv")      # CSV
prior_ids_tsv      <- file.path(out_dir, "new_GRN_Lung_GSE254951_filtered_ids.tsv")  # TSV (no header)

# ---------- 1) 读取表达矩阵并提取基因名 ----------
# 文件可能是 tab 或逗号分隔，fread 通常可自动识别
expr_dt <- fread(expr_path, check.names = FALSE)

# 第一列为 GeneSymbol（基因名）
genes <- unique(expr_dt[[1]])
genes <- genes[!is.na(genes) & genes != ""]

# 为保证可复现：排序后编号
genes_ord <- sort(genes)

# ---------- 2) 建立基因编号并保存 TXT ----------
gene_map <- data.table(
  name = genes_ord,
  ids  = seq.int(0, length(genes_ord) - 1)
)

# 保存：name<TAB>ids（带表头）
fwrite(gene_map, gene_map_txt, sep = "\t", quote = FALSE, col.names = TRUE)

# 构建 name -> id 的查找向量
id_lookup <- gene_map$ids
names(id_lookup) <- gene_map$name

# ---------- 3) 读取先验调控对（兼容有/无表头，tab/逗号） ----------
prior_try <- fread(prior_path, header = TRUE, check.names = FALSE)

if (!all(c("Gene1", "Gene2") %in% colnames(prior_try))) {
  # 无表头
  prior_dt <- fread(prior_path, header = FALSE, check.names = FALSE)
  stopifnot(ncol(prior_dt) >= 2)
  setnames(prior_dt, 1:3, c("Gene1", "Gene2", "Type"), skip_absent = TRUE)
  prior_dt <- prior_dt[, .(Gene1, Gene2, Type)]
} else {
  # 有表头
  prior_dt <- as.data.table(prior_try)
  if (!("Type" %in% colnames(prior_dt))) {
    if (ncol(prior_dt) >= 3) {
      setnames(prior_dt, 1:3, c("Gene1", "Gene2", "Type"))
    } else {
      prior_dt[, Type := NA_character_]
    }
  }
  prior_dt <- prior_dt[, .(Gene1, Gene2, Type)]
}

# 清理空值 + 极端情况：删除表头被当数据的一行
prior_dt <- prior_dt[!is.na(Gene1) & Gene1 != "" & !is.na(Gene2) & Gene2 != ""]
prior_dt <- prior_dt[!(Gene1 == "Gene1" & Gene2 == "Gene2")]

# ---------- 4) 按表达矩阵基因集合筛选先验对，并保存 CSV ----------
prior_filtered <- prior_dt[Gene1 %in% genes_ord & Gene2 %in% genes_ord]
fwrite(prior_filtered, prior_filtered_csv, sep = ",", quote = FALSE)

# ---------- 5) 先验对按基因编号转为 (id1, id2)，保存 TSV（无表头） ----------
prior_ids <- data.table(
  id1 = unname(id_lookup[prior_filtered$Gene1]),
  id2 = unname(id_lookup[prior_filtered$Gene2])
)

# 保险：去掉 NA
prior_ids <- prior_ids[!is.na(id1) & !is.na(id2)]

# 输出 TSV：两列数字，无表头
fwrite(prior_ids, prior_ids_tsv, sep = "\t", col.names = FALSE, quote = FALSE)

# ---------- 6) 汇总 ----------
cat("✅ 完成！输出目录：", out_dir, "\n\n")
cat("表达矩阵基因数：", length(genes_ord), "\n")
cat("先验对（筛选后）条数：", nrow(prior_filtered), "\n")
cat("先验对（编号后）条数：", nrow(prior_ids), "\n\n")
cat("输出文件：\n")
cat(" - 基因ID编号(TXT): ", gene_map_txt, "\n")
cat(" - 先验对筛选(CSV): ", prior_filtered_csv, "\n")
cat(" - 先验对编号(TSV): ", prior_ids_tsv, "\n")

library(data.table)

# =========================
# 0) 路径（按需改）
# =========================
grn_path <- "D:/放假/SCD/其他数据/GSE254951/GSE254951(deepfgrn五折).tsv"
map_path <- "D:/放假/SCD/其他数据/GSE254951/gene_id_map_GSE254951.txt"

out_all_path   <- "D:/放假/SCD/其他数据/GSE254951/GSE254951_deepfgrn_gene_EdgeWeight_all.tsv"
out_top10_path <- "D:/放假/SCD/其他数据/GSE254951/GSE254951_deepfgrn_gene_EdgeWeight_top10pct.tsv"

# =========================
# 1) 读入 GRN 结果 + 基因编号表
# =========================
grn <- fread(grn_path, sep = "\t", header = TRUE)
gm  <- fread(map_path, sep = "\t", header = TRUE)

stopifnot(all(c("TF_id","Target_id","prob_no","prob_act","prob_rep","pred_label") %in% names(grn)))
stopifnot(all(c("name","ids") %in% names(gm)))

# 构建 id -> gene 的查找向量
setnames(gm, c("name","ids"), c("gene","id"))
id2gene <- gm$gene
names(id2gene) <- as.character(gm$id)

# =========================
# 2) 将数字ID映射回基因名
# =========================
grn[, TF_id := as.integer(TF_id)]
grn[, Target_id := as.integer(Target_id)]

grn[, TF     := id2gene[as.character(TF_id)]]
grn[, Target := id2gene[as.character(Target_id)]]

if (anyNA(grn$TF) || anyNA(grn$Target)) {
  warning("存在未映射到基因名的ID：TF NA=", sum(is.na(grn$TF)),
          "；Target NA=", sum(is.na(grn$Target)))
}

# =========================
# 3) 三列概率合并为一列 EdgeWeight（方案B）
#    act: +prob_act
#    rep: -prob_rep
#    no :  (1 - prob_no)
# =========================
grn[, pred_label := as.integer(pred_label)]

grn[, EdgeWeight := fifelse(pred_label == 1L,  prob_act,
                            fifelse(pred_label == 2L, -prob_rep,
                                    (1 - prob_no)))]

# =========================
# 4) 生成三列并改列名：TF  Target  EdgeWeight
# =========================
out_all <- grn[, .(TF, Target, EdgeWeight)]

# 保存全部边
fwrite(out_all, out_all_path, sep = "\t", quote = FALSE)

# =========================
# 5) 筛选 |EdgeWeight| 前 10% 的调控边
# =========================
thr <- quantile(abs(out_all$EdgeWeight), probs = 0.9, na.rm = TRUE, type = 7)
out_top10 <- out_all[abs(EdgeWeight) >= thr]

fwrite(out_top10, out_top10_path, sep = "\t", quote = FALSE)

cat("✅ 完成！\n")
cat("全部边数：", nrow(out_all), "\n")
cat("Top10% 阈值 |EdgeWeight| >= ", thr, "\n", sep = "")
cat("Top10% 边数：", nrow(out_top10), "\n")
cat("输出文件：\n - ", out_all_path, "\n - ", out_top10_path, "\n", sep = "")