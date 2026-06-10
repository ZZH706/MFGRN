# =========================================================
#  Poisson 筛选关键 TF —— 边列表批处理（带中文标签＋TF 白名单）
#  读取四张 “TF Target” 边表，结合 Homo_TF_clean.txt 过滤，
#  输出 *_poisson_all.txt / *_poisson_TF.txt 到指定目录，
#  并在控制台打印中文概要信息。
# =========================================================
library(data.table)
library(igraph)

## ---------- 1. 输入文件列表（命名便于映射） ----------
edge_files <- c(
  c1 = "E:/1数据/≥2/evaluation_results-c1/two_or_more_overlap_c1_TF_TARGET.txt",
  c2 = "E:/1数据/≥2/evaluation_results-c2/two_or_more_overlap_c2_TF_TARGET.txt",
  s1 = "E:/1数据/≥2/evaluation_results-s1/two_or_more_overlap_s1_TF_TARGET.txt",
  s2 = "E:/1数据/≥2/evaluation_results-s2/two_or_more_overlap_s2_TF_TARGET.txt"
)

## ---------- 2. 中文标签 ----------
label_map <- c(
  c1 = "对照组运动前",
  c2 = "对照组运动后",
  s1 = "患病组运动前",
  s2 = "患病组运动后"
)

## ---------- 3. 输出目录 ----------
out_dir <- "E:/1数据/泊松分布筛选关键转录因子"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- 4. 参数 ----------
alpha     <- 0.05      # FDR 阈值
tail_prob <- 0.999     # 度阈值分位

## ---------- 5. 读取转录因子白名单 ----------
tf_list <- fread("E:/1数据/Homo_TF_clean.txt",
                 header = TRUE, sep = "\t", colClasses = "character")[, TF]

## ---------- 6. 主循环 ----------
for (tag in names(edge_files)) {
  
  ef <- edge_files[[tag]]
  
  # ---- 6.1 读边列表 ----
  edges <- fread(ef, header = TRUE, sep = "\t", colClasses = "character")
  setnames(edges, 1:2, c("TF", "Target"))   # 统一列名
  
  # ---- 6.2 建图 ----
  g <- graph_from_data_frame(edges, directed = TRUE)
  
  # ---- 6.3 基本统计量 ----
  n <- vcount(g)               # 节点数
  m <- ecount(g)               # 边数
  lambda <- m / n              # 期望出度 (≈ (n−1)*p)
  
  # ---- 6.4 只对 TF 白名单节点计算出度 ----
  deg_out <- degree(g, mode = "out")
  deg_tf  <- deg_out[names(deg_out) %in% tf_list]
  
  # ---- 6.5 Poisson 右尾 p & FDR ----
  p_raw <- 1 - ppois(deg_tf - 1, lambda = lambda)
  q_val <- p.adjust(p_raw, method = "BH")
  
  # ---- 6.6 判定显著 ----
  deg_cut <- qpois(tail_prob, lambda = lambda)
  sig_idx <- which(q_val < alpha & deg_tf >= deg_cut)
  
  # ---- 6.7 保存结果 ----
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
  
  # ---- 6.8 控制台提示 ----
  cat(sprintf(
    "\n[%s] λ = %.2f, 度阈 = %d, 关键 TF = %d\n",
    label_map[tag], lambda, deg_cut, nrow(res_sig)
  ))
}

cat("\n全部网络已完成，结果保存在：", out_dir, "\n")





# ========== 7. GO 富集分析 ==========
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(data.table)

# 创建 GO 富集结果输出目录
go_dir <- file.path(out_dir, "GO_富集分析")
dir.create(go_dir, showWarnings = FALSE)

# 中文标签
label_map <- c(
  c1 = "对照组运动前",
  c2 = "对照组运动后",
  s1 = "患病组运动前",
  s2 = "患病组运动后"
)

# 英文标签（用于图标题）
label_map_en <- c(
  c1 = "Control_Pre",
  c2 = "Control_Post",
  s1 = "SCA_Pre",
  s2 = "SCA_Post"
)

# 批量对每个网络的关键 TF 做 GO 分析
for (tag in names(edge_files)) {
  
  net_tag <- sub("\\.txt$", "", basename(edge_files[[tag]]))
  tf_file <- file.path(out_dir, paste0(net_tag, "_poisson_TF.txt"))
  
  if (!file.exists(tf_file)) next
  
  sig_tfs <- fread(tf_file)
  tf_symbols <- sig_tfs$Gene
  
  # 转换为 ENTREZ ID（GO 分析需要）
  gene_df <- bitr(tf_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- unique(gene_df$ENTREZID)
  
  # 若数量太少，跳过
  if (length(entrez_ids) < 3) {
    cat(sprintf("[%s] 关键 TF 数太少，跳过 GO 分析\n", label_map[tag]))
    next
  }
  
  # GO Biological Process 富集分析
  ego <- enrichGO(gene         = entrez_ids,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "ENTREZID",
                  ont          = "BP",
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable     = TRUE)
  
  # 保存表格
  out_txt <- file.path(go_dir, paste0(net_tag, "_GO.txt"))
  fwrite(as.data.table(ego), out_txt, sep = "\t", quote = FALSE)
  
  # 保存条形图（标题用英文）
  out_pdf <- file.path(go_dir, paste0(net_tag, "_GO_barplot.pdf"))
  pdf(out_pdf, width = 8, height = 6)
  print(barplot(ego, showCategory = 20, title = paste0("GO Enrichment - ", label_map_en[tag])))
  dev.off()
  
  cat(sprintf("[%s] GO 分析完成，结果保存于 %s\n", label_map[tag], go_dir))
}

cat("\n全部 GO 富集分析已完成，结果保存在：", go_dir, "\n")


library(gridExtra)

# 存储所有 barplot 图对象
go_plots <- list()

# 重新走一遍（或在原分析循环中保存图对象）
for (tag in names(edge_files)) {
  
  net_tag <- sub("\\.txt$", "", basename(edge_files[[tag]]))
  tf_file <- file.path(out_dir, paste0(net_tag, "_poisson_TF.txt"))
  if (!file.exists(tf_file)) next
  
  sig_tfs <- fread(tf_file)
  tf_symbols <- sig_tfs$Gene
  gene_df <- bitr(tf_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- unique(gene_df$ENTREZID)
  
  if (length(entrez_ids) < 3) next
  
  ego <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  
  # 存储 barplot 图对象（这里不输出单独 PDF）
  go_plots[[tag]] <- barplot(ego, showCategory = 20, title = paste0("GO Enrichment - ", label_map_en[tag]))
}

# 合并图像（2 行 2 列）
pdf(file.path(go_dir, "Combined_GO_barplot.pdf"), width = 14, height = 10)
grid.arrange(grobs = go_plots, ncol = 2)
dev.off()




library(gridExtra)
library(ggplot2)

# 存储所有 barplot 图对象
go_plots <- list()

for (tag in names(edge_files)) {
  
  net_tag <- sub("\\.txt$", "", basename(edge_files[[tag]]))
  tf_file <- file.path(out_dir, paste0(net_tag, "_poisson_TF.txt"))
  if (!file.exists(tf_file)) next
  
  sig_tfs <- fread(tf_file)
  tf_symbols <- sig_tfs$Gene
  gene_df <- bitr(tf_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- unique(gene_df$ENTREZID)
  
  if (length(entrez_ids) < 3) next
  
  ego <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  
  go_plots[[tag]] <- barplot(ego, showCategory = 20, title = paste0("GO Enrichment - ", label_map_en[tag])) +
    theme(axis.text.y = element_text(size = 8))  # 调小纵轴标签字体
}

pdf(file.path(go_dir, "Combined_GO_barplot.pdf"), width = 14, height = 10)
grid.arrange(grobs = go_plots, ncol = 2)
dev.off()


library(gridExtra)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)

# 存储所有 barplot 图对象
kegg_plots <- list()

for (tag in names(edge_files)) {
  
  net_tag <- sub("\\.txt$", "", basename(edge_files[[tag]]))
  tf_file <- file.path(out_dir, paste0(net_tag, "_poisson_TF.txt"))
  if (!file.exists(tf_file)) next
  
  sig_tfs <- fread(tf_file)
  tf_symbols <- sig_tfs$Gene
  # 转换基因ID为ENTREZID，KEGG分析需要ENTREZID
  gene_df <- bitr(tf_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- unique(gene_df$ENTREZID)
  
  if (length(entrez_ids) < 3) next
  
  # KEGG富集分析
  ekegg <- enrichKEGG(gene = entrez_ids,
                      organism = "hsa",           # "hsa"表示人类
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
  
  # 如果想转换结果中的ID为基因名（可选）
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  # 保存KEGG条形图，调小字体防止标签重叠
  kegg_plots[[tag]] <- barplot(ekegg, showCategory = 20, title = paste0("KEGG Enrichment - ", label_map_en[tag])) +
    theme(axis.text.y = element_text(size = 6))
}

pdf(file.path(go_dir, "Combined_KEGG_barplot.pdf"), width = 14, height = 10)
grid.arrange(grobs = kegg_plots, ncol = 2)
dev.off()



library(gridExtra)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)

# 存储所有 barplot 图对象
kegg_plots <- list()

for (tag in names(edge_files)) {
  
  net_tag <- sub("\\.txt$", "", basename(edge_files[[tag]]))
  tf_file <- file.path(out_dir, paste0(net_tag, "_poisson_TF.txt"))
  if (!file.exists(tf_file)) next
  
  sig_tfs <- fread(tf_file)
  tf_symbols <- sig_tfs$Gene
  gene_df <- bitr(tf_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- unique(gene_df$ENTREZID)
  
  if (length(entrez_ids) < 3) next
  
  ekegg <- enrichKEGG(gene = entrez_ids,
                      organism = "hsa",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
  
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  kegg_plots[[tag]] <- barplot(ekegg,
                               showCategory = 20,
                               font.size = 11,
                               label_format = 50,  # ⬅️ 自动换行：每行最多40字符
                               title = paste0("KEGG Enrichment - ", label_map_en[tag])) +
    scale_x_continuous(breaks = seq(0, 60, 20)) +
    theme(axis.text.y = element_text(size = 8),
          plot.margin = unit(c(1, 1, 1, 3), "cm"))
  
}

# 保持整张图PDF大小不变
pdf(file.path(go_dir, "Combined_KEGG_barplot.pdf"), width = 14, height = 10)
grid.arrange(grobs = kegg_plots, ncol = 2)
dev.off()










# =========================================================
#  Poisson 筛选关键 TF + GO 三类合并条形图富集分析
# =========================================================

library(data.table)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

## ---------- 1. 输入文件列表（命名便于映射） ----------
edge_files <- c(
  c1 = "E:/1数据/≥2/evaluation_results-c1/two_or_more_overlap_c1_TF_TARGET.txt",
  c2 = "E:/1数据/≥2/evaluation_results-c2/two_or_more_overlap_c2_TF_TARGET.txt",
  s1 = "E:/1数据/≥2/evaluation_results-s1/two_or_more_overlap_s1_TF_TARGET.txt",
  s2 = "E:/1数据/≥2/evaluation_results-s2/two_or_more_overlap_s2_TF_TARGET.txt"
)

## ---------- 2. 中文标签 ----------
label_map <- c(
  c1 = "对照组运动前",
  c2 = "对照组运动后",
  s1 = "患病组运动前",
  s2 = "患病组运动后"
)

## ---------- 3. 输出目录 ----------
out_dir <- "E:/1数据/泊松分布筛选关键转录因子"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- 4. 参数 ----------
alpha     <- 0.05      # FDR 阈值
tail_prob <- 0.999     # 度阈值分位

## ---------- 5. 读取转录因子白名单 ----------
tf_list <- fread("E:/1数据/Homo_TF_clean.txt",
                 header = TRUE, sep = "\t", colClasses = "character")[, TF]

## ---------- 6. Poisson 筛选关键 TF 主循环 ----------
for (tag in names(edge_files)) {
  
  ef <- edge_files[[tag]]
  
  # 读边列表
  edges <- fread(ef, header = TRUE, sep = "\t", colClasses = "character")
  setnames(edges, 1:2, c("TF", "Target"))   # 统一列名
  
  # 建图
  g <- graph_from_data_frame(edges, directed = TRUE)
  
  # 基本统计量
  n <- vcount(g)
  m <- ecount(g)
  lambda <- m / n
  
  # 只对 TF 白名单节点计算出度
  deg_out <- degree(g, mode = "out")
  deg_tf  <- deg_out[names(deg_out) %in% tf_list]
  
  # Poisson 右尾 p & FDR
  p_raw <- 1 - ppois(deg_tf - 1, lambda = lambda)
  q_val <- p.adjust(p_raw, method = "BH")
  
  # 判定显著阈值
  deg_cut <- qpois(tail_prob, lambda = lambda)
  sig_idx <- which(q_val < alpha & deg_tf >= deg_cut)
  
  # 保存结果
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
  
  cat(sprintf("\n[%s] λ = %.2f, 度阈 = %d, 关键 TF = %d\n",
              label_map[tag], lambda, deg_cut, nrow(res_sig)))
}

cat("\n全部网络Poisson筛选已完成，结果保存在：", out_dir, "\n")



## ---------- 7. GO 富集三类合并绘图 ----------

go_dir <- file.path(out_dir, "GO_富集分析")
dir.create(go_dir, showWarnings = FALSE)

label_map_en <- c(
  c1 = "Control_Pre",
  c2 = "Control_Post",
  s1 = "SCA_Pre",
  s2 = "SCA_Post"
)

all_go_results <- list()

for (tag in names(edge_files)) {
  
  net_tag <- sub("\\.txt$", "", basename(edge_files[[tag]]))
  tf_file <- file.path(out_dir, paste0(net_tag, "_poisson_TF.txt"))
  if (!file.exists(tf_file)) next
  
  sig_tfs <- fread(tf_file)
  tf_symbols <- sig_tfs$Gene
  
  gene_df <- bitr(tf_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- unique(gene_df$ENTREZID)
  
  if (length(entrez_ids) < 3) {
    cat(sprintf("[%s] 关键 TF 数太少，跳过 GO 分析\n", label_map[tag]))
    next
  }
  
  ontologies <- c("BP", "CC", "MF")
  go_res_list <- list()
  
  for (ont in ontologies) {
    ego <- enrichGO(gene         = entrez_ids,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    ont          = ont,
                    pAdjustMethod= "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable     = TRUE)
    
    if (is.null(ego) || nrow(ego) == 0) next
    
    df_ego <- as.data.frame(ego)
    df_ego$Ontology <- ont
    df_ego$Group <- label_map_en[tag]
    
    go_res_list[[ont]] <- df_ego
  }
  
  if (length(go_res_list) == 0) next
  
  go_all <- dplyr::bind_rows(go_res_list)
  
  fwrite(as.data.table(go_all),
         file.path(go_dir, paste0(net_tag, "_GO_all_BP_CC_MF.txt")),
         sep = "\t", quote = FALSE)
  
  all_go_results[[tag]] <- go_all
  
  cat(sprintf("[%s] GO分析完成，数据已保存\n", label_map[tag]))
}

go_combined_df <- dplyr::bind_rows(all_go_results)

if (nrow(go_combined_df) == 0) stop("没有有效的GO富集结果，无法绘图！")

go_top <- go_combined_df %>%
  group_by(Group, Ontology) %>%
  slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup()

go_top <- go_top %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

p <- ggplot(go_top, aes(x = -log10(p.adjust), y = Description, fill = Ontology)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ Group, scales = "free_y") +
  scale_fill_manual(values = c(BP = "#1b9e77", CC = "#d95f02", MF = "#7570b3")) +
  labs(x = expression(-log[10](adjusted~p-value)),
       y = NULL,
       fill = "GO Category",
       title = "GO Enrichment Across Groups") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 10),
        legend.position = "bottom")

pdf(file.path(go_dir, "Combined_GO_Barplot_BP_CC_MF.pdf"), width = 14, height = 10)
print(p)
dev.off()

cat("\nGO三类合并条形图已保存至：", file.path(go_dir, "Combined_GO_Barplot_BP_CC_MF.pdf"), "\n")
