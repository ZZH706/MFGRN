## =========================
## 0) 加载依赖（edgeR/limma 用 Bioconductor 安装）
## =========================
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR","limma"), ask=FALSE, update=FALSE)
install.packages(c("data.table","ggplot2","ggrepel"), dependencies=TRUE)

library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)

## =========================
## 1) 输入文件路径
## =========================
f_con_t1 <- "E:/SCD/数据/原始数据/expr_CON_T1.csv"
f_con_t2 <- "E:/SCD/数据/原始数据/expr_CON_T2.csv"
f_sca_t1 <- "E:/SCD/数据/原始数据/expr_SCA_T1.csv"
f_sca_t2 <- "E:/SCD/数据/原始数据/expr_SCA_T2.csv"

tf_file  <- "E:/SCD/数据/构建网络的数据/human_tf_list.txt"

## =========================
## 2) 读取表达矩阵（第一列为基因名，后面为样本列）
## =========================
read_expr <- function(path){
  x <- fread(path)
  gene <- x[[1]]
  mat <- as.matrix(x[, -1, with=FALSE])
  rownames(mat) <- make.unique(as.character(gene))
  colnames(mat) <- make.unique(colnames(mat))
  mat
}

con_t1 <- read_expr(f_con_t1)
con_t2 <- read_expr(f_con_t2)
sca_t1 <- read_expr(f_sca_t1)
sca_t2 <- read_expr(f_sca_t2)

genes <- Reduce(intersect, list(rownames(con_t1), rownames(con_t2), rownames(sca_t1), rownames(sca_t2)))
con_t1 <- con_t1[genes, , drop=FALSE]
con_t2 <- con_t2[genes, , drop=FALSE]
sca_t1 <- sca_t1[genes, , drop=FALSE]
sca_t2 <- sca_t2[genes, , drop=FALSE]

expr <- cbind(con_t1, con_t2, sca_t1, sca_t2)

## =========================
## 2.5) 读取 TF 列表，并对齐命名
## =========================
tf_list <- fread(tf_file, header=FALSE)$V1
tf_list <- unique(toupper(tf_list))

rownames(expr) <- toupper(rownames(expr))
tf_in_expr <- intersect(tf_list, rownames(expr))

message("TFs in list: ", length(tf_list))
message("TFs found in expr: ", length(tf_in_expr))
if(length(tf_in_expr) < 20){
  message("WARNING: 匹配到的 TF 很少；若表达矩阵是ENSG ID，需要先映射到symbol。")
}

## =========================
## 3) limma/voom 拟合（全基因拟合，后面只取TF结果）
## =========================
meta <- data.frame(
  sample = colnames(expr),
  Group  = c(rep("CON", ncol(con_t1)), rep("CON", ncol(con_t2)),
             rep("SCA", ncol(sca_t1)), rep("SCA", ncol(sca_t2))),
  Time   = c(rep("T1", ncol(con_t1)), rep("T2", ncol(con_t2)),
             rep("T1", ncol(sca_t1)), rep("T2", ncol(sca_t2))),
  stringsAsFactors = FALSE
)
meta$Condition <- factor(paste(meta$Group, meta$Time, sep="_"),
                         levels=c("CON_T1","CON_T2","SCA_T1","SCA_T2"))

is_counts <- all(abs(expr - round(expr)) < 1e-6, na.rm=TRUE) && max(expr, na.rm=TRUE) > 50
message("Detected counts? ", is_counts)

design <- model.matrix(~0 + meta$Condition)
colnames(design) <- levels(meta$Condition)

if(is_counts){
  dge <- DGEList(expr)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
} else {
  fit <- lmFit(log2(expr + 1), design)
}

contr <- makeContrasts(
  baseline     = SCA_T1 - CON_T1,
  exercise_CON = CON_T2 - CON_T1,
  exercise_SCA = SCA_T2 - SCA_T1,
  interaction  = (SCA_T2 - SCA_T1) - (CON_T2 - CON_T1),
  levels = design
)

fit2 <- eBayes(contrasts.fit(fit, contr))

res_baseline <- topTable(fit2, coef="baseline",     number=Inf, adjust.method="BH", sort.by="P")
res_con_ex   <- topTable(fit2, coef="exercise_CON", number=Inf, adjust.method="BH", sort.by="P")
res_sca_ex   <- topTable(fit2, coef="exercise_SCA", number=Inf, adjust.method="BH", sort.by="P")
res_inter    <- topTable(fit2, coef="interaction",  number=Inf, adjust.method="BH", sort.by="P")

## 只保留 TF
filter_to_tf <- function(tt, tf_set){
  tt$gene <- rownames(tt)
  keep <- toupper(tt$gene) %in% tf_set
  tt <- tt[keep, , drop=FALSE]
  rownames(tt) <- tt$gene
  tt$gene <- NULL
  tt
}
res_baseline_tf <- filter_to_tf(res_baseline, tf_in_expr)
res_con_ex_tf   <- filter_to_tf(res_con_ex,   tf_in_expr)
res_sca_ex_tf   <- filter_to_tf(res_sca_ex,   tf_in_expr)
res_inter_tf    <- filter_to_tf(res_inter,    tf_in_expr)

## =========================
## 4) 火山图函数（不保存，直接展示）
## =========================
plot_volcano_tf <- function(tt, title,
                            fdr_cut=0.05, lfc_cut=0.5,
                            label_top=20){
  if(nrow(tt) == 0) stop("No TFs to plot for: ", title)
  
  df <- tt
  df$gene <- rownames(df)
  df$neglog10FDR <- -log10(df$adj.P.Val)
  df$Sig <- "NS"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC >=  lfc_cut] <- "Up"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC <= -lfc_cut] <- "Down"
  
  # 标注：优先显著且|lfc|满足阈值，不足则补充最小FDR
  lab <- df[df$adj.P.Val < fdr_cut & abs(df$logFC) >= lfc_cut, , drop=FALSE]
  lab <- lab[order(lab$adj.P.Val), , drop=FALSE]
  if(nrow(lab) < label_top){
    extra <- df[order(df$adj.P.Val), , drop=FALSE]
    extra <- extra[!extra$gene %in% lab$gene, , drop=FALSE]
    lab <- rbind(lab, head(extra, label_top - nrow(lab)))
  } else {http://127.0.0.1:9347/graphics/plot_zoom_png?width=982&height=861
    lab <- head(lab, label_top)
  }
  
  ggplot(df, aes(x=logFC, y=neglog10FDR)) +
    geom_point(aes(color=Sig), alpha=0.85, size=2.2) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr_cut), linetype="dashed") +
    scale_color_manual(values=c(NS="grey70", Up="red3", Down="royalblue3")) +
    labs(title=title, x="log2 Fold Change", y="-log10(FDR)") +
    theme_bw(base_size=12) +
    theme(plot.title = element_text(hjust=0.5)) +
    ggrepel::geom_text_repel(
      data=lab,
      aes(label=gene),
      max.overlaps=Inf,
      size=3
    )
}

## =========================
## 5) 依次展示 4 张 TF 火山图
## =========================
p1 <- plot_volcano_tf(res_baseline_tf, "TF Volcano: SCA_T1 vs CON_T1")
print(p1)

p2 <- plot_volcano_tf(res_con_ex_tf, "TF Volcano: CON_T2 vs CON_T1 (Exercise response)")
print(p2)

p3 <- plot_volcano_tf(res_sca_ex_tf, "TF Volcano: SCA_T2 vs SCA_T1 (Exercise response)")
print(p3)

p4 <- plot_volcano_tf(res_inter_tf, "TF Volcano: Interaction (DiD)")
print(p4)






## =========================
## 0) 加载依赖（edgeR/limma 用 Bioconductor 安装）
## =========================
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR","limma"), ask=FALSE, update=FALSE)
install.packages(c("data.table","ggplot2","ggrepel"), dependencies=TRUE)

library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)

## =========================
## 1) 输入文件路径
## =========================
f_con_t1 <- "E:/SCD/数据/原始数据/expr_CON_T1.csv"
f_con_t2 <- "E:/SCD/数据/原始数据/expr_CON_T2.csv"
f_sca_t1 <- "E:/SCD/数据/原始数据/expr_SCA_T1.csv"
f_sca_t2 <- "E:/SCD/数据/原始数据/expr_SCA_T2.csv"

## =========================
## 2) 你指定的 TF 列表（只画这些）
## =========================
tf_focus <- c(
  "AHR","ATF4","ATF5","BHLHE40","BHLHE41","CLOCK","CREB1","CREM","DBP",
  "EGR1","EGR3","ESR1","FOXO3","HNF4A","JUND","KLF10","KLF9","NFIL3",
  "NFKB2","NFYA","NPAS2","NCOA1","NR1H3","NR5A1","NR5A2","PPARA","PPARG",
  "PROX1","RELB","RORA","RORC","SP1","SREBF1","STAT5A","STAT5B","TEF",
  "TP53","TWIST1","ZFHX3"
)
tf_focus <- unique(toupper(tf_focus))

## =========================
## 3) 读取表达矩阵（第一列为基因名，后面为样本列）
## =========================
read_expr <- function(path){
  x <- fread(path)
  gene <- x[[1]]
  mat <- as.matrix(x[, -1, with=FALSE])
  rownames(mat) <- make.unique(as.character(gene))
  colnames(mat) <- make.unique(colnames(mat))
  mat
}

con_t1 <- read_expr(f_con_t1)
con_t2 <- read_expr(f_con_t2)
sca_t1 <- read_expr(f_sca_t1)
sca_t2 <- read_expr(f_sca_t2)

genes <- Reduce(intersect, list(rownames(con_t1), rownames(con_t2), rownames(sca_t1), rownames(sca_t2)))
con_t1 <- con_t1[genes, , drop=FALSE]
con_t2 <- con_t2[genes, , drop=FALSE]
sca_t1 <- sca_t1[genes, , drop=FALSE]
sca_t2 <- sca_t2[genes, , drop=FALSE]

expr <- cbind(con_t1, con_t2, sca_t1, sca_t2)

## 统一大小写，便于与 tf_focus 匹配
rownames(expr) <- toupper(rownames(expr))

tf_in_expr <- intersect(tf_focus, rownames(expr))
message("TFs requested: ", length(tf_focus))
message("TFs found in expr: ", length(tf_in_expr))
if(length(tf_in_expr) == 0) stop("一个都没匹配到：请检查表达矩阵行名是否为基因symbol（而非ENSG）。")

## =========================
## 4) limma/voom 拟合（全基因拟合，但最后只取指定TF）
## =========================
meta <- data.frame(
  sample = colnames(expr),
  Condition = factor(c(
    rep("CON_T1", ncol(con_t1)),
    rep("CON_T2", ncol(con_t2)),
    rep("SCA_T1", ncol(sca_t1)),
    rep("SCA_T2", ncol(sca_t2))
  ), levels=c("CON_T1","CON_T2","SCA_T1","SCA_T2"))
)

is_counts <- all(abs(expr - round(expr)) < 1e-6, na.rm=TRUE) && max(expr, na.rm=TRUE) > 50
message("Detected counts? ", is_counts)

design <- model.matrix(~0 + meta$Condition)
colnames(design) <- levels(meta$Condition)

if(is_counts){
  dge <- DGEList(expr)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
} else {
  fit <- lmFit(log2(expr + 1), design)
}

contr <- makeContrasts(
  baseline = SCA_T1 - CON_T1,
  levels = design
)

fit2 <- eBayes(contrasts.fit(fit, contr))
res_baseline <- topTable(fit2, coef="baseline", number=Inf, adjust.method="BH", sort.by="P")

## 只保留你指定的 TF
res_baseline$gene <- toupper(rownames(res_baseline))
res_tf <- res_baseline[res_baseline$gene %in% tf_in_expr, , drop=FALSE]
rownames(res_tf) <- res_tf$gene

## =========================
## 5) 火山图：只画指定TF（不保存，只展示）
## =========================
plot_volcano_focus <- function(tt, title, fdr_cut=0.05, lfc_cut=0.5){
  df <- tt
  df$gene <- rownames(df)
  df$neglog10FDR <- -log10(df$adj.P.Val)
  
  df$Sig <- "NS"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC >=  lfc_cut] <- "Up"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC <= -lfc_cut] <- "Down"
  
  # 标注：对这批TF全部标出来（数量不大）
  ggplot(df, aes(x=logFC, y=neglog10FDR)) +
    geom_point(aes(color=Sig), alpha=0.9, size=3) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr_cut), linetype="dashed") +
    scale_color_manual(values=c(NS="grey70", Up="red3", Down="royalblue3")) +
    labs(title=title, x="log2 Fold Change (SCA_T1 / CON_T1)", y="-log10(FDR)") +
    theme_bw(base_size=12) +
    theme(plot.title = element_text(hjust=0.5)) +
    ggrepel::geom_text_repel(
      aes(label=gene),
      max.overlaps=Inf,
      size=3
    )
}

p <- plot_volcano_focus(
  res_tf,
  title = "TF Volcano (Focused set): SCA_T1 vs CON_T1",
  fdr_cut = 0.05,
  lfc_cut = 0.5
)
print(p)



## =========================
## 0) 加载依赖（edgeR/limma 用 Bioconductor 安装）
## =========================
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR","limma"), ask=FALSE, update=FALSE)
install.packages(c("data.table","ggplot2","ggrepel"), dependencies=TRUE)

library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)

## =========================
## 1) 输入文件路径
## =========================
f_con_t1 <- "E:/SCD/数据/原始数据/expr_CON_T1.csv"
f_con_t2 <- "E:/SCD/数据/原始数据/expr_CON_T2.csv"
f_sca_t1 <- "E:/SCD/数据/原始数据/expr_SCA_T1.csv"
f_sca_t2 <- "E:/SCD/数据/原始数据/expr_SCA_T2.csv"

## =========================
## 2) 你指定的 TF 列表（只画这些）
## =========================
tf_focus <- c(
  "CEBPA","CREB1","ELK1","EPAS1","FOSL2","GATA6","GLI1","GLI3","FOXA1",
  "HOXA5","HES1","RBPJ","MSX1","MYCN","NFIB","PROX1","RARA","RARG","SIX1",
  "SOX11","SP3","SREBF1","SRF","TBX2","TCF21","THRA","THRB","WT1","HMGA2",
  "MSC","KLF2","SPDEF","LEF1","FOXP2"
)
tf_focus <- unique(toupper(tf_focus))

## =========================
## 3) 读取表达矩阵（第一列为基因名，后面为样本列）
## =========================
read_expr <- function(path){
  x <- fread(path)
  gene <- x[[1]]
  mat <- as.matrix(x[, -1, with=FALSE])
  rownames(mat) <- make.unique(as.character(gene))
  colnames(mat) <- make.unique(colnames(mat))
  mat
}

con_t1 <- read_expr(f_con_t1)
con_t2 <- read_expr(f_con_t2)
sca_t1 <- read_expr(f_sca_t1)
sca_t2 <- read_expr(f_sca_t2)

genes <- Reduce(intersect, list(rownames(con_t1), rownames(con_t2), rownames(sca_t1), rownames(sca_t2)))
con_t1 <- con_t1[genes, , drop=FALSE]
con_t2 <- con_t2[genes, , drop=FALSE]
sca_t1 <- sca_t1[genes, , drop=FALSE]
sca_t2 <- sca_t2[genes, , drop=FALSE]

expr <- cbind(con_t1, con_t2, sca_t1, sca_t2)

## 统一大小写，便于与 tf_focus 匹配
rownames(expr) <- toupper(rownames(expr))

tf_in_expr <- intersect(tf_focus, rownames(expr))
message("TFs requested: ", length(tf_focus))
message("TFs found in expr: ", length(tf_in_expr))
if(length(tf_in_expr) == 0) stop("一个都没匹配到：请检查表达矩阵行名是否为基因symbol（而非ENSG）。")

## =========================
## 4) limma/voom 拟合（全基因拟合，但最后只取指定TF）
## =========================
meta <- data.frame(
  sample = colnames(expr),
  Condition = factor(c(
    rep("CON_T1", ncol(con_t1)),
    rep("CON_T2", ncol(con_t2)),
    rep("SCA_T1", ncol(sca_t1)),
    rep("SCA_T2", ncol(sca_t2))
  ), levels=c("CON_T1","CON_T2","SCA_T1","SCA_T2"))
)

is_counts <- all(abs(expr - round(expr)) < 1e-6, na.rm=TRUE) && max(expr, na.rm=TRUE) > 50
message("Detected counts? ", is_counts)

design <- model.matrix(~0 + meta$Condition)
colnames(design) <- levels(meta$Condition)

if(is_counts){
  dge <- DGEList(expr)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
} else {
  fit <- lmFit(log2(expr + 1), design)
}

contr <- makeContrasts(
  baseline = SCA_T1 - CON_T1,
  levels = design
)

fit2 <- eBayes(contrasts.fit(fit, contr))
res_baseline <- topTable(fit2, coef="baseline", number=Inf, adjust.method="BH", sort.by="P")

## 只保留你指定的 TF
res_baseline$gene <- toupper(rownames(res_baseline))
res_tf <- res_baseline[res_baseline$gene %in% tf_in_expr, , drop=FALSE]
rownames(res_tf) <- res_tf$gene

## =========================
## 5) 火山图：只画指定TF（不保存，只展示）
## =========================
plot_volcano_focus <- function(tt, title, fdr_cut=0.05, lfc_cut=0.5){
  df <- tt
  df$gene <- rownames(df)
  df$neglog10FDR <- -log10(df$adj.P.Val)
  
  df$Sig <- "NS"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC >=  lfc_cut] <- "Up"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC <= -lfc_cut] <- "Down"
  
  ggplot(df, aes(x=logFC, y=neglog10FDR)) +
    geom_point(aes(color=Sig), alpha=0.9, size=3) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr_cut), linetype="dashed") +
    scale_color_manual(values=c(NS="grey70", Up="red3", Down="royalblue3")) +
    labs(title=title, x="log2 Fold Change (SCA_T1 / CON_T1)", y="-log10(FDR)") +
    theme_bw(base_size=12) +
    theme(plot.title = element_text(hjust=0.5)) +
    ggrepel::geom_text_repel(
      aes(label=gene),
      max.overlaps=Inf,
      size=3
    )
}

p <- plot_volcano_focus(
  res_tf,
  title = "TF Volcano (Focused set): SCA_T1 vs CON_T1",
  fdr_cut = 0.05,
  lfc_cut = 0.5
)
print(p)


































library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)

## =========================
## 1) 输入文件路径（GSE168532）
## =========================
f_hd  <- "E:/SCD/其他数据/GSE168532/GSE168532_HD_Gene.csv"
f_scd <- "E:/SCD/其他数据/GSE168532/GSE168532_SCD_Gene.csv"

## =========================
## 2) 你指定的 TF 列表（只画这些）
## =========================
tf_focus <- c(
  "AHR","ATF4","ATF5","BHLHE40","BHLHE41","CLOCK","CREB1","CREM","DBP",
  "EGR1","EGR3","ESR1","FOXO3","HNF4A","JUND","KLF10","KLF9","NFIL3",
  "NFKB2","NFYA","NPAS2","NCOA1","NR1H3","NR5A1","NR5A2","PPARA","PPARG",
  "PROX1","RELB","RORA","RORC","SP1","SREBF1","STAT5A","STAT5B","TEF",
  "TP53","TWIST1","ZFHX3"
)
tf_focus <- unique(toupper(tf_focus))

## =========================
## 3) 读取表达矩阵（默认：第一列是基因名，后面是样本）
## =========================
read_expr <- function(path){
  x <- fread(path)
  gene <- x[[1]]
  mat <- as.matrix(x[, -1, with=FALSE])
  rownames(mat) <- make.unique(as.character(gene))
  colnames(mat) <- make.unique(colnames(mat))
  storage.mode(mat) <- "numeric"
  mat
}

hd  <- read_expr(f_hd)
scd <- read_expr(f_scd)

## 取交集基因并合并
genes <- intersect(rownames(hd), rownames(scd))
hd  <- hd[genes, , drop=FALSE]
scd <- scd[genes, , drop=FALSE]
expr <- cbind(hd, scd)

## 统一大小写，便于与 TF 列表匹配
rownames(expr) <- toupper(rownames(expr))
tf_in_expr <- intersect(tf_focus, rownames(expr))
message("TFs requested: ", length(tf_focus))
message("TFs found in expr: ", length(tf_in_expr))
if(length(tf_in_expr) == 0) stop("一个都没匹配到：请检查表达矩阵行名是否为 gene symbol（而非ENSG等）。")

## =========================
## 4) limma/voom 拟合（两组：HD vs SCD）
## =========================
meta <- data.frame(
  sample = colnames(expr),
  Condition = factor(c(rep("HD", ncol(hd)), rep("SCD", ncol(scd))), levels=c("HD","SCD"))
)

design <- model.matrix(~0 + meta$Condition)
colnames(design) <- levels(meta$Condition)

## 自动判断数据类型：counts / 非负归一化 / 已log(可能含负值)
min_val <- min(expr, na.rm=TRUE)
is_counts <- (min_val >= 0) && all(abs(expr - round(expr)) < 1e-6, na.rm=TRUE) && max(expr, na.rm=TRUE) > 50

message("Min value in expr: ", min_val)
message("Detected counts? ", is_counts)

if(is_counts){
  dge <- DGEList(expr)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
} else if(min_val < 0){
  # 可能已是log/标准化数据（包含负值），直接limma
  fit <- lmFit(expr, design)
} else {
  # 非负但非整数：TPM/FPKM/CPM等，log2(x+1)
  fit <- lmFit(log2(expr + 1), design)
}

contr <- makeContrasts(SCD_vs_HD = SCD - HD, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))
res <- topTable(fit2, coef="SCD_vs_HD", number=Inf, adjust.method="BH", sort.by="P")

## 只保留指定 TF
res$gene <- toupper(rownames(res))
res_tf <- res[res$gene %in% tf_in_expr, , drop=FALSE]
rownames(res_tf) <- res_tf$gene

## =========================
## 5) 火山图：只画指定TF（不保存，只展示）
## =========================
plot_volcano_focus <- function(tt, title, fdr_cut=0.05, lfc_cut=0.5){
  df <- tt
  df$gene <- rownames(df)
  df$neglog10FDR <- -log10(df$adj.P.Val)
  
  df$Sig <- "NS"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC >=  lfc_cut] <- "Up"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC <= -lfc_cut] <- "Down"
  
  ggplot(df, aes(x=logFC, y=neglog10FDR)) +
    geom_point(aes(color=Sig), alpha=0.9, size=3) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr_cut), linetype="dashed") +
    scale_color_manual(values=c(NS="grey70", Up="red3", Down="royalblue3")) +
    labs(title=title, x="log2 Fold Change (SCD / HD)", y="-log10(FDR)") +
    theme_bw(base_size=12) +
    theme(plot.title = element_text(hjust=0.5)) +
    ggrepel::geom_text_repel(aes(label=gene), max.overlaps=Inf, size=3)
}

p <- plot_volcano_focus(
  res_tf,
  title = "TF Volcano (GSE168532): SCD vs HD",
  fdr_cut = 0.05,
  lfc_cut = 0.5
)
print(p)



library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)

## =========================
## 1) 输入文件路径（GSE168532）
## =========================
f_hd  <- "E:/SCD/其他数据/GSE168532/GSE168532_HD_Gene.csv"
f_scd <- "E:/SCD/其他数据/GSE168532/GSE168532_SCD_Gene.csv"

## =========================
## 2) 你指定的 TF 列表（只画这些）
## =========================
tf_focus <- c(
  "CEBPA","CREB1","ELK1","EPAS1","FOSL2","GATA6","GLI1","GLI3","FOXA1",
  "HOXA5","HES1","RBPJ","MSX1","MYCN","NFIB","PROX1","RARA","RARG","SIX1",
  "SOX11","SP3","SREBF1","SRF","TBX2","TCF21","THRA","THRB","WT1","HMGA2",
  "MSC","KLF2","SPDEF","LEF1","FOXP2"
)
tf_focus <- unique(toupper(tf_focus))

## =========================
## 3) 读取表达矩阵（默认：第一列是基因名，后面是样本）
## =========================
read_expr <- function(path){
  x <- fread(path)
  gene <- x[[1]]
  mat <- as.matrix(x[, -1, with=FALSE])
  rownames(mat) <- make.unique(as.character(gene))
  colnames(mat) <- make.unique(colnames(mat))
  storage.mode(mat) <- "numeric"
  mat
}

hd  <- read_expr(f_hd)
scd <- read_expr(f_scd)

## 取交集基因并合并
genes <- intersect(rownames(hd), rownames(scd))
hd  <- hd[genes, , drop=FALSE]
scd <- scd[genes, , drop=FALSE]
expr <- cbind(hd, scd)

## 统一大小写，便于与 TF 列表匹配
rownames(expr) <- toupper(rownames(expr))
tf_in_expr <- intersect(tf_focus, rownames(expr))
message("TFs requested: ", length(tf_focus))
message("TFs found in expr: ", length(tf_in_expr))
if(length(tf_in_expr) == 0) stop("一个都没匹配到：请检查表达矩阵行名是否为 gene symbol（而非ENSG等）。")

## =========================
## 4) limma/voom 拟合（两组：HD vs SCD）
## =========================
meta <- data.frame(
  sample = colnames(expr),
  Condition = factor(c(rep("HD", ncol(hd)), rep("SCD", ncol(scd))), levels=c("HD","SCD"))
)

design <- model.matrix(~0 + meta$Condition)
colnames(design) <- levels(meta$Condition)

## 自动判断数据类型：counts / 非负归一化 / 已log(可能含负值)
min_val <- min(expr, na.rm=TRUE)
is_counts <- (min_val >= 0) && all(abs(expr - round(expr)) < 1e-6, na.rm=TRUE) && max(expr, na.rm=TRUE) > 50

message("Min value in expr: ", min_val)
message("Detected counts? ", is_counts)

if(is_counts){
  dge <- DGEList(expr)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
} else if(min_val < 0){
  fit <- lmFit(expr, design)          # 已log/标准化（可能含负值）
} else {
  fit <- lmFit(log2(expr + 1), design) # TPM/FPKM/CPM等
}

contr <- makeContrasts(SCD_vs_HD = SCD - HD, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))
res <- topTable(fit2, coef="SCD_vs_HD", number=Inf, adjust.method="BH", sort.by="P")

## 只保留指定 TF
res$gene <- toupper(rownames(res))
res_tf <- res[res$gene %in% tf_in_expr, , drop=FALSE]
rownames(res_tf) <- res_tf$gene

## =========================
## 5) 火山图：只画指定TF（不保存，只展示）
## =========================
plot_volcano_focus <- function(tt, title, fdr_cut=0.05, lfc_cut=0.5){
  df <- tt
  df$gene <- rownames(df)
  df$neglog10FDR <- -log10(df$adj.P.Val)
  
  df$Sig <- "NS"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC >=  lfc_cut] <- "Up"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC <= -lfc_cut] <- "Down"
  
  ggplot(df, aes(x=logFC, y=neglog10FDR)) +
    geom_point(aes(color=Sig), alpha=0.9, size=3) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr_cut), linetype="dashed") +
    scale_color_manual(values=c(NS="grey70", Up="red3", Down="royalblue3")) +
    labs(title=title, x="log2 Fold Change (SCD / HD)", y="-log10(FDR)") +
    theme_bw(base_size=12) +
    theme(plot.title = element_text(hjust=0.5)) +
    ggrepel::geom_text_repel(aes(label=gene), max.overlaps=Inf, size=3)
}

p <- plot_volcano_focus(
  res_tf,
  title = "TF Volcano (GSE168532): SCD vs HD",
  fdr_cut = 0.05,
  lfc_cut = 0.5
)
print(p)



library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)

## =========================
## 1) 输入文件路径（GSE168532）
## =========================
f_hd   <- "E:/SCD/其他数据/GSE168532/GSE168532_HD_Gene.csv"
f_scd  <- "E:/SCD/其他数据/GSE168532/GSE168532_SCD_Gene.csv"
tf_file <- "E:/SCD/数据/构建网络的数据/human_tf_list.txt"  # 1列，列名TF

## =========================
## 2) 读取表达矩阵（第一列为基因名，后面为样本列）
## =========================
read_expr <- function(path){
  x <- fread(path)
  gene <- x[[1]]
  mat <- as.matrix(x[, -1, with=FALSE])
  rownames(mat) <- make.unique(as.character(gene))
  colnames(mat) <- make.unique(colnames(mat))
  storage.mode(mat) <- "numeric"
  mat
}

hd  <- read_expr(f_hd)
scd <- read_expr(f_scd)

## 取交集基因并合并
genes <- intersect(rownames(hd), rownames(scd))
hd  <- hd[genes, , drop=FALSE]
scd <- scd[genes, , drop=FALSE]
expr <- cbind(hd, scd)

## 统一基因名大小写，便于匹配
rownames(expr) <- toupper(rownames(expr))

## =========================
## 3) 读取 TF 列表（全部TF）
## =========================
tf_dt <- fread(tf_file)  # 期望列名为 TF
if(!("TF" %in% colnames(tf_dt))) stop("TF文件未检测到列名 'TF'，请确认human_tf_list.txt第一行列名为TF。")

tf_list <- unique(toupper(tf_dt$TF))
tf_in_expr <- intersect(tf_list, rownames(expr))

message("TFs in file: ", length(tf_list))
message("TFs found in expr: ", length(tf_in_expr))
if(length(tf_in_expr) == 0) stop("一个TF都没匹配到：请检查表达矩阵行名是否为gene symbol（而非ENSG等）。")

## =========================
## 4) limma/voom 拟合（两组：HD vs SCD）
## =========================
meta <- data.frame(
  sample = colnames(expr),
  Condition = factor(c(rep("HD", ncol(hd)), rep("SCD", ncol(scd))), levels=c("HD","SCD"))
)

design <- model.matrix(~0 + meta$Condition)
colnames(design) <- levels(meta$Condition)

## 自动判断数据类型：counts / 非负归一化 / 已log(可能含负值)
min_val <- min(expr, na.rm=TRUE)
is_counts <- (min_val >= 0) &&
  all(abs(expr - round(expr)) < 1e-6, na.rm=TRUE) &&
  max(expr, na.rm=TRUE) > 50

message("Min value in expr: ", min_val)
message("Detected counts? ", is_counts)

if(is_counts){
  dge <- DGEList(expr)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
} else if(min_val < 0){
  fit <- lmFit(expr, design)            # 已log/标准化（可能含负值）
} else {
  fit <- lmFit(log2(expr + 1), design)  # TPM/FPKM/CPM等
}

contr <- makeContrasts(SCD_vs_HD = SCD - HD, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))

res <- topTable(fit2, coef="SCD_vs_HD", number=Inf, adjust.method="BH", sort.by="P")

## 只保留 TF
res$gene <- toupper(rownames(res))
res_tf <- res[res$gene %in% tf_in_expr, , drop=FALSE]
rownames(res_tf) <- res_tf$gene

message("TFs in DE result: ", nrow(res_tf))

## =========================
## 5) 火山图：只画TF（不保存，只展示）
## =========================
plot_volcano_tf <- function(tt, title, fdr_cut=0.05, lfc_cut=0.5, label_top=25){
  df <- tt
  df$gene <- rownames(df)
  df$neglog10FDR <- -log10(df$adj.P.Val)
  
  df$Sig <- "NS"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC >=  lfc_cut] <- "Up"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC <= -lfc_cut] <- "Down"
  
  # 标注：优先标注显著TF；如果不够，再补最小FDR的TF
  lab <- df[df$adj.P.Val < fdr_cut & abs(df$logFC) >= lfc_cut, , drop=FALSE]
  lab <- lab[order(lab$adj.P.Val), , drop=FALSE]
  if(nrow(lab) < label_top){
    extra <- df[order(df$adj.P.Val), , drop=FALSE]
    extra <- extra[!extra$gene %in% lab$gene, , drop=FALSE]
    lab <- rbind(lab, head(extra, label_top - nrow(lab)))
  } else {
    lab <- head(lab, label_top)
  }
  
  ggplot(df, aes(x=logFC, y=neglog10FDR)) +
    geom_point(aes(color=Sig), alpha=0.85, size=2.2) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr_cut), linetype="dashed") +
    scale_color_manual(values=c(NS="grey70", Up="red3", Down="royalblue3")) +
    labs(title=title, x="log2 Fold Change (SCD / HD)", y="-log10(FDR)") +
    theme_bw(base_size=12) +
    theme(plot.title = element_text(hjust=0.5)) +
    ggrepel::geom_text_repel(
      data = lab,
      aes(label=gene),
      max.overlaps = Inf,
      size = 3
    )
}

p <- plot_volcano_tf(
  res_tf,
  title = "TF Volcano (GSE168532): SCD vs HD",
  fdr_cut = 0.05,
  lfc_cut = 0.5,
  label_top = 25
)
print(p)



























library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)

## =========================
## 1) 输入文件路径（GSE102881）
## =========================
f_hd  <- "E:/SCD/其他数据/GSE102881/GSE102881_HD_Gene.csv"
f_scd <- "E:/SCD/其他数据/GSE102881/GSE102881_SCD_Gene.csv"

## =========================
## 2) 你指定的 TF 列表（只画这些）
## =========================
tf_focus <- c(
  "AHR","ATF4","ATF5","BHLHE40","BHLHE41","CLOCK","CREB1","CREM","DBP",
  "EGR1","EGR3","ESR1","FOXO3","HNF4A","JUND","KLF10","KLF9","NFIL3",
  "NFKB2","NFYA","NPAS2","NCOA1","NR1H3","NR5A1","NR5A2","PPARA","PPARG",
  "PROX1","RELB","RORA","RORC","SP1","SREBF1","STAT5A","STAT5B","TEF",
  "TP53","TWIST1","ZFHX3"
)
tf_focus <- unique(toupper(tf_focus))

## =========================
## 3) 读取表达矩阵（默认：第一列是基因名，后面是样本）
## =========================
read_expr <- function(path){
  x <- fread(path)
  gene <- x[[1]]
  mat <- as.matrix(x[, -1, with=FALSE])
  rownames(mat) <- make.unique(as.character(gene))
  colnames(mat) <- make.unique(colnames(mat))
  storage.mode(mat) <- "numeric"
  mat
}

hd  <- read_expr(f_hd)
scd <- read_expr(f_scd)

## 取交集基因并合并
genes <- intersect(rownames(hd), rownames(scd))
hd  <- hd[genes, , drop=FALSE]
scd <- scd[genes, , drop=FALSE]
expr <- cbind(hd, scd)

## 统一大小写，便于与 TF 列表匹配
rownames(expr) <- toupper(rownames(expr))
tf_in_expr <- intersect(tf_focus, rownames(expr))
message("TFs requested: ", length(tf_focus))
message("TFs found in expr: ", length(tf_in_expr))
if(length(tf_in_expr) == 0) stop("一个都没匹配到：请检查表达矩阵行名是否为 gene symbol（而非ENSG等）。")

## =========================
## 4) limma/voom 拟合（两组：HD vs SCD）
## =========================
meta <- data.frame(
  sample = colnames(expr),
  Condition = factor(c(rep("HD", ncol(hd)), rep("SCD", ncol(scd))), levels=c("HD","SCD"))
)

design <- model.matrix(~0 + meta$Condition)
colnames(design) <- levels(meta$Condition)

## 自动判断数据类型：counts / 非负归一化 / 已log(可能含负值)
min_val <- min(expr, na.rm=TRUE)
is_counts <- (min_val >= 0) &&
  all(abs(expr - round(expr)) < 1e-6, na.rm=TRUE) &&
  max(expr, na.rm=TRUE) > 50

message("Min value in expr: ", min_val)
message("Detected counts? ", is_counts)

if(is_counts){
  dge <- DGEList(expr)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
} else if(min_val < 0){
  fit <- lmFit(expr, design)           # 已log/标准化（可能含负值）
} else {
  fit <- lmFit(log2(expr + 1), design) # TPM/FPKM/CPM等
}

contr <- makeContrasts(SCD_vs_HD = SCD - HD, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))
res <- topTable(fit2, coef="SCD_vs_HD", number=Inf, adjust.method="BH", sort.by="P")

## 只保留指定 TF
res$gene <- toupper(rownames(res))
res_tf <- res[res$gene %in% tf_in_expr, , drop=FALSE]
rownames(res_tf) <- res_tf$gene

## =========================
## 5) 火山图：只画指定TF（不保存，只展示）
## =========================
plot_volcano_focus <- function(tt, title, fdr_cut=0.05, lfc_cut=0.5){
  df <- tt
  df$gene <- rownames(df)
  df$neglog10FDR <- -log10(df$adj.P.Val)
  
  df$Sig <- "NS"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC >=  lfc_cut] <- "Up"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC <= -lfc_cut] <- "Down"
  
  ggplot(df, aes(x=logFC, y=neglog10FDR)) +
    geom_point(aes(color=Sig), alpha=0.9, size=3) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr_cut), linetype="dashed") +
    scale_color_manual(values=c(NS="grey70", Up="red3", Down="royalblue3")) +
    labs(title=title, x="log2 Fold Change (SCD / HD)", y="-log10(FDR)") +
    theme_bw(base_size=12) +
    theme(plot.title = element_text(hjust=0.5)) +
    ggrepel::geom_text_repel(aes(label=gene), max.overlaps=Inf, size=3)
}

p <- plot_volcano_focus(
  res_tf,
  title = "TF Volcano (GSE102881): SCD vs HD",
  fdr_cut = 0.05,
  lfc_cut = 0.5
)
print(p)





library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)

## =========================
## 1) 输入文件路径（GSE102881）
## =========================
f_hd  <- "E:/SCD/其他数据/GSE102881/GSE102881_HD_Gene.csv"
f_scd <- "E:/SCD/其他数据/GSE102881/GSE102881_SCD_Gene.csv"

## =========================
## 2) 你指定的 TF 列表（只画这些）
## =========================
tf_focus <- c(
  "CEBPA","CREB1","ELK1","EPAS1","FOSL2","GATA6","GLI1","GLI3","FOXA1",
  "HOXA5","HES1","RBPJ","MSX1","MYCN","NFIB","PROX1","RARA","RARG","SIX1",
  "SOX11","SP3","SREBF1","SRF","TBX2","TCF21","THRA","THRB","WT1","HMGA2",
  "MSC","KLF2","SPDEF","LEF1","FOXP2"
)
tf_focus <- unique(toupper(tf_focus))

## =========================
## 3) 读取表达矩阵（默认：第一列是基因名，后面是样本）
## =========================
read_expr <- function(path){
  x <- fread(path)
  gene <- x[[1]]
  mat <- as.matrix(x[, -1, with=FALSE])
  rownames(mat) <- make.unique(as.character(gene))
  colnames(mat) <- make.unique(colnames(mat))
  storage.mode(mat) <- "numeric"
  mat
}

hd  <- read_expr(f_hd)
scd <- read_expr(f_scd)

## 取交集基因并合并
genes <- intersect(rownames(hd), rownames(scd))
hd  <- hd[genes, , drop=FALSE]
scd <- scd[genes, , drop=FALSE]
expr <- cbind(hd, scd)

## 统一大小写，便于与 TF 列表匹配
rownames(expr) <- toupper(rownames(expr))
tf_in_expr <- intersect(tf_focus, rownames(expr))
message("TFs requested: ", length(tf_focus))
message("TFs found in expr: ", length(tf_in_expr))
if(length(tf_in_expr) == 0) stop("一个都没匹配到：请检查表达矩阵行名是否为 gene symbol（而非ENSG等）。")

## =========================
## 4) limma/voom 拟合（两组：HD vs SCD）
## =========================
meta <- data.frame(
  sample = colnames(expr),
  Condition = factor(c(rep("HD", ncol(hd)), rep("SCD", ncol(scd))), levels=c("HD","SCD"))
)

design <- model.matrix(~0 + meta$Condition)
colnames(design) <- levels(meta$Condition)

## 自动判断数据类型：counts / 非负归一化 / 已log(可能含负值)
min_val <- min(expr, na.rm=TRUE)
is_counts <- (min_val >= 0) && all(abs(expr - round(expr)) < 1e-6, na.rm=TRUE) && max(expr, na.rm=TRUE) > 50

message("Min value in expr: ", min_val)
message("Detected counts? ", is_counts)

if(is_counts){
  dge <- DGEList(expr)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
} else if(min_val < 0){
  fit <- lmFit(expr, design)           # 已log/标准化（可能含负值）
} else {
  fit <- lmFit(log2(expr + 1), design) # TPM/FPKM/CPM等
}

contr <- makeContrasts(SCD_vs_HD = SCD - HD, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))
res <- topTable(fit2, coef="SCD_vs_HD", number=Inf, adjust.method="BH", sort.by="P")

## 只保留指定 TF
res$gene <- toupper(rownames(res))
res_tf <- res[res$gene %in% tf_in_expr, , drop=FALSE]
rownames(res_tf) <- res_tf$gene

## =========================
## 5) 火山图：只画指定TF（不保存，只展示）
## =========================
plot_volcano_focus <- function(tt, title, fdr_cut=0.05, lfc_cut=0.5){
  df <- tt
  df$gene <- rownames(df)
  df$neglog10FDR <- -log10(df$adj.P.Val)
  
  df$Sig <- "NS"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC >=  lfc_cut] <- "Up"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC <= -lfc_cut] <- "Down"
  
  ggplot(df, aes(x=logFC, y=neglog10FDR)) +
    geom_point(aes(color=Sig), alpha=0.9, size=3) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr_cut), linetype="dashed") +
    scale_color_manual(values=c(NS="grey70", Up="red3", Down="royalblue3")) +
    labs(title=title, x="log2 Fold Change (SCD / HD)", y="-log10(FDR)") +
    theme_bw(base_size=12) +
    theme(plot.title = element_text(hjust=0.5)) +
    ggrepel::geom_text_repel(aes(label=gene), max.overlaps=Inf, size=3)
}

p <- plot_volcano_focus(
  res_tf,
  title = "TF Volcano (GSE102881): SCD vs HD",
  fdr_cut = 0.05,
  lfc_cut = 0.5
)
print(p)



library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)

## =========================
## 1) 输入文件路径（GSE102881）
## =========================
f_hd   <- "E:/SCD/其他数据/GSE102881/GSE102881_HD_Gene.csv"
f_scd  <- "E:/SCD/其他数据/GSE102881/GSE102881_SCD_Gene.csv"
tf_file <- "E:/SCD/数据/构建网络的数据/human_tf_list.txt"  # 1列，列名TF

## =========================
## 2) 读取表达矩阵（第一列为基因名，后面为样本列）
## =========================
read_expr <- function(path){
  x <- fread(path)
  gene <- x[[1]]
  mat <- as.matrix(x[, -1, with=FALSE])
  rownames(mat) <- make.unique(as.character(gene))
  colnames(mat) <- make.unique(colnames(mat))
  storage.mode(mat) <- "numeric"
  mat
}

hd  <- read_expr(f_hd)
scd <- read_expr(f_scd)

## 取交集基因并合并
genes <- intersect(rownames(hd), rownames(scd))
hd  <- hd[genes, , drop=FALSE]
scd <- scd[genes, , drop=FALSE]
expr <- cbind(hd, scd)

## 统一基因名大小写，便于匹配
rownames(expr) <- toupper(rownames(expr))

## =========================
## 3) 读取 TF 列表（全部TF）
## =========================
tf_dt <- fread(tf_file)  # 期望列名为 TF
if(!("TF" %in% colnames(tf_dt))) stop("TF文件未检测到列名 'TF'，请确认human_tf_list.txt第一行列名为TF。")

tf_list <- unique(toupper(tf_dt$TF))
tf_in_expr <- intersect(tf_list, rownames(expr))

message("TFs in file: ", length(tf_list))
message("TFs found in expr: ", length(tf_in_expr))
if(length(tf_in_expr) == 0) stop("一个TF都没匹配到：请检查表达矩阵行名是否为gene symbol（而非ENSG等）。")

## =========================
## 4) limma/voom 拟合（两组：HD vs SCD）
## =========================
meta <- data.frame(
  sample = colnames(expr),
  Condition = factor(c(rep("HD", ncol(hd)), rep("SCD", ncol(scd))), levels=c("HD","SCD"))
)

design <- model.matrix(~0 + meta$Condition)
colnames(design) <- levels(meta$Condition)

## 自动判断数据类型：counts / 非负归一化 / 已log(可能含负值)
min_val <- min(expr, na.rm=TRUE)
is_counts <- (min_val >= 0) &&
  all(abs(expr - round(expr)) < 1e-6, na.rm=TRUE) &&
  max(expr, na.rm=TRUE) > 50

message("Min value in expr: ", min_val)
message("Detected counts? ", is_counts)

if(is_counts){
  dge <- DGEList(expr)
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
} else if(min_val < 0){
  fit <- lmFit(expr, design)            # 已log/标准化（可能含负值）
} else {
  fit <- lmFit(log2(expr + 1), design)  # TPM/FPKM/CPM等
}

contr <- makeContrasts(SCD_vs_HD = SCD - HD, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))

res <- topTable(fit2, coef="SCD_vs_HD", number=Inf, adjust.method="BH", sort.by="P")

## 只保留 TF
res$gene <- toupper(rownames(res))
res_tf <- res[res$gene %in% tf_in_expr, , drop=FALSE]
rownames(res_tf) <- res_tf$gene

message("TFs in DE result: ", nrow(res_tf))

## =========================
## 5) 火山图：只画TF（不保存，只展示）
## =========================
plot_volcano_tf <- function(tt, title, fdr_cut=0.05, lfc_cut=0.5, label_top=25){
  df <- tt
  df$gene <- rownames(df)
  df$neglog10FDR <- -log10(df$adj.P.Val)
  
  df$Sig <- "NS"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC >=  lfc_cut] <- "Up"
  df$Sig[df$adj.P.Val < fdr_cut & df$logFC <= -lfc_cut] <- "Down"
  
  # 标注：优先标注显著TF；如果不够，再补最小FDR的TF
  lab <- df[df$adj.P.Val < fdr_cut & abs(df$logFC) >= lfc_cut, , drop=FALSE]
  lab <- lab[order(lab$adj.P.Val), , drop=FALSE]
  if(nrow(lab) < label_top){
    extra <- df[order(df$adj.P.Val), , drop=FALSE]
    extra <- extra[!extra$gene %in% lab$gene, , drop=FALSE]
    lab <- rbind(lab, head(extra, label_top - nrow(lab)))
  } else {
    lab <- head(lab, label_top)
  }
  
  ggplot(df, aes(x=logFC, y=neglog10FDR)) +
    geom_point(aes(color=Sig), alpha=0.85, size=2.2) +
    geom_vline(xintercept=c(-lfc_cut, lfc_cut), linetype="dashed") +
    geom_hline(yintercept=-log10(fdr_cut), linetype="dashed") +
    scale_color_manual(values=c(NS="grey70", Up="red3", Down="royalblue3")) +
    labs(title=title, x="log2 Fold Change (SCD / HD)", y="-log10(FDR)") +
    theme_bw(base_size=12) +
    theme(plot.title = element_text(hjust=0.5)) +
    ggrepel::geom_text_repel(
      data = lab,
      aes(label=gene),
      max.overlaps = Inf,
      size = 3
    )
}

p <- plot_volcano_tf(
  res_tf,
  title = "TF Volcano (GSE102881): SCD vs HD",
  fdr_cut = 0.05,
  lfc_cut = 0.5,
  label_top = 25
)
print(p)
