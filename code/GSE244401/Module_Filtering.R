library(data.table)
library(openxlsx)
library(miRspongeR)

set.seed(123)

## =========================
## 0. 参数与输出目录
## =========================
min_size <- 10
max_size <- 300

out_dir <- "E:/SCD/下游分析/功能模块/模块"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## =========================
## 1. 读取四个网络
## =========================
net_files <- list(
  c1 = "E:/SCD/数据/合并后的数据/c1_TF_TARGET_cleaned.txt",
  c2 = "E:/SCD/数据/合并后的数据/c2_TF_TARGET_cleaned.txt",
  s1 = "E:/SCD/数据/合并后的数据/s1_TF_TARGET_cleaned.txt",
  s2 = "E:/SCD/数据/合并后的数据/s2_TF_TARGET_cleaned.txt"
)

# 兼容不同列名/列数：至少保证前两列是边的两端；第三列(如有)作为 weight
normalize_linklist <- function(dt) {
  if (!is.data.table(dt)) dt <- as.data.table(dt)
  if (ncol(dt) < 2) stop("网络边文件至少需要两列（TF/Source, Target）")
  
  # 统一列名（不改变内容）
  if (ncol(dt) == 2) {
    setnames(dt, names(dt)[1:2], c("Source", "Target"))
  } else {
    setnames(dt, names(dt)[1:3], c("Source", "Target", "Weight"))
  }
  dt
}

linklists <- lapply(net_files, function(fp) {
  if (!file.exists(fp)) stop("网络文件不存在：", fp)
  normalize_linklist(fread(fp))
})

## =========================
## 2. 仅用 infomap 做模块识别 + 过滤模块大小 10~300
## =========================
run_infomap_and_filter <- function(linklist_dt, min_size = 10, max_size = 300) {
  # netModule 的 modulesize 是“下限”；上限用后处理过滤实现
  clu <- netModule(linklist_dt, method = "infomap", modulesize = min_size)
  
  if (is.null(clu) || length(clu) == 0) return(clu)
  
  # 给模块命名（如果没有名字）
  if (is.null(names(clu)) || any(names(clu) == "")) {
    names(clu) <- paste0("M", seq_along(clu))
  }
  
  # 过滤上限
  sizes <- sapply(clu, function(x) length(unique(x)))
  clu <- clu[sizes >= min_size & sizes <= max_size]
  
  # 重新编号/命名（可选：更整齐）
  if (length(clu) > 0) names(clu) <- paste0("M", seq_along(clu))
  clu
}

clusters <- lapply(linklists, run_infomap_and_filter, min_size = min_size, max_size = max_size)

# 保存模块对象（可选）
save(clusters, file = file.path(out_dir, sprintf("Infomap_modules_%s_%s_all.Rdata", min_size, max_size)))

cat("\nInfomap 模块数量（过滤 10~300 后）：\n")
for (nm in names(clusters)) {
  cat(sprintf("%-3s : %d 个模块\n", nm, length(clusters[[nm]])))
}

## =========================
## 3. 读取四个关键转录因子列表
## =========================
tf_files <- list(
  c1 = "E:/SCD/下游分析/关键转录因子/two_or_more_overlap_c1_TF_TARGET_poisson_TF_Gene_only.txt",
  c2 = "E:/SCD/下游分析/关键转录因子/two_or_more_overlap_c2_TF_TARGET_poisson_TF_Gene_only.txt",
  s1 = "E:/SCD/下游分析/关键转录因子/two_or_more_overlap_s1_TF_TARGET_poisson_TF_Gene_only.txt",
  s2 = "E:/SCD/下游分析/关键转录因子/two_or_more_overlap_s2_TF_TARGET_poisson_TF_Gene_only.txt"
)

read_tf_list <- function(fp) {
  if (!file.exists(fp)) stop("TF 列表文件不存在：", fp)
  x <- readLines(fp, warn = FALSE, encoding = "UTF-8")
  x <- trimws(x)
  x <- x[x != ""]
  unique(x)
}

tf_lists <- lapply(tf_files, read_tf_list)

cat("\n关键 TF 列表大小：\n")
for (nm in names(tf_lists)) {
  cat(sprintf("%-3s : %d 个 TF\n", nm, length(tf_lists[[nm]])))
}

## =========================
## 4. 筛选“前十个包含关键 TF 的模块”
##    排序：KeyTF_Count（降序）-> ModuleSize（降序）-> ModuleID（升序）
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
## 5. 对 c1/c2/s1/s2 各自筛选并导出 Excel
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
  
  cat(sprintf("\n[%s] 含关键TF的模块数（过滤后）：%d；导出Top10：%d\n",
              cond,
              ifelse(is.null(clusters[[cond]]), 0, length(clusters[[cond]])),
              nrow(res$summary)))
}

out_xlsx <- file.path(out_dir, sprintf("Top10_modules_with_keyTF_infomap_%s_%s.xlsx", min_size, max_size))
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("\n全部完成！结果已导出到：\n", out_xlsx, "\n")
