##========================================
## 1. 环境设置 & 读入数据
##========================================
library(data.table)

# 设置工作目录（结果会保存在这里）
setwd("E:/SCD/下游分析/关键转录因子")

# 读入对照组和 SCD 组 TF 度中心性结果
c1 <- fread("two_or_more_overlap_c1_TF_TARGET_poisson_TF.txt")
s1 <- fread("two_or_more_overlap_s1_TF_TARGET_poisson_TF.txt")

#（可选）看一下列名是否是：Gene, Degree, p_value, q_value
print(colnames(c1))
print(colnames(s1))

##========================================
## 2. 按基因合并两组结果
##========================================
tf_merge <- merge(
  c1, s1,
  by       = "Gene",
  suffixes = c("_c1", "_s1"),
  all      = TRUE   # 保留只在一组出现的 TF
)

# 把缺失的 Degree 当作 0（只在另一组出现）
tf_merge[is.na(Degree_c1), Degree_c1 := 0]
tf_merge[is.na(Degree_s1), Degree_s1 := 0]

# 为了后面用 q_value 过滤，把 NA 的 p/q_value 设为 1（即不显著）
na_pq_cols <- c("p_value_c1", "p_value_s1", "q_value_c1", "q_value_s1")
for (col in na_pq_cols) {
  if (!col %in% colnames(tf_merge)) next
  tf_merge[is.na(get(col)), (col) := 1]
}

##========================================
## 3. 计算度差异 & log2 倍数变化
##========================================
# 度差：SCD - 对照
tf_merge[, Degree_diff     := Degree_s1 - Degree_c1]
tf_merge[, Degree_abs_diff := abs(Degree_diff)]

# 度的 log2 fold change（+1 避免除以 0）
tf_merge[, Degree_log2FC := log2((Degree_s1 + 1) / (Degree_c1 + 1))]

##========================================
## 4. 按差异大小排序，找差异较大的 TF
##========================================
# 按绝对差值从大到小排序
tf_sorted <- tf_merge[order(-Degree_abs_diff)]

# 差异最大的前 50 个 TF
tf_top50 <- tf_sorted[1:100]

# 设一个“差异较大”的阈值，比如 |Degree_diff| ≥ 5
tf_bigdiff <- tf_merge[Degree_abs_diff >= 5]

# 同时考虑显著性：差异大且至少一组 q_value < 0.05
tf_bigdiff_sig <- tf_merge[
  Degree_abs_diff >= 5 &
    (q_value_c1 < 0.05 | q_value_s1 < 0.05)
]

# 更偏向 SCD 的 TF：SCD 中度高、SCD 显著、对照不显著
tf_SCD_specific <- tf_merge[
  Degree_abs_diff >= 5 &
    q_value_s1 < 0.05 & q_value_c1 > 0.1
][order(-Degree_diff)]   # 从 SCD >> c1 的排序

##========================================
## 5. 导出结果到文件（方便用 Excel 看）
##========================================
# 1) 全部 TF 对比结果
fwrite(
  tf_merge,
  file = "c1_vs_s1_TF_degree_compare_all.txt",
  sep  = "\t", quote = FALSE
)

# 2) 差异绝对值最大的前 50 个 TF
fwrite(
  tf_top50,
  file = "c1_vs_s1_TF_degree_top50_absdiff.txt",
  sep  = "\t", quote = FALSE
)

# 3) 差异大且至少一组显著的 TF
fwrite(
  tf_bigdiff_sig,
  file = "c1_vs_s1_TF_degree_bigdiff_sig.txt",
  sep  = "\t", quote = FALSE
)

# 4) 更偏向 SCD 的候选关键 TF
fwrite(
  tf_SCD_specific,
  file = "c1_vs_s1_TF_degree_SCD_specific.txt",
  sep  = "\t", quote = FALSE
)

##========================================
## 6. （可选）画散点图看看整体差异
##========================================
# 如果你在 RStudio，可以直接看到图形
png("c1_vs_s1_TF_degree_scatter.png", width = 1200, height = 1200, res = 150)
plot(
  tf_merge$Degree_c1,
  tf_merge$Degree_s1,
  xlab = "Degree (Control, c1)",
  ylab = "Degree (SCD, s1)",
  main = "TF Degree: c1 vs s1"
)
abline(0, 1, lty = 2)
dev.off()





##========================================
## 1. 环境设置 & 读入数据
##========================================
library(data.table)

# 设置工作目录（结果会保存在这里）
setwd("E:/SCD/下游分析/关键转录因子")

# 读入 s1（运动前） 和 s2（运动后） TF 度中心性结果
s1 <- fread("two_or_more_overlap_s1_TF_TARGET_poisson_TF.txt")
s2 <- fread("two_or_more_overlap_s2_TF_TARGET_poisson_TF.txt")

#（可选）看一下列名是否是：Gene, Degree, p_value, q_value
print(colnames(s1))
print(colnames(s2))

##========================================
## 2. 按基因合并两组结果
##========================================
tf_merge <- merge(
  s1, s2,
  by       = "Gene",
  suffixes = c("_s1", "_s2"),   # 改成运动前/运动后
  all      = TRUE
)

# NA 的 Degree 设为 0（只出现在某一组）
tf_merge[is.na(Degree_s1), Degree_s1 := 0]
tf_merge[is.na(Degree_s2), Degree_s2 := 0]

# NA 的 p/q_value 设为 1（表示不显著）
na_pq_cols <- c("p_value_s1", "p_value_s2", "q_value_s1", "q_value_s2")
for (col in na_pq_cols) {
  if (!col %in% colnames(tf_merge)) next
  tf_merge[is.na(get(col)), (col) := 1]
}

##========================================
## 3. 计算度差异 & log2 FC
##========================================
# 度差：运动后 - 运动前
tf_merge[, Degree_diff     := Degree_s2 - Degree_s1]
tf_merge[, Degree_abs_diff := abs(Degree_diff)]

# log2 fold change
tf_merge[, Degree_log2FC := log2((Degree_s2 + 1) / (Degree_s1 + 1))]

##========================================
## 4. 排序和筛选差异 TF
##========================================
# 按绝对差值排序
tf_sorted <- tf_merge[order(-Degree_abs_diff)]

# 差异最大的前 100 个 TF
tf_top100 <- tf_sorted[1:100]

# 设一个“差异大”的阈值，例如 ≥ 5
tf_bigdiff <- tf_merge[Degree_abs_diff >= 5]

# 差异大 + 至少一组显著
tf_bigdiff_sig <- tf_merge[
  Degree_abs_diff >= 5 &
    (q_value_s1 < 0.05 | q_value_s2 < 0.05)
]

# 运动后特异增强 TF（运动后显著，运动前不显著）
tf_s2_specific <- tf_merge[
  Degree_abs_diff >= 5 &
    q_value_s2 < 0.05 & q_value_s1 > 0.1
][order(-Degree_diff)]

##========================================
## 5. 导出结果文件
##========================================
fwrite(tf_merge,
       file = "s1_vs_s2_TF_degree_compare_all.txt",
       sep = "\t", quote = FALSE)

fwrite(tf_top100,
       file = "s1_vs_s2_TF_degree_top100_absdiff.txt",
       sep = "\t", quote = FALSE)

fwrite(tf_bigdiff_sig,
       file = "s1_vs_s2_TF_degree_bigdiff_sig.txt",
       sep = "\t", quote = FALSE)

fwrite(tf_s2_specific,
       file = "s1_vs_s2_TF_degree_s2_specific.txt",
       sep = "\t", quote = FALSE)

##========================================
## 6. 画散点图
##========================================
png("s1_vs_s2_TF_degree_scatter.png", width = 1200, height = 1200, res = 150)
plot(
  tf_merge$Degree_s1,
  tf_merge$Degree_s2,
  xlab = "Degree (before exercise, s1)",
  ylab = "Degree (after exercise, s2)",
  main = "TF Degree: s1 vs s2 (Before vs After Exercise)"
)
abline(0, 1, lty = 2)
dev.off()

















##========================================
## 1. 环境设置 & 读入数据
##========================================
library(data.table)

# 设置工作目录（结果会保存在这里）
setwd("E:/SCD/下游分析/关键转录因子")

# 读入 对照组运动后（c2） 和 患病组运动后（s2）
c2 <- fread("two_or_more_overlap_c2_TF_TARGET_poisson_TF.txt")
s2 <- fread("two_or_more_overlap_s2_TF_TARGET_poisson_TF.txt")

# 看一下列名是否正确
print(colnames(c2))
print(colnames(s2))

##========================================
## 2. 按基因合并两组结果
##========================================
tf_merge <- merge(
  c2, s2,
  by       = "Gene",
  suffixes = c("_c2", "_s2"),   # 修改为运动后两个组
  all      = TRUE
)

# NA 的 Degree 设为 0
tf_merge[is.na(Degree_c2), Degree_c2 := 0]
tf_merge[is.na(Degree_s2), Degree_s2 := 0]

# NA 的 p/q 值设为 1
na_pq_cols <- c("p_value_c2", "p_value_s2", "q_value_c2", "q_value_s2")
for (col in na_pq_cols) {
  if (!col %in% colnames(tf_merge)) next
  tf_merge[is.na(get(col)), (col) := 1]
}

##========================================
## 3. 计算度差异 & log2 fold change
##========================================
# 度差：患病运动后 - 对照运动后
tf_merge[, Degree_diff     := Degree_s2 - Degree_c2]
tf_merge[, Degree_abs_diff := abs(Degree_diff)]

# log2FC
tf_merge[, Degree_log2FC := log2((Degree_s2 + 1) / (Degree_c2 + 1))]

##========================================
## 4. 按差异大小排序，筛选 TF
##========================================
# 排序
tf_sorted <- tf_merge[order(-Degree_abs_diff)]

# 取前 100 个差异最大 TF
tf_top100 <- tf_sorted[1:100]

# 大差异 TF（阈值可调）
tf_bigdiff <- tf_merge[Degree_abs_diff >= 5]

# 差异大 + 至少一组显著
tf_bigdiff_sig <- tf_merge[
  Degree_abs_diff >= 5 &
    (q_value_c2 < 0.05 | q_value_s2 < 0.05)
]

# 患病组运动后特异增强的 TF
tf_s2_specific <- tf_merge[
  Degree_abs_diff >= 5 &
    q_value_s2 < 0.05 & q_value_c2 > 0.1
][order(-Degree_diff)]

##========================================
## 5. 导出结果
##========================================
fwrite(tf_merge,
       file = "c2_vs_s2_TF_degree_compare_all.txt",
       sep = "\t", quote = FALSE)

fwrite(tf_top100,
       file = "c2_vs_s2_TF_degree_top100_absdiff.txt",
       sep = "\t", quote = FALSE)

fwrite(tf_bigdiff_sig,
       file = "c2_vs_s2_TF_degree_bigdiff_sig.txt",
       sep = "\t", quote = FALSE)

fwrite(tf_s2_specific,
       file = "c2_vs_s2_TF_degree_s2_specific.txt",
       sep = "\t", quote = FALSE)

##========================================
## 6. 绘制散点图
##========================================
png("c2_vs_s2_TF_degree_scatter.png", width = 1200, height = 1200, res = 150)
plot(
  tf_merge$Degree_c2,
  tf_merge$Degree_s2,
  xlab = "Degree (Control after exercise, c2)",
  ylab = "Degree (SCD after exercise, s2)",
  main = "TF Degree: c2 vs s2 (Control vs SCD after exercise)"
)
abline(0, 1, lty = 2)
dev.off()
























library(data.table)

#===========================================================
#  函数：比较两组 TF 度中心性（列名：Gene, Degree, Network）
#===========================================================
compare_TF_degree <- function(file_A, file_B, comparison_name) {
  
  #--------------------------------------------
  # 1. 读入两组数据
  #--------------------------------------------
  A <- fread(file_A)
  B <- fread(file_B)
  
  # 显示列名确认
  print(colnames(A))
  print(colnames(B))
  
  #--------------------------------------------
  # 2. 合并数据
  #--------------------------------------------
  tf <- merge(
    A, B,
    by = "Gene",
    suffixes = c("_A", "_B"),
    all = TRUE
  )
  
  # 缺失 Degree → 0（只在某一组出现）
  tf[is.na(Degree_A), Degree_A := 0]
  tf[is.na(Degree_B), Degree_B := 0]
  
  #--------------------------------------------
  # 3. 计算差异指标
  #--------------------------------------------
  tf[, Degree_diff     := Degree_B - Degree_A]
  tf[, Degree_abs_diff := abs(Degree_diff)]
  tf[, Degree_log2FC   := log2((Degree_B + 1) / (Degree_A + 1))]
  
  #--------------------------------------------
  # 4. Top100 差异（分方向）
  #--------------------------------------------
  Increase_top100 <- tf[Degree_diff > 0][order(-Degree_diff)][1:100]
  Decrease_top100 <- tf[Degree_diff < 0][order(Degree_diff)][1:100]
  
  # 大差异：|diff| ≥ 5（可修改）
  tf_bigdiff <- tf[Degree_abs_diff >= 5]
  
  #--------------------------------------------
  # 5. 输出结果文件
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
  # 6. 绘制散点图
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
  
  # 返回列表
  return(list(
    all = tf,
    increase = Increase_top100,
    decrease = Decrease_top100
  ))
}


#===========================================================
#  三组比较（你指定的三个比较）
#===========================================================

# 你的文件路径
file_CT1 <- "E:/SCD/下游分析/度中心/degree_centrality_Control_T1.txt"
file_CT2 <- "E:/SCD/下游分析/度中心/degree_centrality_Control_T2.txt"
file_ST1 <- "E:/SCD/下游分析/度中心/degree_centrality_SCA_T1.txt"
file_ST2 <- "E:/SCD/下游分析/度中心/degree_centrality_SCA_T2.txt"

#-----------------------------------------------------------
# ① Control T1 vs SCA T1
#-----------------------------------------------------------
compare_TF_degree(
  file_A = file_CT1,
  file_B = file_ST1,
  comparison_name = "Control_T1_vs_SCA_T1"
)

#-----------------------------------------------------------
# ② SCA T1 → SCA T2
#-----------------------------------------------------------
compare_TF_degree(
  file_A = file_ST1,
  file_B = file_ST2,
  comparison_name = "SCA_T1_vs_SCA_T2"
)

#-----------------------------------------------------------
# ③ Control T2 → SCA T2
#-----------------------------------------------------------
compare_TF_degree(
  file_A = file_CT2,
  file_B = file_ST2,
  comparison_name = "Control_T2_vs_SCA_T2"
)

