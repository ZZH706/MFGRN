###############################################################
#  小世界网络判定（Humphries σ 指数 + 网络密度 Density）
# ------------------------------------------------------------
#  批量分析：刚才的四个 overlapped 网络文件
#  输入 : TF Target EdgeWeight（.txt, 制表符分隔）
#  输出 : network_topology_sigma_<文件名>.txt
#  σ = (C_real / ⟨C_rand⟩) / (L_real / ⟨L_rand⟩)，σ > 1 判小世界
###############################################################

# -------------------- 依赖包 --------------------
library(igraph)
library(readr)
library(dplyr)
library(tools)

# -------------------- 参数 ----------------------
edge_files <- c(
  "E:/1数据/≥2/c1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/c2_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s1_TF_TARGET_cleaned.txt",
  "E:/1数据/≥2/s2_TF_TARGET_cleaned.txt"
)
N_rand    <- 100       # 随机网络次数
set.seed(123)         # 结果可复现

# -------------------- 循环处理每个网络 --------------------
for (edge_file in edge_files) {
  
  cat("\n==============================\n")
  cat("📂 读取:", edge_file, "\n")
  
  # 读取完整的 TF–Target–EdgeWeight 表，不做筛选
  edges <- read_tsv(edge_file,
                    col_types = cols(
                      TF         = col_character(),
                      Target     = col_character(),
                      EdgeWeight = col_double()
                    ))
  stopifnot(nrow(edges) > 0)
  
  # 构建有向图并简化，再转成无向
  g_dir <- graph_from_data_frame(edges[, c("TF","Target")],
                                 directed = TRUE) %>%
    simplify(remove.multiple = TRUE,
             remove.loops     = TRUE)
  g <- as.undirected(g_dir, mode = "collapse")
  
  cat("节点数:", vcount(g), "   边数:", ecount(g), "\n")
  
  # 真实网络指标
  C_real <- transitivity(g, "average")
  L_real <- average.path.length(g, unconnected = TRUE)
  D_real <- edge_density(g, loops = FALSE)
  
  cat(sprintf("C = %.6f   L = %.6f   Density = %.6f\n",
              C_real, L_real, D_real))
  
  # 随机网络采样
  get_stats <- function(nv, ne) {
    gr <- sample_gnm(nv, ne, directed = FALSE, loops = FALSE)
    c(C = transitivity(gr, "average"),
      L = average.path.length(gr, unconnected = TRUE))
  }
  cat("⏳ 生成随机网络", N_rand, "次...\n")
  rand_stats <- replicate(N_rand,
                          get_stats(vcount(g), ecount(g)))
  
  C_rand_mu <- mean(rand_stats["C", ])
  L_rand_mu <- mean(rand_stats["L", ])
  cat(sprintf("⟨C_rand⟩ = %.6f   ⟨L_rand⟩ = %.6f\n",
              C_rand_mu, L_rand_mu))
  
  # 计算 σ 指数与判定
  sigma <- (C_real / C_rand_mu) / (L_real / L_rand_mu)
  is_sw <- sigma > 1
  cat(sprintf("σ = %.3f  → %s小世界网络\n",
              sigma, ifelse(is_sw, "符合", "不符合")))
  
  # 保存结果
  out <- tibble(
    网络        = file_path_sans_ext(basename(edge_file)),
    节点数      = vcount(g),
    边数        = ecount(g),
    Density     = round(D_real, 6),
    C_real      = round(C_real, 6),
    L_real      = round(L_real, 6),
    C_rand_mean = round(C_rand_mu, 6),
    L_rand_mean = round(L_rand_mu, 6),
    sigma       = round(sigma, 3),
    小世界判定   = ifelse(is_sw, "是", "否"),
    随机次数     = N_rand
  )
  
  out_file <- paste0("network_topology_sigma_",
                     file_path_sans_ext(basename(edge_file)),
                     ".txt")
  write.table(out,
              file      = out_file,
              sep       = "\t",
              quote     = FALSE,
              row.names = FALSE,
              fileEncoding = "UTF-8")
  cat("✅ 已保存:", out_file, "\n")
}
