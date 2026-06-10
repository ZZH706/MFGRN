# ============================================================
# Ablation Study Visualization
# 单方法 vs 两两融合 vs 三者融合
# 作者：ChatGPT
# 输出目录：
# E:/SCD/下游分析/消融实验/
# ============================================================

# ============================================================
# 0. 清空环境
# ============================================================

rm(list = ls())

# ============================================================
# 1. 加载R包
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================
# 2. 创建输出目录
# ============================================================

outdir <- "E:/SCD/下游分析/消融实验"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# ============================================================
# 3. 输入数据
# ============================================================

df <- data.frame(
  
  Group = c(
    rep("Control_Pre", 7),
    rep("Control_Post", 7),
    rep("SCA_Pre", 7),
    rep("SCA_Post", 7)
  ),
  
  Method = c(
    
    # -------------------------
    # Control_Pre
    # -------------------------
    "3DCEMA",
    "DeepFGRN",
    "DeepSEM",
    "DS+DF",
    "DS+3D",
    "3D+DF",
    "All_Fusion",
    
    # -------------------------
    # Control_Post
    # -------------------------
    "3DCEMA",
    "DeepFGRN",
    "DeepSEM",
    "DS+DF",
    "DS+3D",
    "3D+DF",
    "All_Fusion",
    
    # -------------------------
    # SCA_Pre
    # -------------------------
    "3DCEMA",
    "DeepFGRN",
    "DeepSEM",
    "DS+DF",
    "DS+3D",
    "3D+DF",
    "All_Fusion",
    
    # -------------------------
    # SCA_Post
    # -------------------------
    "3DCEMA",
    "DeepFGRN",
    "DeepSEM",
    "DS+DF",
    "DS+3D",
    "3D+DF",
    "All_Fusion"
  ),
  
  AVG_AUROC = c(
    
    # =====================================================
    # Control_Pre
    # =====================================================
    
    0.6598,
    0.5543,
    0.5907,
    
    mean(c(0.6573,0.5961,0.5887), na.rm = TRUE),
    
    mean(c(0.6456,0.4851,0.6020), na.rm = TRUE),
    
    0.6218,
    
    0.7608,
    
    # =====================================================
    # Control_Post
    # =====================================================
    
    0.5320,
    0.5891,
    0.5613,
    
    mean(c(0.6035,0.5238,0.5391), na.rm = TRUE),
    
    mean(c(0.7255,0.4999,0.5632), na.rm = TRUE),
    
    mean(c(0.7052,0.5177,0.5561), na.rm = TRUE),
    
    0.6361,
    
    # =====================================================
    # SCA_Pre
    # =====================================================
    
    0.5984,
    0.7097,
    0.5391,
    
    mean(c(0.4928,0.5829,0.5144), na.rm = TRUE),
    
    mean(c(0.6073,0.6091,0.5467), na.rm = TRUE),
    
    mean(c(0.6127,0.6949,0.8128), na.rm = TRUE),
    
    0.7109,
    
    # =====================================================
    # SCA_Post
    # =====================================================
    
    0.6306,
    0.7357,
    0.5742,
    
    mean(c(0.6080,0.5789,0.5029), na.rm = TRUE),
    
    mean(c(0.5195,0.6697,0.5155), na.rm = TRUE),
    
    mean(c(0.5889,0.9824,0.4988), na.rm = TRUE),
    
    0.7851
  )
)

# ============================================================
# 4. 添加类别
# ============================================================

df$Category <- rep(
  c(
    rep("Single", 3),
    rep("Pairwise", 3),
    "Triple"
  ),
  times = 4
)

# ============================================================
# 5. 设置顺序
# ============================================================

df$Method <- factor(
  df$Method,
  levels = c(
    "3DCEMA",
    "DeepFGRN",
    "DeepSEM",
    "DS+DF",
    "DS+3D",
    "3D+DF",
    "All_Fusion"
  )
)

# ============================================================
# 6. 设置组别顺序
# ============================================================

df$Group <- factor(
  df$Group,
  levels = c(
    "Control_Pre",
    "Control_Post",
    "SCA_Pre",
    "SCA_Post"
  )
)

# ============================================================
# 7. 绘制主图（消融实验柱状图）
# ============================================================

p1 <- ggplot(
  df,
  aes(
    x = Method,
    y = AVG_AUROC,
    fill = Category
  )
) +
  
  geom_bar(
    stat = "identity",
    width = 0.75,
    color = "black"
  ) +
  
  facet_wrap(
    ~Group,
    nrow = 2
  ) +
  
  geom_text(
    aes(label = round(AVG_AUROC, 3)),
    vjust = -0.4,
    size = 4
  ) +
  
  scale_fill_manual(
    values = c(
      "Single" = "#8DA0CB",
      "Pairwise" = "#FC8D62",
      "Triple" = "#66C2A5"
    )
  ) +
  
  ylim(0, 0.8) +
  
  labs(
    title = "Ablation Study of GRN Fusion Strategies",
    subtitle = "Single-method vs Pairwise Fusion vs Triple Fusion",
    x = "",
    y = "Average AUROC"
  ) +
  
  theme_bw(base_size = 15) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 11
    ),
    
    axis.text.y = element_text(
      face = "bold"
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 14
    ),
    
    strip.text = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 12
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 13,
      hjust = 0.5
    )
  )

# 显示
print(p1)

# ============================================================
# 8. 保存主图
# ============================================================

ggsave(
  filename = file.path(outdir, "Ablation_Barplot.pdf"),
  plot = p1,
  width = 14,
  height = 8
)

ggsave(
  filename = file.path(outdir, "Ablation_Barplot.png"),
  plot = p1,
  width = 14,
  height = 8,
  dpi = 600
)

# ============================================================
# 9. 构建趋势图数据
# ============================================================

trend_df <- data.frame(
  
  Group = c(
    "Control_Pre",
    "Control_Post",
    "SCA_Pre",
    "SCA_Post"
  ),
  
  Single_Best = c(
    0.6598,
    0.5891,
    0.7097,
    0.7357
  ),
  
  Pairwise_Mean = c(
    
    mean(c(
      mean(c(0.6573,0.5961,0.5887)),
      mean(c(0.6456,0.4851,0.6020)),
      0.6218
    )),
    
    mean(c(
      mean(c(0.6035,0.5238,0.5391)),
      mean(c(0.7255,0.4999,0.5632)),
      mean(c(0.7052,0.5177,0.5561))
    )),
    
    mean(c(
      mean(c(0.4928,0.5829,0.5144)),
      mean(c(0.6073,0.6091,0.5467)),
      mean(c(0.6127,0.6949,0.8128))
    )),
    
    mean(c(
      mean(c(0.6080,0.5789,0.5029)),
      mean(c(0.5195,0.6697,0.5155)),
      mean(c(0.5889,0.9824,0.4988))
    ))
  ),
  
  Triple_Fusion = c(
    0.7608,
    0.6361,
    0.7109,
    0.7851
  )
)

# ============================================================
# 10. 转换成长格式
# ============================================================

trend_long <- trend_df %>%
  pivot_longer(
    cols = -Group,
    names_to = "Fusion_Level",
    values_to = "AUROC"
  )

# ============================================================
# 11. 设置顺序
# ============================================================

trend_long$Fusion_Level <- factor(
  trend_long$Fusion_Level,
  levels = c(
    "Single_Best",
    "Pairwise_Mean",
    "Triple_Fusion"
  )
)

# ============================================================
# 12. 绘制趋势图
# ============================================================

p2 <- ggplot(
  trend_long,
  aes(
    x = Fusion_Level,
    y = AUROC,
    group = Group,
    color = Group
  )
) +
  
  geom_line(
    linewidth = 1.5
  ) +
  
  geom_point(
    size = 4
  ) +
  
  geom_text(
    aes(label = round(AUROC, 3)),
    vjust = -0.7,
    size = 4
  ) +
  
  ylim(0.45, 0.85) +
  
  labs(
    title = "Performance Trend Across Fusion Levels",
    subtitle = "Ablation analysis of fusion depth",
    x = "",
    y = "Average AUROC"
  ) +
  
  theme_bw(base_size = 15) +
  
  theme(
    
    axis.text = element_text(
      face = "bold"
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 14
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 12
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 13,
      hjust = 0.5
    )
  )

# 显示
print(p2)

# ============================================================
# 13. 保存趋势图
# ============================================================

ggsave(
  filename = file.path(outdir, "Fusion_Trend.pdf"),
  plot = p2,
  width = 10,
  height = 7
)

ggsave(
  filename = file.path(outdir, "Fusion_Trend.png"),
  plot = p2,
  width = 10,
  height = 7,
  dpi = 600
)

# ============================================================
# 14. 保存数据表
# ============================================================

write.csv(
  df,
  file = file.path(outdir, "Ablation_Data.csv"),
  row.names = FALSE
)

write.csv(
  trend_long,
  file = file.path(outdir, "Fusion_Trend_Data.csv"),
  row.names = FALSE
)

# ============================================================
# 15. 完成提示
# ============================================================

cat("========================================\n")
cat("消融实验可视化完成！\n")
cat("输出目录：\n")
cat(outdir, "\n")
cat("========================================\n")