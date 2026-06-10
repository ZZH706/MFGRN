## 安装并载入包（只需安装一次）
# install.packages("fmsb")
library(fmsb)

##-----------------------------
## 1. 录入所有数据
##-----------------------------
# 为了统一刻度，先把所有组的数据放在一起求最大/最小
all_C <- c(0.2751, 0.2295, 0.2938, 0.4227,  # 预测：聚类系数
           0.0020, 0.0020, 0.0022, 0.0021) # 随机：聚类系数
all_L <- c(2.3400, 2.3890, 2.1763, 2.2344, # 预测：平均路径长度
           2.99154, 2.9242, 2.8975, 2.9017)   # 随机：平均路径长度
all_D <- c(0.0020, 0.0020, 0.0022, 0.0021, # 预测：网络密度
           0.0020, 0.0020, 0.0022, 0.0021) # 随机：网络密度

# 统一的最大/最小行（第一行为最大，第二行为最小）
maxmin <- data.frame(
  聚类系数     = c(max(all_C) * 1.1, min(all_C) * 0.9),
  平均路径长度 = c(max(all_L) * 1.1, min(all_L) * 0.9),
  网络密度     = c(max(all_D) * 1.1, min(all_D) * 0.9)
)
rownames(maxmin) <- c("max", "min")

##-----------------------------
## 2. 为每一组构造数据框
##-----------------------------
# 对照组运动前
dat_ctrl_pre <- rbind(
  maxmin,
  `预测网络` = c(0.2751, 2.3400, 0.0020),
  `随机网络` = c(0.0020, 2.99154, 0.0020)
)

# 对照组运动后
dat_ctrl_post <- rbind(
  maxmin,
  `预测网络` = c(0.2295, 2.3890, 0.0020),
  `随机网络` = c(0.0020, 2.9242, 0.0020)
)

# 患病组运动前
dat_dis_pre <- rbind(
  maxmin,
  `预测网络` = c(0.2938, 2.1763, 0.0022),
  `随机网络` = c(0.0022, 2.8975, 0.0022)
)

# 患病组运动后
dat_dis_post <- rbind(
  maxmin,
  `预测网络` = c(0.4227, 2.2344, 0.0021),
  `随机网络` = c(0.0021, 2.9017, 0.0021)
)

##-----------------------------
## 3. 绘制 4 个三角雷达图（去掉百分比刻度标签）
##-----------------------------
par(mfrow = c(2, 2), mar = c(1, 2, 3, 2))  # 2×2 排版

radarchart(dat_ctrl_pre,
           axistype = 1,
           caxislabels = rep("", 5),    # 不显示百分比刻度
           pcol = c("red", "blue"),
           pfcol = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)),
           plwd = 2,
           title = "对照组运动前")
legend("topright", legend = c("预测网络", "随机网络"),
       col = c("red","blue"), pch = 19, bty = "n")

radarchart(dat_ctrl_post,
           axistype = 1,
           caxislabels = rep("", 5),
           pcol = c("red", "blue"),
           pfcol = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)),
           plwd = 2,
           title = "对照组运动后")
legend("topright", legend = c("预测网络", "随机网络"),
       col = c("red","blue"), pch = 19, bty = "n")

radarchart(dat_dis_pre,
           axistype = 1,
           caxislabels = rep("", 5),
           pcol = c("red", "blue"),
           pfcol = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)),
           plwd = 2,
           title = "患病组运动前")
legend("topright", legend = c("预测网络", "随机网络"),
       col = c("red","blue"), pch = 19, bty = "n")

radarchart(dat_dis_post,
           axistype = 1,
           caxislabels = rep("", 5),
           pcol = c("red", "blue"),
           pfcol = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)),
           plwd = 2,
           title = "患病组运动后")
legend("topright", legend = c("预测网络", "随机网络"),
       col = c("red","blue"), pch = 19, bty = "n")




## 安装并载入包（只需安装一次）
# install.packages("fmsb")
library(fmsb)

set.seed(123)

##-----------------------------
## 1. 录入四组“预测网络”指标（不再用随机网络）
##-----------------------------
# 三个指标：聚类系数 / 平均路径长度 / 网络密度
ctrl_pre  <- c(0.2751, 2.3400, 0.0020)  # 对照组运动前
ctrl_post <- c(0.2295, 2.3890, 0.0020)  # 对照组运动后
dis_pre   <- c(0.2938, 2.1763, 0.0022)  # 患病组运动前
dis_post  <- c(0.4227, 2.2344, 0.0021)  # 患病组运动后

all_mat <- rbind(ctrl_pre, ctrl_post, dis_pre, dis_post)
colnames(all_mat) <- c("聚类系数", "平均路径长度", "网络密度")

## 统一的最大/最小行（第一行为最大，第二行为最小）
maxmin <- data.frame(
  聚类系数     = c(max(all_mat[, "聚类系数"]) * 1.1,     min(all_mat[, "聚类系数"]) * 0.9),
  平均路径长度 = c(max(all_mat[, "平均路径长度"]) * 1.1, min(all_mat[, "平均路径长度"]) * 0.9),
  网络密度     = c(max(all_mat[, "网络密度"]) * 1.1,     min(all_mat[, "网络密度"]) * 0.9)
)
rownames(maxmin) <- c("max", "min")

##-----------------------------
## 2. 构造三组对比数据框（每张图两条线：A vs B）
##-----------------------------
dat_cmp1 <- rbind(
  maxmin,
  `对照组运动前` = ctrl_pre,
  `患病组运动前` = dis_pre
)

dat_cmp2 <- rbind(
  maxmin,
  `对照组运动后` = ctrl_post,
  `患病组运动后` = dis_post
)

dat_cmp3 <- rbind(
  maxmin,
  `患病组运动前` = dis_pre,
  `患病组运动后` = dis_post
)

##-----------------------------
## 3. 绘制 3 个雷达图（不显示刻度标签）
##-----------------------------
par(mfrow = c(1, 3), mar = c(1, 2, 3, 2))  # 1×3 排版

# 统一绘图参数
axis_n <- 5
cols_line <- c("red", "blue")
cols_fill <- c(rgb(1, 0, 0, 0.25), rgb(0, 0, 1, 0.25))

radarchart(dat_cmp1,
           axistype = 1,
           caxislabels = rep("", axis_n),  # 不显示刻度标签
           pcol = cols_line,
           pfcol = cols_fill,
           plwd = 2,
           title = "对照组运动前 vs 患病组运动前")
legend("topright", legend = c("对照组运动前", "患病组运动前"),
       col = cols_line, pch = 19, bty = "n", cex = 0.9)

radarchart(dat_cmp2,
           axistype = 1,
           caxislabels = rep("", axis_n),
           pcol = cols_line,
           pfcol = cols_fill,
           plwd = 2,
           title = "对照组运动后 vs 患病组运动后")
legend("topright", legend = c("对照组运动后", "患病组运动后"),
       col = cols_line, pch = 19, bty = "n", cex = 0.9)

radarchart(dat_cmp3,
           axistype = 1,
           caxislabels = rep("", axis_n),
           pcol = cols_line,
           pfcol = cols_fill,
           plwd = 2,
           title = "患病组运动前 vs 患病组运动后")
legend("topright", legend = c("患病组运动前", "患病组运动后"),
       col = cols_line, pch = 19, bty = "n", cex = 0.9)
