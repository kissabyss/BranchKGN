# split_data
# gene_expression
# gas_data_merge


# 提取四类表达量数据
rna_data_merge <- list()

for (i in seq_along(index)) { # 根据索引规则自动处理列拼接
  group <- index[[i]]  # 获取当前组的索引
  print(group)
  
  if (length(group) == 1) {
    # 如果只有一个索引，直接提取并转换为普通矩阵
    rna_data_merge[[as.character(i)]] <- as.matrix(split_data[[group]]@assays[["SCT"]]@data)
  } else {
    # 如果有多个索引，进行列拼接
    group_data <- lapply(group, function(j) {
      as.matrix(split_data[[j]]@assays[["SCT"]]@data)  # 转换为普通矩阵
    })
    
    # 按列合并矩阵
    rna_data_merge[[as.character(i)]] <- do.call(cbind, group_data)
  }
}

for (i in seq_along(index)) {
  group <- index[[i]]
  csv_name <- paste0(i, "_split_rna.csv")
  if (length(group) == 1) {
    expression_matrix <- GetAssayData(object = split_data[[i]], assay = 'RNA', slot = "data")
  } else {
  # 将GAS@Dimnames[[1]]作为行名，并将expression_matrix转换为数据框
    group_data <- lapply(group, function(j) {
      as.matrix(split_data[[j]]@assays[["SCT"]]@data)  # 转换为普通矩阵
    })
    expression_matrix <- do.call(cbind, group_data)
  }
  expression_matrix_df <- as.data.frame(expression_matrix, row.names = rownames(lymph_obj@assays[["RNA"]]))
  # 将数据框写入CSV文件
  fwrite(expression_matrix_df, file = file.path(result_path, csv_name),row.names = TRUE)
}

select_gene_exp_data<-list() # 行是基因，列是细胞
# 遍历每个数据集
for (i in 1:length(rna_data_merge)) {
  print(paste("Processing dataset:", i))
  select_gene_exp_data[[i]] <- rna_data_merge[[i]][gene,] 
}


# 提取同一基因在不同簇上的表达数据
expression_data <- data.frame()
for (i in 1:length(rna_data_merge)) {
  cluster_data <- rna_data_merge[[i]][gene, ]  # 提取第 i 个簇的该基因表达数据
  cluster_df <- data.frame(
    Expression = as.numeric(cluster_data),  # 表达值转为数值型
    Cluster = paste0("Cluster_", i)         # 添加簇编号
  )
  expression_data <- rbind(expression_data, cluster_df)  # 合并到总数据框
}

# 查看整理后的数据框
head(expression_data)

# 绘制同一基因在不同簇上的表达情况
p<-ggplot(expression_data, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(outlier.size = 1, outlier.color = "red", alpha = 0.7) +  # 箱线图
  geom_jitter(color = "black", size = 0.5, alpha = 0.5, width = 0.2) + # 添加抖动点
  labs(
    title = paste("Expression of", gene, "in Different Clusters"),
    x = "Cluster",
    y = "Expression Level"
  ) +
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # 旋转x轴标签
    legend.position = "none"  # 去掉图例
  )

print(p)


##  提取这些基因的表达数据
data_matrix <- as.matrix(rna_data_merge[["2"]])
heatmap_data <- data_matrix[gene, ]
# 加载 pheatmap 包
library(pheatmap)

any(is.na(heatmap_data))       # 检查是否有 NA
any(is.nan(heatmap_data))      # 检查是否有 NaN
any(is.infinite(heatmap_data)) # 检查是否有 Inf

# 打印异常值的数量
cat("NA count:", sum(is.na(heatmap_data)), "\n")
cat("NaN count:", sum(is.nan(heatmap_data)), "\n")
cat("Inf count:", sum(is.infinite(heatmap_data)), "\n")

# 加载包
library(ComplexHeatmap)
library(circlize)  # 用于颜色梯度

# 自定义颜色梯度
heatmap_colors <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# 绘制热图
Heatmap(
  heatmap_data,
  name = "Expression",
  col = heatmap_colors,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_title = "Cells",
  row_title = "Genes",
  heatmap_legend_param = list(title = "Expression Level")
)



# 加载必要的包
library(ComplexHeatmap)
library(circlize)


###  创建一个空列表存储热图
heatmap_list <- list()
gene<-genes$gene_intersection
heatmap_colors <- colorRamp2(c(-5, 0, 5), c( "blue","white", "red"))
# 遍历 rna_data_merge 的四层数据
for (i in seq_along(rna_data_merge)) {
  # 提取当前层数据
  data_matrix <- as.matrix(rna_data_merge[[i]])

  # 提取目标基因表达数据
  heatmap_data <- data_matrix[gene, ]

  # 检查数据维度是否正确
  if (is.null(dim(heatmap_data))) {
    heatmap_data <- matrix(heatmap_data, nrow = 1)  # 如果只有一行，将其转换为矩阵格式
    rownames(heatmap_data) <- gene
  }

  # 创建当前层的热图
  heatmap <- Heatmap(
    heatmap_data,
    name = "Expression",  # 设置为统一的图例名称
    col = heatmap_colors,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = paste("Cluster", i),  # 设置标题为对应的层
    row_title = "Genes",
  )

  # 将热图加入到列表中
  heatmap_list[[i]] <- heatmap
}
# 保存热图到 PDF 文件中，确保图例位置正确
output_path <- file.path(result_path, "heatmaps.pdf")
pdf(output_path, width = 12, height = 10,)  # 设置合适的宽度和高度
draw(
  heatmap_list[[1]] + heatmap_list[[2]] + heatmap_list[[3]] + heatmap_list[[4]],
  heatmap_legend_side = "right",  # 图例放在右侧
  annotation_legend_side = "right",  # 确保注释与图例在同侧
)
# 添加主标题
# grid.text(
#   "Expression heatmap of selected genes in differentiated clusters", 
#   x = 0.5, y = 0.96, gp = gpar(fontsize = 12, fontface = "bold")
# )

dev.off()  # 关闭 PDF 设备



merged_data <- do.call(cbind, rna_data_merge)  # 按列合并所有矩阵

genes_of_interest <- genes$gene_intersection  # 替换为你的目标基因名
selected_data <- merged_data[rownames(merged_data) %in% genes_of_interest, ]


# 数据标准化（按行进行 Z-score 标准化）
scaled_data <- t(scale(t(selected_data)))

# 定义样本分组信息
# 为每个矩阵的列指定一个分组标签
sample_groups <- c(
  rep("Cluster1", ncol(rna_data_merge[['1']])),
  rep("Cluster2", ncol(rna_data_merge[['2']])),
  rep("Cluster3", ncol(rna_data_merge[['3']])),
  rep("Cluster4", ncol(rna_data_merge[['4']]))
)

# 构造分组注释数据框
sample_annotation <- data.frame(Cluster = sample_groups)
rownames(sample_annotation) <- colnames(merged_data)  # 设置样本名为行名

# 绘制热力图，按分组信息显示列注释
pheatmap(
  scaled_data,
  cluster_rows = FALSE,        # 对行（基因）进行聚类
  cluster_cols = FALSE,        # 对列（样本）进行聚类
  annotation_col = sample_annotation,  # 添加列注释（分组信息）
  show_rownames = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("white",  "red","blue"))(50),  # 自定义颜色
  fontsize = 10,              # 字体大小
  main = "Gene Expression Heatmap by Cluster"  # 图标题
)



set.seed(1234) 
# 提取表达数据
gene<-genes[["gene_intersection"]]
# gene_expression_1 <- GetAssayData(object = lymph_obj, assay = 'SCT', slot = "count")
gene_expression_1 <- GetAssayData(object = lymph_obj, assay = 'RNA', slot = "count")
gene_expression<-gene_expression_1[gene,]

sce1 <- SingleCellExperiment(
  assays = list(counts = gene_expression)  # 将 gene_expression 作为 counts 数据
)
pca1<- prcomp(t(gene_expression), scale. = FALSE)
# 提取前两个主成分
rd11 <- pca1$x[, 1:2]
plot(rd11, col = rgb(0,0,0,.5), pch=16, asp = 1)
reducedDims(sce1) <- SimpleList(PCA = rd11)
cl11 <- Mclust(rd11)$classification
colData(sce1)$GMM <- cl11
plot(rd11, col = brewer.pal(9,"Set1")[cl11], pch=8, asp = 0.4)

gmm_classification1 <- colData(sce1)$GMM
sce1 <- slingshot(sce1, clusterLabels = 'GMM', reducedDim = 'PCA')
# 提取 PCA 降维结果和颜色
output_file <- file.path(result_path, "slingshot_class_select-gene11.png")s
# 在 RStudio 窗口中绘制图像
plot(reducedDims(sce1)$PCA, 
     col = brewer.pal(9, 'Set1')[sce1$GMM], 
     pch = 16, 
     asp = 0.8,
     main = "Drawing trajectories using selected gene expression data",
     # xlim = c(-50, 150), ylim = c(-30, 80),  # 调整x轴和y轴范围
     xlab = "PCA 1", 
     ylab = "PCA 2")

# 绘制 Slingshot 轨迹线
slingshot_data1 <- SlingshotDataSet(sce1)
# lines(slingshot_data1, lwd = 2, type = 'lineages', col = 'blue')  # 画 class lineage 轨迹
# legend("bottomright", legend = levels(as.factor(sce1$GMM)), fill = brewer.pal(9,'Set1'), title = "GMM Clusters",cex = 0.8,bty = "n")
legend(x=68,y=20, legend = levels(as.factor(sce1$GMM)), fill = brewer.pal(9,'Set1'), title = "GMM Clusters",cex = 0.8,bty = "n")

# 保存到文件
dev.copy(png, filename = output_file, width = 8, height = 6, units = "in", res = 300)
dev.off()

# 提取 PCA 降维结果
pca_coords <- reducedDims(sce1)$PCA
# write.csv(pca_coords,file.path(result_path,"pca.csv"))

# 提取所有轨迹的拟时数据
pseudotime_list <- slingPseudotime(sce1)  # 获取所有轨迹的拟时
if (is.null(pseudotime_list)) {
  stop("Pseudotime information is not available in the SCE object.")
}

# 为每条轨迹绘制时间，使用不同颜色表示
num_trajectories <- ncol(pseudotime_list)  # 轨迹数量
num_trajectories
colors <- brewer.pal(min(9, num_trajectories), "Set1")[1:num_trajectories]   # 分配轨迹颜色
names(colors) <- colnames(pseudotime_list)  # 轨迹颜色标记
# 设置输出文件路径
output_file <- file.path(result_path, "slingshot_pseudotime_new1.png")

# 在 RStudio 窗口中绘制图像
plot(pca_coords, 
     col = rgb(0, 0, 0, 0.5),  # 初始点为透明黑色
     pch = 8, 
     asp = 0.4,
     xlim = c(-50, 25), ylim = c(-45, 15),  # 调整x轴和y轴范围
     main = "Cell Pseudotime in PCA Space",
     xlab = " ", 
     ylab = " ",
     cex.main = 1.2,  # 主标题字体大小
     cex.lab = 1.1,   # 轴标签字体大小
     cex.axis = 0.9)  # 轴刻度字体大小

# 自定义从绿色到蓝色的渐变
green_to_blue <- colorRampPalette(c("green", "blue"))(100)

for (i in seq_len(num_trajectories)) {
  pseudotime <- pseudotime_list[, i]
  
  # 处理拟时的颜色梯度
  pseudotime_scaled <- as.numeric(cut(pseudotime, breaks = 100))
  cell_colors <- green_to_blue[pseudotime_scaled]
  
  # 在 PCA 图上绘制点，颜色对应拟时
  points(pca_coords, 
         col = cell_colors, 
         pch = 16, cex = 0.8)
}


# 遍历每条轨迹，分轨迹绘制拟时
for (i in seq_len(num_trajectories)) {
  pseudotime <- pseudotime_list[, i]
  
  # 处理拟时的颜色梯度
  pseudotime_colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
  pseudotime_scaled <- as.numeric(cut(pseudotime, breaks = 100))
  cell_colors <- pseudotime_colors[pseudotime_scaled]
  
  # 在 PCA 图上绘制点，颜色对应拟时
  points(pca_coords, 
         col = cell_colors, 
         pch = 16,cex = 0.8)
}
# 添加图例
# 使用 legend 函数添加水平图例
legend(x = 30, y = 50,  # 图例的坐标
       legend = c("Low", "High"),
       fill = c(pseudotime_colors[1], pseudotime_colors[100]),
       # col = brewer.pal(10, "YlOrRd"),
       col = green_to_blue,
       horiz = FALSE,
       cex = 0.8,
       bty = "n",
       title = "Pseudotime")

# 添加连续颜色图例
# 使用 image.plot 添加连续颜色图例
# image.plot(legend.only = TRUE, 
#            zlim = range(pseudotime_list, na.rm = TRUE), 
#            col = pseudotime_colors, 
#            horizontal = FALSE,    # 设置图例为水平显示
#            legend.lab = "Pseudotime", 
#            legend.line = 10,      # 调整图例位置
#            legend.cex = 0.8,      # 调整图例文字大小
#            legend.width = 0.8,    # 调整图例宽度
#            legend.shrink = 0.9,   # 缩小图例比例
#            add = TRUE)  # 将图例添加到现有图形中

# 保存到文件
dev.copy(png, filename = output_file, width = 8, height = 6, units = "in", res = 300)
dev.off()
cat("Visualization of all trajectories' pseudotime saved to", output_file, "\n")


## 计算所选基因的L2FC，p_value,FDR
gene_list <- genes[["gene_intersection"]]
# 提取每个时期的基因表达矩阵
expr_t1 <- as.matrix(rna_data_merge[[1]])  # 第一个时期
expr_t2 <- as.matrix(rna_data_merge[[2]])  # 第二个时期
expr_t3 <- as.matrix(rna_data_merge[[3]])  # 第三个时期
expr_t4 <- as.matrix(rna_data_merge[[4]])  # 第四个时期
# 将所有时期的表达矩阵存储到列表中
expr_list <- list(t1 = expr_t1, t2 = expr_t2, t3 = expr_t3, t4 = expr_t4)
# 创建一个空数据框保存结果
results <- data.frame(
  Gene = character(),
  Comparison = character(),
  L2FC = numeric(),
  p_value = numeric(),
  FDR = numeric(),
  stringsAsFactors = FALSE
)

# 获取所有两两组合
comparisons <- combn(names(expr_list), 2, simplify = FALSE)
# 遍历每对比较
for (comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  
  # 提取两个组的表达矩阵
  expr1 <- expr_list[[group1]]
  expr2 <- expr_list[[group2]]
  
  # 遍历目标基因
  for (gene in gene_list) {
    # 检查基因是否存在于所有时期的数据中
    if (gene %in% rownames(expr1) & gene %in% rownames(expr2)) {
      # 获取基因在两个时期的表达值
      values1 <- expr1[gene, ]
      values2 <- expr2[gene, ]
      
      # 计算 L2FC (log2 fold change)
      mean_expr1 <- mean(values1, na.rm = TRUE)
      mean_expr2 <- mean(values2, na.rm = TRUE)
      l2fc <- log2((mean_expr2 + 1e-6) / (mean_expr1 + 1e-6))  # 添加伪计数
      
      # 计算 p-value (使用 t 检验)
      p_value <- t.test(values1, values2)$p.value
      
      # 保存结果到数据框
      results <- rbind(
        results,
        data.frame(
          Gene = gene,
          Comparison = paste(group1, "vs", group2),
          L2FC = l2fc,
          p_value = p_value,
          FDR = NA  # FDR 暂时置为空
        )
      )
    }
  }
}
# 计算 FDR (使用 Benjamini-Hochberg 校正)
results$FDR <- p.adjust(results$p_value, method = "BH")

# 筛选 Comparison 包含 "t2" 且 p_value 显著的基因
significant_genes <- results %>%
  filter(
    grepl("t2", Comparison),  # 筛选 Comparison 中包含 "t2"
    p_value < 0.05            # 筛选 p_value 小于 0.05 的基因
  )

# 打印筛选后的结果
print(significant_genes)

library(limma)

# 1. 创建分组信息
group <- factor(c(rep("t1", ncol(expr_t1)), 
                  rep("t2", ncol(expr_t2)), 
                  rep("t3", ncol(expr_t3)), 
                  rep("t4", ncol(expr_t4))))

# 2. 合并所有时期的表达矩阵
expr_all <- cbind(expr_t1, expr_t2, expr_t3, expr_t4)

# 3. 创建设计矩阵
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 4. 差异表达分析
fit <- lmFit(expr_all, design)

# 5. 创建对比矩阵（例如，每个时期相对于 t2）
contrast_matrix <- makeContrasts(
  t1_vs_t2 = t1 - t2,
  t3_vs_t2 = t3 - t2,
  t4_vs_t2 = t4 - t2,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 6. 提取结果
# 针对 gene_list 筛选基因
results_t1 <- topTable(fit2, coef = "t1_vs_t2", number = Inf)
results_t3 <- topTable(fit2, coef = "t3_vs_t2", number = Inf)
results_t4 <- topTable(fit2, coef = "t4_vs_t2", number = Inf)

# 筛选 gene_list 中的基因
results_t1 <- results_t1[rownames(results_t1) %in% gene_list, ]
results_t3 <- results_t3[rownames(results_t3) %in% gene_list, ]
results_t4 <- results_t4[rownames(results_t4) %in% gene_list, ]

# 7. 合并结果
all_results <- list(
  t1 = results_t1,
  t3 = results_t3,
  t4 = results_t4
)

# 显示结果
all_results

library(ggplot2)

results_t1 <- results_t3 %>%
  mutate(significant = ifelse(abs(logFC) > 0.5 & P.Value < 0.1, "Significant", "Non-significant"))

# ggplot(results_t1, aes(x = logFC, y = -log10(P.Value), color = significant)) +
#   geom_point(alpha = 0.6) +
#   geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "red") +
#   geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "red") +
#   labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(P-value)") +
#   scale_color_manual(values = c("grey", "red")) +
#   theme_minimal()

gene_matrix <- expr_all[rownames(expr_all) %in% gene_list, ]

pheatmap(gene_matrix, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("blue", "white", "red"))(20),
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Heatmap of Differential Expression")

# 绘制热图，禁用聚类
pheatmap(gene_matrix, 
         clustering_distance_rows = "euclidean",  # 设置行聚类距离度量
         clustering_distance_cols = "euclidean",  # 设置列聚类距离度量
         cluster_rows = FALSE,  # 禁用行聚类
         cluster_cols = FALSE,  # 禁用列聚类
         color = colorRampPalette(c("blue", "white", "red"))(20),  # 颜色梯度
         show_rownames = FALSE,  # 不显示行名
         show_colnames = FALSE,  # 不显示列名
         main = "Heatmap of Differential Expression")  # 主标题

# 提取 logFC 数据
logFC_matrix <- sapply(all_results, function(x) x$logFC)
# 对 logFC_matrix 进行 Z-score 归一化
logFC_matrix_normalized <- scale(logFC_matrix)
rownames(logFC_matrix_normalized) <- rownames(results_t1)
library(pheatmap)

# 绘制热图
pheatmap(logFC_matrix_normalized, 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         main = "Log2 Fold Change of Gene Set Across Time Points",
         color = colorRampPalette(c("blue", "white", "red"))(50),  # 自定义颜色映射
         show_rownames = FALSE,  # 显示基因名
         show_colnames = TRUE,  # 显示时间点
         annotation_col = NULL,  # 如果需要，可以添加列注释
         annotation_row = NULL   # 如果需要，可以添加行注释
)
 

