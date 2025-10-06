write.csv(lymph_obj@meta.data[["seurat_clusters"]], file = "seurat_clusters.csv", row.names = TRUE)
write.csv(lymph_obj@meta.data[["GMM"]], file = "gmm_clusters.csv", row.names = TRUE)

# 定义映射规则
mapping <- list(
  c("1", "3"),  # 分化前：1 和 3
  c("5"),       # 分化中：5
  c("6"),       # 分化后1：6
  c("2", "4")   # 分化后2：2 和 4
)

# 创建一个新的向量，用于存储新的分类
branch <- rep(NA, length = nrow(lymph_obj@meta.data))

# 遍历映射规则，将GMM的值映射到新的分类
for (i in seq_along(mapping)) {
  branch[lymph_obj$GMM %in% mapping[[i]]] <- i
}

# 将新的分类添加到Seurat对象的元数据中
lymph_obj <- AddMetaData(lymph_obj, metadata = branch, col.name = "branch")

seurat_obj=lymph_obj
seurat_obj <- SetIdent(seurat_obj, value = "branch")
# 将branch列转换为因子类型
lymph_obj@meta.data$branch <- as.factor(lymph_obj@meta.data$branch)
seurat_obj@meta.data$branch <- as.factor(seurat_obj@meta.data$branch)



# # 获取标记为1或2的细胞的标识符
# cells_to_keep <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$GMM == 1 | seurat_obj@meta.data$GMM == 2]
# 
# # 使用subset函数根据细胞标识符筛选Seurat对象
# subset_data <- subset(seurat_obj, cells = cells_to_keep)
# 
# # 查看筛选后的元数据
# head(subset_data@meta.data)

# 分支感知成对检验 1-2,2-3,2-4
# 创建一个列表来存储差异基因
diff_genes_list <- list()


# 找出在簇标识为1的细胞中相对于簇标识为2的细胞差异表达的基因
markers_1_vs_2 <- FindMarkers(seurat_obj, 
                              ident.1 = 1, 
                              ident.2 = 2, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25,
                              group.by = "branch")

markers_2_vs_3 <- FindMarkers(seurat_obj, 
                              ident.1 = 2, 
                              ident.2 = 3, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25,
                              group.by = "branch")

markers_2_vs_4 <- FindMarkers(seurat_obj, 
                              ident.1 = 2, 
                              ident.2 = 4, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25,
                              group.by = "branch")
# 查看差异表达基因结果
head(markers_1_vs_2)
head(markers_2_vs_3)
head(markers_2_vs_4)

diff_genes_list <- list()
diff_genes_list[['markers_1_vs_2']] <- markers_1_vs_2
diff_genes_list[['markers_2_vs_3']] <- markers_2_vs_3
diff_genes_list[['markers_2_vs_4']] <- markers_2_vs_4

# 查看结果
diff_genes_list


library(pheatmap)

# # 定义一个函数来绘制热图
# plot_heatmap <- function(markers, seurat_obj, title) {
#   # 选择差异基因
#   top_genes <- markers %>% rownames_to_column(var = "gene") %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>% select(gene)
#   
#   # 提取这些基因的表达数据
#   gene_expr <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")[top_genes$gene, ]
#   
#   # 绘制热图
#   pheatmap(gene_expr, clustering_distance_rows = "euclidean",
#            clustering_distance_cols = "euclidean",
#            cluster_rows = TRUE, cluster_cols = TRUE,
#            show_rownames = TRUE, show_colnames = TRUE,
#            color = colorRampPalette(c("blue", "white", "red"))(50),
#            main = title)
# }

gene_expr_data <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")

# 提取特定基因的表达数据
expr_genes_1 <- gene_expr_data[genes_1, ]
expr_genes_2 <- gene_expr_data[genes_2, ]
expr_genes_3 <- gene_expr_data[genes_3, ]
expr_genes_4 <- gene_expr_data[genes_4, ]
# 绘制热图
pheatmap(expr_genes_1, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Expression of Genes in Group 1")

gene_other <- c("UBE2D2","ANAPC16", "LYST", "MAML2", "SIPA1L1", "USP15", "USP3", "CHD3","DOCK8")

# 绘制小提琴图
VlnPlot(lymph_obj, features = gene_other, assay = "RNA", slot = "counts",pt.size = 0.1, ncol = 3, idents = NULL, group.by = "branch")



# 读取数据
data <- read.csv("D:/YYH/gas/Data/GSE60361/mmc3.csv", row.names = 1, check.names = FALSE)
genes_mouse <- read.csv("D:/YYH/gas/Data/GSE60361/all_genes.csv") # 310
genes_mouse <- read.csv("D:/YYH/gas/Data/GSE60361/cluster_genes.csv") # 131
genes_mouse <- read.csv("D:/YYH/gas/Data/GSE60361/DE_genes.csv") # 94
gene_17 <- c('TNC', 'PTPRZ1', 'FAM107A', 'HOPX', 'LIFR', 'IL6ST', 'NOG', 'ITGB5', 'SDC3', 'HS6ST1', 'PRKCA', 'PAX6', 'SOX2', 'EOMES', 'STAT3', 'PTN', 'BMP7')

seurat_mouse <- CreateSeuratObject(counts = data, project = "GSE60361", min.cells = 3, min.features = 200)
HGT_mouse <- run_HGT(GAS = as.matrix(data),result_dir='D:/YYH/gas/Result/RNA_ATAC',data_type = 'scRNA_scATAC')
cell_hgt_matrix_mouse <- HGT_mouse[['cell_hgt_matrix']]
df_mouse <- get_gene_module(seurat_mouse, as.matrix(data), HGT_mouse[["attention"]])
# 检查重复条目
duplicated_rows <- df_mouse[duplicated(df_mouse[, c("node1", "node2")]), ]
print("重复条目：")
print(duplicated_rows)
gas_data_mouse <- acast(df_mouse, node1 ~ node2, value.var = "weight", fill = 0)


D <- as.matrix(genes_mouse)
D <- c(gene_17, as.matrix(genes_mouse))

n=2000
gene_hgt <- getTopGenesIntersection(gas_data_mouse,n,0.3)
gene_hgt_mouse <-gene_hgt[["gene_intersection"]]

# 使用vst方法
seurat_mouse <- NormalizeData(seurat_mouse)
seurat_mouse <- FindVariableFeatures(seurat_mouse, selection.method = "vst", nfeatures = n)
A<-VariableFeatures(seurat_mouse)

# 使用mvp方法
seurat_mouse <- NormalizeData(seurat_mouse)
seurat_mouse <- FindVariableFeatures(seurat_mouse, selection.method = "mvp", nfeatures = n)
B<-VariableFeatures(seurat_mouse)

# 使用disp方法
seurat_mouse <- NormalizeData(seurat_mouse)
seurat_mouse <- FindVariableFeatures(seurat_mouse, selection.method = "disp", nfeatures = n)
C <- VariableFeatures(seurat_mouse)



length(A)
length(B)
length(C)
length(D)
length(gene_hgt_mouse)

# 计算交集
A_D <- intersect(A, D)
B_D <- intersect(B, D)
C_D <- intersect(C, D)
H_D <- intersect(gene_hgt_mouse, D)

# 计算准确率
accuracy_A <- length(A_D) / length(D) * 100
accuracy_B <- length(B_D) / length(D) * 100
accuracy_C <- length(C_D) / length(D) * 100
accuracy_H <- length(H_D) / length(D) * 100
# 创建数据框
accuracy_data <- data.frame(
  Method = c("vst", "mvp", "disp","BranchKGN"),
  Accuracy = c(accuracy_A, accuracy_B, accuracy_C,accuracy_H)
)

# 绘制柱状图
ggplot(accuracy_data, aes(x = Method, y = Accuracy)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = sprintf("%.2f%%", Accuracy)), vjust = -0.5) +
  labs(title = "Accuracy of Different Methods",
       x = "Method",
       y = "Accuracy (%)") +
  theme_minimal()

# 创建一个列表，包含所有集合
sets_list <- list(
  VST = A,
  MVP = B,
  Dispersion = C,
  BranchKGN = gene_hgt_mouse,
  H = D
)

# 绘制交集图
upset(fromList(sets_list), order.by = "freq", nsets = 5, nintersects = 31, 
      mainbar.y.label = "Number of Genes", sets.x.label = "Number of Sets", 
      point.size = 3.5, line.size = 2, set.metadata = NULL, 
      )

# # 计算交集数量
n_A <- length(A)
n_B <- length(B)
n_C <- length(C)
n_D <- length(D)
n_AB <- length(intersect(A, B))
n_AC <- length(intersect(A, C))
n_AD <- length(intersect(A, D))
n_BC <- length(intersect(B, C))
n_BD <- length(intersect(B, D))
n_CD <- length(intersect(C, D))
n_ABC <- length(Reduce(intersect, list(A, B, C)))
n_ABD <- length(Reduce(intersect, list(A, B, D)))
n_ACD <- length(Reduce(intersect, list(A, C, D)))
n_BCD <- length(Reduce(intersect, list(B, C, D)))
n_ABCD <- length(Reduce(intersect,list(A, B, C, D)))

# 关闭当前图形窗口并创建一个新的窗口
grid.newpage()

# 绘制四集合韦恩图
venn.plot <- draw.quad.venn(
  area1 = n_A,
  area2 = n_B,
  area3 = n_C,
  area4 = n_D,
  n12 = n_AB,
  n13 = n_AC,
  n14 = n_AD,
  n23 = n_BC,
  n24 = n_BD,
  n34 = n_CD,
  n123 = n_ABC,
  n124 = n_ABD,
  n134 = n_ACD,
  n234 = n_BCD,
  n1234 = n_ABCD,
  category = c("VST", "MVP", "Dispersion", "Mouse"),
  fill = c("red", "blue", "green", "purple"),
  lty = "blank",
  cex = 1.2,
  cat.cex = 1.2,
  cat.col = c("red", "blue", "green", "purple")
)

