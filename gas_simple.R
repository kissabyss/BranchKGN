# Load function
set.seed(12345) ## for reproducibility
source("D:\\YYH\\gas\\gas.r")

# set python environment
Sys.setenv(RETICULATE_PYTHON = "D:\\Software\\Anaconda3\\envs\\deepmaps_env\\python.exe")
use_python("D:\\Software\\Anaconda3\\envs\\deepmaps_env\\python.exe")
py_config()


# file_name='pbmc_granulocyte_sorted_3k'
# file_name='pbmc_unsorted_10k'
# file_name='pbmc_unsorted_3k'
file_name='e18_mouse_brain_fresh_5k'
result_path=file.path(getwd(),'Result',file_name)
if (!file.exists(result_path)) {
  dir.create(result_path, recursive = TRUE)
}

# read data
# h5Path = "lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5"
# h5Path = "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
# h5Path = "pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5"
# h5Path = "pbmc_unsorted_10k_filtered_feature_bc_matrix.h5"
# h5Path = "pbmc_unsorted_3k_filtered_feature_bc_matrix.h5"
h5Path = "e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5"
# 创建 a seurat object
lymph_obj <- ReadData(file_name = file_name, min_cell = 0.01)
VlnPlot(lymph_obj, features = c("nCount_RNA", "nFeature_RNA",  "nCount_ATAC", "nFeature_ATAC"), ncol = 2,
        log = TRUE, pt.size = 0) + NoLegend()

# 聚类分析，wnn构建矩阵
lymph_obj<-performAnalysis(lymph_obj)


# 可视化 RNA UMAP 聚类
umap_rna <- DimPlot(lymph_obj, reduction = 'umap.rna', group.by = 'ident') + ggtitle("RNA UMAP Clusters")
ggsave(file.path(result_path, "umap_rna.png"), plot = umap_rna, width = 8, height = 6,dpi = 300)
ggsave(file.path(result_path, "umap_rna.pdf"), plot = umap_rna, width = 8, height = 6, units = "in")

# # 可视化 ATAC UMAP 聚类
umap_atac <- DimPlot(lymph_obj, reduction = 'umap.atac', group.by = 'ident') + ggtitle("ATAC UMAP Clusters")
ggsave(file.path(result_path, "umap_atac.png"), plot = umap_atac, width = 8, height = 6,dpi = 300)
ggsave(file.path(result_path, "umap_atac.pdf"), plot = umap_atac, width = 8, height = 6, units = "in")

# 可视化加权 UMAP 聚类
umap_wnn <- DimPlot(lymph_obj, reduction = 'wnn.umap', group.by = 'ident') + ggtitle("WNN UMAP Clusters")
print(umap_wnn)
ggsave(file.path(result_path, "umap_wnn.png"), plot = umap_wnn, width = 8, height = 6,dpi = 300)
ggsave(file.path(result_path, "umap_wnn.pdf"), plot = umap_wnn, width = 8, height = 6, units = "in")

# 将图形组合在一起
combined_plot <- umap_rna + umap_atac + umap_wnn
print(combined_plot)
ggsave(file.path(result_path, "combined_plot.png"), plot = combined_plot, width = 8, height = 6,dpi = 300)
ggsave(file.path(result_path, "combined_plot.pdf"), plot = combined_plot, width = 8, height = 6, units = "in")


# obj_exps <- as.data.frame(as.matrix(GetAssayData(lymph_obj, assay = "RNA", layer = "data")))
# obj_atac<-as.data.frame(as.matrix(GetAssayData(lymph_obj, assay = "ATAC", layer = "data")))
# write.csv(obj_exps,file.path(result_path,"expssions.csv"))
# write.csv(obj_atac,file.path(result_path,"atacvalue.csv"))

# write.csv(lymph_obj@reductions[["wnn.umap"]]@cell.embeddings,file.path(result_path,"wnn_umap.csv"))
# write.csv(lymph_obj@reductions[["wnn.umap"]]@cell.embeddings,file.path(result_path,"wnn_umap.csv"))


## ----轨迹分析------------------------------------------------------------
sce <- as.SingleCellExperiment(lymph_obj, assay = "SCT") 
pca<- prcomp(t(GetAssayData(lymph_obj, assay = "SCT", slot = "data")), scale. = FALSE)
# 提取前两个主成分
rd1 <- pca$x[, 1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
reducedDims(sce) <- SimpleList(PCA = rd1)
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

gmm_classification <- colData(sce)$GMM
sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')

# 将 GMM 分类结果添加到 Seurat 对象的元数据
lymph_obj$GMM <- gmm_classification
head(lymph_obj$GMM)

# 按照 GMM 分类划分数据集
split_data <- SplitObject(lymph_obj, split.by = "GMM")
# 提取 PCA 降维结果
pca_coords <- reducedDims(sce)$PCA
write.csv(pca_coords,file.path(result_path,"pca.csv"))

# 提取所有轨迹的拟时数据
pseudotime_list <- slingPseudotime(sce)  # 获取所有轨迹的拟时
if (is.null(pseudotime_list)) {
  stop("Pseudotime information is not available in the SCE object.")
}

# 查看 Slingshot 轨迹
# 定义绘图函数
plot_slingshot <- function(output_file, format = "pdf", width = 8, height = 6, res = 300) {
  # 提取PCA降维结果和颜色
  pca_coords <- reducedDims(sce)$PCA
  cluster_colors <- brewer.pal(9, 'Set1')[sce$GMM]
  gmm_levels <- levels(as.factor(sce$GMM))
  
  # 根据文件格式选择图形设备
  if (format == "pdf") {
    pdf(output_file, width = width, height = height)
  } else if (format == "png") {
    png(output_file, width = width, height = height, units = "in", res = res)
  } else if (format=="tiff"){
    tiff(output_file, width = width, height = height, units = "in", res = res)
  }else{
    stop("Unsupported format. Please choose 'pdf' or 'png'.")
  }
  
  # 绘制PCA点图，颜色根据GMM聚类
  plot(pca_coords, 
       col = cluster_colors, 
       pch = 16, 
       asp = 1,
       main = "Slingshot Trajectory with GMM Clusters",
       xlab = "PCA 1", 
       ylab = "PCA 2")  # 添加x和y轴标签
  
  # 绘制Slingshot轨迹线
  slingshot_data <- SlingshotDataSet(sce)
  lines(slingshot_data, lwd = 2, type = 'lineages', col = 'blue')  # 画class lineage轨迹
  
  # 添加图例
  legend(
    x = "topright",
    legend = levels(as.factor(sce$GMM)), 
    fill = brewer.pal(9, 'Set1'), 
    title = "GMM Clusters", 
    cex = 0.8, 
    bty = "n")
  
  # 关闭图形设备
  dev.off()
}


# 调用绘图函数
# plot_slingshot(file.path(result_path, "slingshot_class.pdf"), format = "pdf")
plot_slingshot(file.path(result_path, "slingshot_class.png"), format = "png")
plot_slingshot(file.path(result_path, "slingshot_class.tif"), format = "tiff")




# 定义绘制轨迹拟时的函数
plot_trajectories_pseudotime <- function(output_file, format = "png", width = 8, height = 6, res = 500) {
  # 提取PCA降维结果
  pca_coords <- reducedDims(sce)$PCA
  
  # 获取轨迹数量和颜色
  num_trajectories <- ncol(pseudotime_list)
  colors <- brewer.pal(min(9, num_trajectories), "Set1")[1:num_trajectories]
  names(colors) <- colnames(pseudotime_list)
  
  # 设置图形设备
  if (format == "tiff") {
    tiff(output_file, width = width, height = height,res=res,units = "px", pointsize = 12)
  } else if (format == "png") {
    png(output_file, width = width, height = height, units = "in", res = res)
  } else {
    stop("Unsupported format. Please choose 'tiff' or 'png'.")
  }
  
  # 设置绘图参数
  par(mar = c(5, 5, 4, 6))  # 增加右边距以留出空间放置图例
  
  # 绘制PCA图，初始点为透明黑色
  plot(pca_coords, 
       col = rgb(0, 0, 0, 0.5), 
       pch = 16, 
       asp = 1,
       main = "Cell Pseudotime in PCA Space",
       xlab = "", 
       ylab = "")
  
  # 遍历每条轨迹，分轨迹绘制拟时
  for (i in seq_len(num_trajectories)) {
    pseudotime <- pseudotime_list[, i]
    
    # 处理拟时的颜色梯度
    pseudotime_colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
    pseudotime_scaled <- as.numeric(cut(pseudotime, breaks = 100))
    cell_colors <- pseudotime_colors[pseudotime_scaled]
    
    # 在PCA图上绘制点，颜色对应拟时
    points(pca_coords, 
           col = cell_colors, 
           pch = 16)
  }
  
  # 绘制颜色图例
  library(fields)  # 用于添加图例 (image.plot)
  
  # 调整图例位置和样式
  image.plot(legend.only = TRUE, 
             zlim = range(pseudotime_list, na.rm = TRUE), 
             col = colorRampPalette(brewer.pal(9, "YlOrRd"))(100), 
             horizontal = FALSE,    # 设置图例为水平显示
             legend.lab = "Pseudotime", 
             legend.line = 2,       # 调整标签与图例的距离
             legend.cex = 0.8,      # 调整图例文字大小
             legend.width = 1.2,    # 调整图例宽度
             legend.shrink = 0.8)   # 缩小图例比例
  
  # 关闭图形设备
  dev.off()
  cat("Visualization of all trajectories' pseudotime saved to", output_file, "\n")
}

# 设置保存路径
output_png_file <- file.path(result_path, "slingshot_all_trajectories_pseudotime.png")
output_pdf_file <- file.path(result_path, "slingshot_all_trajectories_pseudotime.tif")
plot_trajectories_pseudotime(output_png_file, format = "png")
plot_trajectories_pseudotime(output_pdf_file, format = "tiff")


# 合并两个轨迹的时间
pseudotime_data <- as.data.frame(pseudotime_list)
pseudotime_data$combined_time <- apply(pseudotime_data, 1, function(row) {
  if (!is.na(row[1]) && !is.na(row[2])) {
    return(mean(c(row[1], row[2])))
  } else if (!is.na(row[1])) {
    return(row[1])
  } else {
    return(row[2])
  }
})
pseudotime_data$diff_time <- ifelse(is.na(pseudotime_data[,1]) | is.na(pseudotime_data[,2]), 0, abs(pseudotime_data[,1] - pseudotime_data[,2]))
# 查看并保存结果
head(pseudotime_data)
fwrite(pseudotime_data, file = file.path(result_path,"time.csv"))


# Calculate peak-gene regulatory potential matrix
ATAC_gene_peak <- CalGenePeakScore(peak_count_matrix = lymph_obj@assays$ATAC@counts,organism = "GRCh38")
GAS_obj<- calculate_GAS(ATAC_gene_peak = ATAC_gene_peak, obj = lymph_obj)                                                                                                                       
GAS <- GAS_obj[[1]]

# 不保存可以注释掉，不然很耗时
# write.csv(GAS,file.path(result_path,"gas.csv"))
lymph_atac <- GAS_obj[[2]]

# 神经网络计算
HGT_result <- run_HGT(GAS = as.matrix(GAS),result_dir='D:/YYH/gas/Result/RNA_ATAC',data_type = 'scRNA_scATAC')


# cell簇的可视化
cell_hgt_matrix <- HGT_result[['cell_hgt_matrix']]
rownames(cell_hgt_matrix) <- colnames(GAS)
lymph_obj <- lymph_obj[, colnames(GAS)]
cell_hgt_matrix <- cell_hgt_matrix[colnames(GAS),]
colnames(cell_hgt_matrix)<- paste("HGT", 1:128, sep = "_")
write.csv(cell_hgt_matrix,file.path(result_path,"cell_hgt_matrix.csv"))


HGT_embedding <-CreateDimReducObject(embeddings = cell_hgt_matrix,
                                     key = "HGT_",
                                     assay = "RNA")
lymph_obj@reductions[['HGT']] <- HGT_embedding
lymph_obj <- FindVariableFeatures(lymph_obj, selection.method = "vst", nfeatures = 2000)
# lymph_obj <- ScaleData(lymph_obj, features = VariableFeatures(lymph_obj))
lymph_obj <-RunUMAP(
  lymph_obj,
  reduction = 'HGT',
  dims = 1:ncol(cell_hgt_matrix),
  reduction.name = "umap.rna",
  reduction.key = "rnaUMAP_"
)
lymph_obj <-FindNeighbors(lymph_obj,reduction = "HGT",dims = 1:ncol(cell_hgt_matrix))
DefaultAssay(lymph_obj) <- "RNA"
lymph_obj <- FindClusters(lymph_obj, resolution = 1)

png(file.path(result_path, "DimPlot_lymph_obj_rna.png"), width = 8, height = 6, units = "in", res = 300) 
DimPlot(lymph_obj, reduction = 'umap.rna',label=TRUE,repel=TRUE)
dev.off()



## 我们的下游
# 计算attention
df <- get_gene_module(lymph_obj, GAS, HGT_result[["attention"]])
# 检查重复条目
duplicated_rows <- df[duplicated(df[, c("node1", "node2")]), ]
print("重复条目：")
print(duplicated_rows)

gas_data <- acast(df, node1 ~ node2, value.var = "weight", fill = 0)
write.csv(t(gas_data),file.path(result_path,"gas_data_transpose.csv"))



# 确保分类信息中的细胞名顺序与 gas_data 的列名一致
gmm_classification <- data.frame(Cell=rownames(as.data.frame(colData(sce)$GMM)),Classification =colData(sce)$GMM,index=FALSE)
gmm_classification <- gmm_classification %>%
  filter(Cell %in% colnames(gas_data))  # 只保留 gas_data 中存在的细胞

# 按分类提取 gas_data 的子矩阵，并保留基因信息
gas_data_list <- gmm_classification %>%
  group_by(Classification) %>%
  summarise(Data = list(gas_data[, Cell, drop = FALSE]), .groups = "drop") %>%
  pull(Data)
names(gas_data_list) <- unique(gmm_classification$Classification)

# fwrite(as.data.frame(gas_data_list[["3"]], file = file.path(result_path,"3_rna.csv")))

for (i in 1:length(split_data)) {
  expression_matrix <- GetAssayData(object = split_data[[i]], assay = 'RNA', slot = "data")
  csv_name <- paste0(i, "_rna.csv")
  # 将GAS@Dimnames[[1]]作为行名，并将expression_matrix转换为数据框
  expression_matrix_df <- as.data.frame(expression_matrix, row.names = rownames(lymph_obj@assays[["RNA"]]))
  # 将数据框写入CSV文件
  fwrite(expression_matrix_df, file = file.path(result_path, csv_name),row.names = TRUE)
}

#### 划分分化前中后数据
# unsort_3k :(1,3) (5)(6)(2,4)

gas_data_merge <- list()# 创建 gas_data_merge 对象

# pbmc_unsort_3k
index<- list( # 定义索引规则
  c("1", "3"),  # 分化前：1 和 3
  c("5"),       # 分化中：5
  c("6"),       # 分化后1：6
  c("2", "4")   # 分化后2：2 和 4
)

# pbmc_sort_3k
# index<- list( # 定义索引规则
#   c("2", "5","7"),  # 分化前：2,5和7
#   c("8","3"),       # 分化中：8,3
#   c("1"),       # 分化后1：1
#   c("4", "6","9")   # 分化后2：4,6和9
# )

# pbmc_10k
#index<- list( # 定义索引规则
#   c("1","7"),  # 分化前：1和7
#   c("6"),       # 分化中：6
#   c("5"),       # 分化后1：5
#   c("2", "2","4")   # 分化后2：2,3和4
# )
for (i in seq_along(index)) {# 根据索引规则自动处理列拼接
  group <- index[[i]]  # 获取当前组的索引
  
  if (length(group) == 1) {
    # 如果只有一个索引，直接赋值
    gas_data_merge[[as.character(i)]] <- gas_data_list[[group]]
  } else {
    # 如果有多个索引，使用 cbind 进行列拼接
    gas_data_merge[[as.character(i)]] <- do.call(cbind, gas_data_list[group])
  }
}

fwrite(as.data.frame(gas_data_merge), file = file.path(result_path,"gas_data_merge.csv"))

# 获取在每个细胞中得分排名前500，在50%的细胞中都有表达的基因
# genes <- getTopGenesIntersection(gas_data_merge[["2"]],500,0.5)
genes <- getTopGenesIntersection(gas_data_merge[["2"]],500,0.5)
print(genes)


fwrite(as.data.frame(genes[["gene_intersection"]]), file = file.path(result_path,"genes.csv"))

# 保存lymph_obj
saveRDS(lymph_obj, file = file.path(result_path, paste0(file_name, ".rds")))
# 保存当前环境中的所有变量到文件
save(list = ls(), file = file.path(result_path, paste0(file_name, ".RData")))
#=========================================================

gene_expression <- data.frame()
gene<-genes[["gene_intersection"]]
# 遍历所有子集数据
for (i in 1:length(split_data)) {
  print(paste("Processing dataset:", i))
  
  # 确保基因名匹配，检查是否存在于当前对象中
  valid_genes <- intersect(gene, rownames(split_data[[i]]))
  if (length(valid_genes) == 0) {
    warning(paste("No valid genes found in dataset", i))
    next
  }
  # 提取表达矩阵
  expression_matrix <- GetAssayData(object = split_data[[i]], assay = 'RNA', slot = "data")
  # 计算基因的平均表达值（按行求均值）
  gene_avg <- rowMeans(expression_matrix[valid_genes, , drop = FALSE])
  gene_expression <- rbind(gene_expression, gene_avg)
}
colnames(gene_expression) <- gene



####### enrichR 包
dbs <- listEnrichrDbs()
head(dbs)
dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021","Reactome_Pathways_2024")

# 输入基因列表
gene <-genes[["gene_intersection"]]
gene.entrez <- bitr(gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db,
                    drop = FALSE)
# 运行富集分析
results <- enrichr(gene, dbs)
# 提取 KEGG 数据并选择前10个显著通路
data_path="KEGG_2021_Human"# "Reactome_Pathways_2024" 
enrichr_results <- results[[data_path]]
enrichr_top20 <- results[[data_path]][1:20, ]
plotEnrich(enrichr_results, showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value")
ggsave("KEGG_Enrichment.png", plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300,path=result_path)

go_result<-results[["GO_Biological_Process_2021"]]
go_top20 <- go_result[1:20, ]
plotEnrich(go_result, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
ggsave("GO_Enrichment.png", plot = last_plot(), width = 10, height = 8, units = "in", dpi = 300,path=result_path)



# 绘制气泡图
# 左侧子图：KEGG 富集分析结果
p1 <- plotEnrich(enrichr_top20, showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("KEGG Enrichment Analysis")

# 右侧子图：GO 富集分析结果
p2 <- plotEnrich(go_top20, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("GO Biological Process Enrichment Analysis")

# 合并左右子图
library(patchwork)
p_combined <- p1 + p2 + plot_layout(ncol = 2)

# 保存气泡图
ggsave("Enrichment_Analysis.png", plot = p_combined, width = 20, height = 8, units = "in", dpi = 300, path = result_path)


## 根据选择的基因在不同簇中的表达量差异，筛选出每个簇中显著高表达的基因

# 计算每个簇中基因的平均表达量
cluster_means <- sapply(names(split_data), function(cluster) {
  rowMeans(split_data[[cluster]])  # 计算当前簇的基因平均表达量
})

# 找到每个簇中显著高表达的基因（例如，表达量排名前5的基因）
top_genes_per_cluster <- apply(cluster_means, 2, function(means) {
  names(sort(means, decreasing = TRUE))[1:20]  # 取表达量排名前5的基因
})

# 计算每个簇中与 gene_list 重叠的基因占比
overlap_ratio <- apply(top_genes_per_cluster, 2, function(cluster_genes) {
  # 计算重叠基因的数量
  overlap_count <- sum(cluster_genes %in% gene)
  # 计算占比
  overlap_ratio <- overlap_count / length(cluster_genes)
  return(overlap_ratio)
})
overlap_ratio


## 绘制选定的基因整体表达情况
gene<-genes[["gene_intersection"]]
gene_expression_1 <- GetAssayData(object = lymph_obj, assay = 'RNA', slot = "data")
gene_expression<-as.matrix(gene_expression_1[gene,])
gene_expression<-t(gene_expression)
# 将数据框转换为长格式，以便ggplot2可以更容易地处理
gene_expression_long <- reshape2::melt(gene_expression, id.vars = "Gene")

ggplot(gene_expression, aes(x = pca_coords[,1], 
                            y = pca_coords[,2],
                            color = gene_expression[,1])) +
  geom_point(alpha = 0.3, size = 2) + # 绘制点，透明度为0.7，增加点的大小
  scale_color_gradient(low = "#1f77b4", high = "#d62728") + # 使用更专业的颜色渐变
  theme_minimal(base_size = 12, base_family = "Times New Roman") + # 使用简洁的主题，调整字体大小和字体类型
  theme(panel.grid.major = element_line(color = "grey90"), # 调整网格线颜色
        panel.grid.minor = element_blank(), # 移除次要网格线
        axis.ticks = element_blank(), # 移除坐标轴刻度线
        plot.title = element_text(hjust = 0.5)) + # 居中标题
  labs(title = "PCA of Gene Expression with BCL11B Expression Color", 
       x = "PCA1", 
       y = "PCA2") + # 添加标签
  theme(legend.position = "right") # 将图例放在右侧
ggsave("PCA_Gene_Expression.png", 
       plot = last_plot(), 
       width = 10, 
       height = 8, 
       units = "in",
       path=result_path,
       dpi = 300)

# 设置全局参数
output_dir <- file.path(result_path, "Mygene_new")
image_width <- 8
image_height <- 6
dpi_value <- 300

# 检查输出目录是否存在，如果不存在则创建
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 遍历基因表达数据的每一列
for (i in 1:dim(gene_expression)[2]) {
  # 生成UMAP图
  p <- ggplot(gene_expression, aes(x = pca_coords[,1], 
                                   y = pca_coords[,2],
                                   color = gene_expression[,i])) +
    geom_point(alpha = 0.7) + # 绘制点，透明度为0.7
    scale_color_gradient(low = "#1f77b4", high = "red", name = "gene_exp") + # 颜色渐变，从低到高，设置图例标题
    theme_minimal(base_size = 12, base_family = "Times New Roman") + # 使用简洁的主题
    theme(panel.background = element_rect(fill = "white"), # 设置白色背景
          plot.title = element_text(hjust = 0.5), # 居中标题
          legend.position = "right") + # 将图例放在右侧
    labs(title = paste("PCA of Gene Expression with ", colnames(gene_expression)[i], " Expression Color", sep=""),
         x = "PCA1", y = "PCA2") # 添加标签
  
  # 保存图像
  ggsave(
    filename = paste0(colnames(gene_expression)[i], "_gene_expression.png"), # 文件名
    plot = p, # 使用当前的图形对象
    width = image_width, # 宽度 (以英寸为单位)
    height = image_height, # 高度 (以英寸为单位)
    path = output_dir, # 输出目录
    dpi = dpi_value # 分辨率
  )
}



# #  查询相关基因
# library(biomaRt)
# 
# # 选择Ensembl数据库
# # ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
# # 获取与基因列表相关的基因
# related_genes <- getBM(attributes = c("hgnc_symbol", "description", "chromosome_name", "start_position", "end_position"),
#                        filters = 'hgnc_symbol',
#                        values = gene,
#                        mart = ensembl)
# 
# # 查看结果
# print(related_genes)



### 输入基因列表
gene <- genes[["gene_intersection"]]
# gene <- c('IAH1', '2510006D16RIK', 'WDR42A', 'ENTPD5', 'FBXO3', 'TRAP1', '4732466D17RIK', 'OXSM', 'CRADD', 'ADHFE1', '2410018G20RIK', 'PLEKHA7', '1810044D09RIK', 'ZFYVE20', 'TMEM80', 'LASS2', 'GSTZ1', 'CDAN1', 'PSMC6', 'TMEM77', 'ASB9', 'YME1L1', 'ASF1A', 'TFAM', '4921530L18RIK', 'TLN1', '4933407N01RIK', 'ALDH8A1', '39509', '2310026E23RIK', 'NLRX1', '2700046G09RIK', 'EI24', 'D4BWG0951E', 'STXBP2', 'IPP', 'CISD1', 'NFS1', 'GPR155', 'RAB11FIP2', 'LOC100044139', 'KDR', 'SEPHS2', 'EXOSC4', 'RILP', 'HOXA7', 'B3BP', 'AFAP1L1', 'TMEM86A', 'CNTD1', 'KLF12', '1700123L14RIK', 'APOOL', 'FZD5', 'TGDS', 'AI316807', 'VPS13B', 'MCAT', 'OVOL1', 'D130020L05RIK', 'PARP16', 'ELL3', 'COQ10A', 'WBSCR18', 'PLSCR2', 'VAMP8', 'TLCD1', 'TCN2', 'GORASP1', 'A930041I02RIK', 'SFXN5', 'FBXL3', 'SF1', 'FBXL6', 'PKIG', 'BC038156', '5730414N17RIK', 'CEP68', 'ATP6V1B2', 'UNC119B', 'HPN', 'OSGEPL1', 'RFESD', 'ZFP748', '3110048L19RIK', 'PTTG1IP', 'FARS2', 'TMEM186', 'ADH5', '4932432K03RIK', 'MRPL9', '1700023B02RIK', 'ZC3H12C', 'ARSK', 'SCP2', 'BC016495', 'CDK5RAP1', 'ARSG', 'LIPT1', 'PRKACA', 'C330002I19RIK', 'ZRSR1', 'LYPLA1', 'SLC33A1', 'NSUN3', '2410012H22RIK', '4632404H12RIK', 'IFT122', 'LOC100048387', 'LOC100047292', '2610036D13RIK', '1110032A03RIK', 'WDR34', 'ITFG1', 'LYRM2', 'CLDN10', 'PRPF18', 'ALDH6A1', 'NME7', 'LYRM5', 'UFC1', 'RAB1', 'ATPAF1', '9030420J04RIK', 'ABHD3', 'INTU', 'ADK', 'TFB1M', 'PCSK7', 'VLDLR', 'MRPL35', 'LRRC1', 'WDR24', 'KALRN', 'MYNN', 'UBIE', 'TMEM166', 'THTPA', '1700034H14RIK', 'PSMC3IP', '9230114K14RIK', 'MED14', 'LOC100047604', 'ATXN2', 'LOC623451', 'MYO6', '5430437P03RIK', 'TMED4', 'TASP1', 'AP4S1', 'SMYD4', '4921508D12RIK', 'PCMTD2', 'WDR20A', '2410166I05RIK', 'DHRS1', 'UBE2E1', 'GNMT', '2210016F16RIK', 'MPP7', 'MIPOL1', 'LOC674449', 'AGBL3', 'AI931714', 'COL4A4', 'TSR2', 'NOTUM', 'ANKRD42', 'NUPL2', 'RP23-195K8.6', 'ZFP148', 'ENSMUSG00000074286', 'LRRC19', 'MTMR14', 'GBE1', 'NUDT6', 'PHF7', 'ZDHHC5', '5730403B10RIK', 'SIP1', 'NUDT12', 'GYS2', 'ESM1', 'SIPA1L1', 'SPTLC1', 'ZFP11', 'PSMB1', 'NPY', 'METTL7A', 'RIOK2', 'METTL8', 'MDH1', 'AW209491', '2700038C09RIK', 'MGAT1', 'PAIP1', 'CD55', 'KLHDC4', '2610019F03RIK', 'TOR1A', 'FKBPL', 'A530050D06RIK', 'HSD3B2', 'TNFSF5IP1', 'DOLPP1', 'ATAD3A', 'LIFR', 'LRRC61', 'FAH', 'SLC7A6OS', 'SLC25A39', 'SLC9A6', 'UBOX5', 'LOC100046168', '7530414M10RIK', '9630013D21RIK', 'D630023B12RIK', '0610013E23RIK', '4732435N03RIK', 'SLC25A40', 'ZFAND1', 'ACAA1A', 'CREBL2', 'VWCE', 'YARS2', '4932438A13RIK', 'NOL7', 'CLCC1', '2810432D09RIK', '6530415H11RIK', 'PMPCB', 'NDUFV1', 'SIAE', 'PMS1', 'ARHGEF12', 'SLC30A6', 'TRPC2', 'CCDC16', 'LRRC40', 'POLRMT', 'TIMM44', 'LRRC44', 'KMO', 'DNAJC18', 'DNAJC19', 'DALRD3', 'DMXL1', 'CACNB4', 'NAT9', 'SMO', '4833426J09RIK', '4933403G14RIK', 'LOC100044400', '2310068J16RIK', 'ACBD4', 'FGFR4', 'RNF167', 'ZCCHC3', '1810049H13RIK', 'ZFP106', 'CNO', '2610528K11RIK', 'TPMT', '1200014M14RIK', 'BRI3', 'A930005H10RIK', 'ARHGAP18', 'C330018D20RIK') 
gene.entrez <- bitr(gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db,
                    drop = FALSE)

### GO 富集分析
GO <- enrichGO(gene.entrez$ENTREZID,
               OrgDb = 'org.Hs.eg.db',
               keyType = "ENTREZID",
               ont = "ALL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)



### KEGG 富集分析
KEGG <- enrichKEGG(gene.entrez$ENTREZID,
                   keyType = "kegg",
                   organism = 'hsa',
                   pvalueCutoff = 0.2)

### 设置全局主题
theme_set(
  theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold"),  # 调整标题字体大小
      axis.title = element_text(size = 16),                # 调整轴标题字体大小
      axis.text = element_text(size = 14),                 # 调整轴刻度标签字体大小
      legend.title = element_text(size = 16),              # 调整图例标题字体大小
      legend.text = element_text(size = 14)                # 调整图例标签字体大小
    )
)


### KEGG 富集分析绘图
# 柱状图
barplot(KEGG, showCategory = 10, title = 'KEGG Pathway')
ggsave(filename = "KEGG_Barplot.png", plot = last_plot(), width = 10, height = 8, dpi = 300, path = result_path)

# 点状图
dotplot(KEGG)
ggsave(filename = "KEGG_Dotplot.png", plot = last_plot(), width = 10, height = 8, dpi = 300, path = result_path)

# 气泡图
cnetplot(KEGG, circular = FALSE, colorEdge = TRUE)
ggsave(filename = "KEGG_Cnetplot.png", plot = last_plot(), width = 10, height = 8, dpi = 300, path = result_path)

# 热图
heatplot(KEGG, showCategory = 50)
ggsave(filename = "KEGG_Heatplot.png", plot = last_plot(), width = 10, height = 8, dpi = 300, path = result_path)


kegg_result <- as.data.frame(KEGG)
kegg_significant <- kegg_result[kegg_result$p.adjust < 0.05, ]

# 准备数据 - 选择前10个最显著的通路
top_pathways <- head(kegg_significant, 10)

# 创建用于ggalluvial的数据格式
alluvial_data <- data.frame()

for(i in 1:nrow(top_pathways)) {
  genes <- unlist(strsplit(top_pathways$geneID[i], "/"))
  pathway <- top_pathways$Description[i]
  
  # 转换ENTREZ ID为基因符号
  gene_symbols <- mapIds(org.Hs.eg.db, 
                         keys = genes,
                         column = "SYMBOL", 
                         keytype = "ENTREZID",
                         multiVals = "first")
  
  # 移除NA值
  gene_symbols <- gene_symbols[!is.na(gene_symbols)]
  
  temp_df <- data.frame(
    gene = gene_symbols,
    pathway = pathway,
    count = top_pathways$Count[i],
    pvalue = top_pathways$p.adjust[i],
    freq = 1
  )
  
  alluvial_data <- rbind(alluvial_data, temp_df)
}

# 为了更好的视觉效果，限制显示的基因数量
gene_freq <- table(alluvial_data$gene)
top_genes <- names(head(sort(gene_freq, decreasing = TRUE), 30))
filtered_data <- alluvial_data[alluvial_data$gene %in% top_genes, ]
library(ggalluvial)
# 方法1: 基础桑基图
basic_sankey <- ggplot(filtered_data, 
                       aes(axis1 = pathway, axis2 =gene )) +
  geom_alluvium(aes(fill = pathway), width = 1/12, alpha = 0.7) +
  geom_stratum(width = 1/12, fill = "lightblue", color = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 3, color = "black") +
  scale_x_discrete(limits = c("Gene", "Pathway"), expand = c(.05, .05)) +
  scale_fill_viridis_d(name = "KEGG pathway") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  labs(
    title = "The KEGG pathway enriches the Sankey plot",
  )

print(basic_sankey)

# 生成彩虹色向量
rainbow_colors <- rainbow(length(unique(filtered_data$pathway)))

# 方法1: 基础桑基图
basic_sankey <- ggplot(filtered_data, 
                       aes(axis1 = pathway, axis2 = gene)) +
  geom_alluvium(aes(fill = pathway), width = 1/12, alpha = 0.7) +
  geom_stratum(width = 1/12, fill = "lightblue", color = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 3, color = "black") +
  scale_x_discrete(limits = c("Gene", "Pathway"), expand = c(.05, .05)) +
  scale_fill_manual(name = "KEGG pathway", values = rainbow_colors) + # 使用彩虹色
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  labs(
    title = "The KEGG pathway enriches the Sankey plot",
    subtitle = "Pathways are colored with a rainbow gradient"
  )

print(basic_sankey)

# 保存为 PNG 文件
ggsave("sankey_plot.png", plot = basic_sankey, device = "png", width = 8, height = 6, dpi = 300, path = result_path)

# 保存为 TIFF 文件
ggsave("sankey_plot.tif", plot = basic_sankey, device = "tiff", width = 8, height = 6, dpi = 300, path = result_path)
# 方法2: 高级美化桑基图
# 创建自定义颜色调色板
n_pathways <- length(unique(filtered_data$pathway))
colors <- rainbow(n_pathways, alpha = 0.7)
advanced_sankey <- ggplot(filtered_data, aes(axis1 = pathway, axis2 = gene)) +
  geom_alluvium(aes(fill = pathway), width = 1/12,  alpha = 0.8,curve_type = "cubic") +
  geom_stratum(width = 1/12, aes(fill = pathway),color = "white", size = 0.5) +
  #geom_stratum(width = 1/12, fill = "lightblue", color = "white", size = 0.5) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, color = "black",fontface = "bold") +
  scale_x_discrete(limits = c("基因", "通路"), expand = c(.1, .1)) +
  scale_fill_manual(values = colors, name = "KEGG通路") +  # 统一设置颜色映射
  theme_void() +
  theme(
    legend.position = "none",  # 隐藏图例节省空间
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "KEGG通路富集分析桑基图",
    subtitle = paste("显示", length(unique(filtered_data$gene)), 
                     "个基因与", length(unique(filtered_data$pathway)), 
                     "个显著富集通路的关联"),
    caption = "数据来源: KEGG数据库富集分析"
  )

print(advanced_sankey)


# 
# ## 获取通路图
# library(pathview)
# for(i in KEGG@result[["ID"]]) {
#   # 构建输出文件名，包含ID和后缀
#   output_file <- file.path(result_path, paste0(i, "_pathview.pdf"))
# 
#   # 调用pathview函数，指定输出文件路径
#   pv.out <- pathview(cpd.data = i, pathway.id = i, species = "hsa", out.suffix = "pathview", out.file = output_file)
# }
