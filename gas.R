library(MAESTRO)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(scater)
library(Seurat)
library(Signac)
library(cluster)
library(bluster)
library(GenomeInfoDb)
library(igraph)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tinytex)
library(tidyverse)
library(rtracklayer)
library(reticulate)
library(parallel)
library(igraph)
library(data.table)
library(Hmisc)
library(rhdf5)
library(reshape2)
library(slingshot)
library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(RColorBrewer)
library(mclust, quietly = TRUE)
library(grDevices)
library(fields) # 用于添加图例 (image.plot)
library(tidyr)
library(enrichR)

library(ReactomePA)
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library("org.Hs.eg.db")
process_single_cell_data <- function(lymph_obj) {

  # 将 Seurat 对象转换为 SingleCellExperiment 对象
  sce <- as.SingleCellExperiment(lymph_obj, assay = "SCT")
  pca<- prcomp(t(GetAssayData(lymph_obj, assay = "SCT", slot = "data")), scale. = FALSE)
  rd1 <- pca$x[, 1:2]
  
  # 将 PCA 降维结果添加到 SCE 对象
  reducedDims(sce) <- SimpleList(PCA = rd1)
  
  # 使用 Mclust 进行聚类
  cl1 <- Mclust(rd1)$classification
  colData(sce)$GMM <- cl1
  
  # 将 GMM 分类结果添加到 Seurat 对象的元数据
  lymph_obj$GMM <- colData(sce)$GMM
  
  # 按照 GMM 分类划分 Seurat 对象
  split_data <- SplitObject(lymph_obj, split.by = "GMM")
  
  # 使用 Slingshot 进行轨迹推断
  sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
  
  # 返回结果
  return(list(sce = sce, split_data = split_data))
}

# Creat a seuart object
ReadData <- function(file_name = NULL, min_cell = 0.1, nmad = 3, gene.filter = TRUE, cell.filter = TRUE) {
  if (gene.filter == FALSE) {
    min_cell <- 0
  }
  h5_file <- paste0(file_name, '_filtered_feature_bc_matrix.h5')
  frag_file<-paste0(file_name, '_atac_fragments.tsv.gz')
  
  h5Path <- file.path(getwd(), 'Data', file_name, h5_file)
  fragPath<-file.path(getwd(), 'Data', file_name, frag_file)
  
  obj <- Read10Xdata(h5Path = h5Path,fragPath=fragPath, min_cell = min_cell)
  if (cell.filter == TRUE) {
    obj <- filterCell(obj, nmad = nmad)
  }
  return(obj)
}

# Read 10X data
##' Read matched scRNA + scATAC data from H5 file
# input:
#  h5Path: the path of h5 file
#  min_cell: the peak / gene will be removed if the value in the gene / peak with more than min_cell cell is equal to zero
# output:
#  a seurat object
Read10Xdata <-
  function(h5Path, fragPath=NULL,min_cell = 0.01,spices='hg38') {
    inputdata.10x <- Read10X_h5(h5Path)
    rna_counts <- inputdata.10x$`Gene Expression`
    atac_counts <- inputdata.10x$Peaks
    
    grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]
    # annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    # seqlevelsStyle(annotations) <- 'UCSC'
    # genome(annotations) <- spices
    
    chrom_assay <- CreateChromatinAssay(
      counts = atac_counts,
      sep = c(":", "-"),
      #genome = 'hg38',
      fragments = fragPath,
      min.cells = ncol(atac_counts) * min_cell,
      # annotation = annotations
    )

    tmp_obj <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
    exp_assay <- CreateAssayObject(counts = rna_counts, min.cells = ncol(rna_counts) * min_cell)
    tmp_obj[["RNA"]] <- exp_assay
    DefaultAssay(tmp_obj) <- "RNA"
    tmp_obj[["percent.mt"]] <-PercentageFeatureSet(tmp_obj, pattern = "^MT-")
    return(tmp_obj)
  }


# Filter abnormal cells
# input:
#  obj: a seurat object
#  nmad: a numeric scalar, specifying the minimum number of MADs away from median required for a value to be called an outlier
# output:
#  a seurat object

filterCell <- function(obj, nmad = 3) {
  atac <- isOutlier(obj$nCount_ATAC,
    nmads = nmad,
    log = F,
    type = "both"
  )
  rna <- isOutlier(obj$nCount_RNA,
    nmads = nmad,
    log = F,
    type = "both"
  )

  mito <- isOutlier(obj$percent.mt,
    nmads = nmad,
    log = F,
    type = "both"
  )
  obj <-
    AddMetaData(obj, atac, col.name = "atac")
  obj <-
    AddMetaData(obj, rna, col.name = "rna")
  obj <-
    AddMetaData(obj, mito, col.name = "mito")
  obj <- subset(
    x = obj,
    subset = atac == F &
      rna == F &
      mito == F
  )

  return(obj)
}

performAnalysis <- function(data) {
  # 设置默认Assay为RNA
  DefaultAssay(data) <- "RNA"
  
  # 增加全局选项的最大大小
  options(future.globals.maxSize = 600 * 1024^2)  # 600 MiB
  
  # 执行SCTransform、PCA和UMAP
  data <- SCTransform(data, verbose = FALSE) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  

  # 设置默认Assay为ATAC
  DefaultAssay(data) <- "ATAC"
  
  # 执行TF-IDF、查找顶级特征、SVD和UMAP
  data <- RunTFIDF(data)
  data <- FindTopFeatures(data, min.cutoff = 'q0')
  data <- RunSVD(data)
  data <- RunUMAP(data, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  
  # 确保 assay 设置为 SCT
  DefaultAssay(data) <- "SCT"
  # 使用 Seurat 的 FindVariableFeatures 提取高度可变基因
  data <- FindVariableFeatures(data, assay = "SCT", selection.method = "vst", nfeatures = 2000)
  # 检查高度可变基因数量
  variable_genes <- VariableFeatures(data)
  if (length(variable_genes) == 0) stop("没有找到高度可变基因，请检查数据或调整参数！")
  # 2. PCA 降维 --------------------------------------------------------
  # 提取高度可变基因的表达矩阵
  hvg_data <- GetAssayData(data, assay = "SCT", slot = "data")[variable_genes, ]
  
  # 计算多模态邻居并进行加权 UMAP 和聚类
  data <- FindMultiModalNeighbors(data, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  data <- RunUMAP(data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  data <- FindClusters(data, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  return(data)
}


##' Calculate gene active score matrix
# input:
#  peak_count_matrix: a peak_count matrix from scATAC-seq with peak * cell which return from filterCell function
#  organism: species type GRCh38 / GRCm38
# output:
#  a gene * peak matrix, the elements represent the regulatory potential for peak to gene

CalGenePeakScore <-
  function(peak_count_matrix, organism = "GRCh38") {
    pbmc_peak <- peak_count_matrix
    n <- nrow(pbmc_peak)
    dia <- diag(n)
    rownames(dia) <- rownames(pbmc_peak)
    colnames(dia) <- 1:ncol(dia)
    gene_peak <-
      ATACCalculateGenescore(dia,
        organism = organism,
        decaydistance = 10000,
        model = "Enhanced"
      )
    colnames(gene_peak) <- rownames(peak_count_matrix)
    return(gene_peak)
  }


ATACCalculateGenescore <- function(inputMat, project = "MAESTRO.scATAC", organism = "GRCh38",
                                   decaydistance = 10000, model = "Enhanced") {
  source_python(paste(system.file(package = "MAESTRO"), "ATACCalculateGenescore.py",
    sep = "/"
  ))
  peaks_list <- rownames(inputMat)
  peaks_list <- gsub(pattern = "\\W", replacement = "_", x = peaks_list)
  if (class(inputMat[1]) != "dgCMatrix") {
    inputMat <- as(as.matrix(inputMat), "dgCMatrix")
  }
  if (model == "Simple") {
    if (organism == "GRCh38") {
      data(GRCh38.refgenes.genescore.simple)
      refgenes.genescore <- GRCh38.refgenes.genescore.simple
    } else {
      data(GRCm38.refgenes.genescore.simple)
      # refgenes.genescore = GRCm38.refgenes.genescore.simple
    }
    genes_list <- refgenes.genescore[, 4]
    refgenes.genescore <- refgenes.genescore[, -4]
    rp_result <- calculate_RP_score(
      cell_peaks = inputMat,
      peaks_list = peaks_list, gene_bed_df = refgenes.genescore,
      genes_list = genes_list, decay = decaydistance,
      model = model
    )
  }
  if (model == "Enhanced") {
    if (organism == "GRCh38") {
      data(GRCh38.refgenes.genescore.adjusted)
      refgenes.genescore <- GRCh38.refgenes.genescore.adjusted
    } else {
      data(GRCm38.refgenes.genescore.adjusted)
      refgenes.genescore <- GRCm38.refgenes.genescore.adjusted
    }
    rp_result <- calculate_RP_score(
      cell_peaks = inputMat,
      peaks_list = peaks_list, gene_bed_df = refgenes.genescore,
      genes_list = NULL, decay = decaydistance, model = model
    )
  }
  rp_matrix <- rp_result[[1]]
  rownames(rp_matrix) <- rp_result[[2]]
  colnames(rp_matrix) <- colnames(inputMat)
  return(rp_matrix)
}



##' Calculate gene active score matrix
# input:
#  ATAC_gene_peak: a matrix with gene * peak which return from CalGenePeakScore fucntion
#  obj: a seurat object after data preprocessing which return from filterCell function
#  method: the method to integrate scRNA-seq and scATAC-seq velo (velocity) / WNN (weighted nearest neighbor)
#  return.weight: if return.weight = T, return modality integrated weight, else return GAS
# output:
#  GAS matrix with gene * peak, the elements represent the gene activity score in each cell
#  obj: a seurat object with obj[['ATAC_active']] a gene * peak matrix, the elements represent the regulatory potential for peak to gene
calculate_GAS <-
  function(ATAC_gene_peak,obj,method = "wnn") {
    WA <- ATAC_gene_peak
    obj <- lymph_obj
    method <- "wnn"
    peak_count <- obj@assays$ATAC@counts
    gene_count <- obj@assays$RNA@counts
    peak_count[peak_count > 0] <- 1
    WA <- ATAC_gene_peak %*% peak_count
    colnames(WA) <- colnames(peak_count)
    rownames(WA) <- rownames(ATAC_gene_peak)
    write(as.matrix(WA),'WA.csv')
    print(colnames(WA)[1:3])
    print(colnames(obj)[1:3])
    
    # 去除所有元素为 0 的行。
    WA <- WA[which(rowSums(as.matrix(WA)) > 0), ]
    gene_count <- gene_count[which(rowSums(as.matrix(gene_count)) > 0), ]
    commongene <- intersect(x = rownames(WA), y = rownames(gene_count))
    WA <- as.matrix(WA)
    WA <- WA[commongene, ]
    atac_active <- CreateAssayObject(
      counts = WA,
      min.cells = 0
    )
    obj[["ATAC_active"]] <- atac_active
    gene_count <- gene_count[commongene, ]
    gene_rowsum <- rowSums(gene_count)
    peak_rowsum <- rowSums(WA)
    norm_gene_count <- gene_count / rowSums(gene_count)
    norm_WBinary <- WA / rowSums(WA)
    # norm_gene_count<-NormalizeData(CreateSeuratObject(counts = gene_count))$ RNA@data
    gene_count <- norm_gene_count
    # norm_WBinary<-NormalizeData(CreateSeuratObject(counts = WA))$RNA@data
    peak_count <- norm_WBinary
    # print(str(peak_count))
    if (method == "wnn") {
      obj <- obj[, colnames(gene_count)]
      DefaultAssay(obj) <- "RNA"
      obj <-
        FindVariableFeatures(obj,
          selection.method = "vst",
          nfeatures = 2000
        )
      obj <- ScaleData(obj, features = VariableFeatures(obj))
      obj <- RunPCA(obj)

      DefaultAssay(obj) <- "ATAC"
      obj <- RunTFIDF(obj)
      obj <- FindTopFeatures(obj, min.cutoff = "q0")
      obj <- RunSVD(obj)
      obj <-
        RunUMAP(
          obj,
          reduction = "lsi",
          dims = 2:50,
          reduction.name = "umap.atac",
          reduction.key = "atacUMAP_"
        )
      obj <-
        FindMultiModalNeighbors(obj,
          reduction.list = list("pca", "lsi"),
          dims.list = list(1:50, 2:50)
        )
      GAS <-
        gene_count * gene_rowsum * obj$RNA.weight + peak_count * obj$ATAC.weight *
          peak_rowsum
    }
      m <- list()
      m[[1]] <- GAS
      m[[2]] <- obj
      return(m)
  }

# required package -- reticulate
# input:
#  GAS: the spliced and normalized matrix obtained from CLR function
#  result_dir: The address for storing the models and optimization results(Type:str)
#  epoch:(Type:int)
#  lr: learning rate(Type:float)
#  n_hid: Number of hidden dimension(Type:int)
#  n_heads: Number of attention head(Type:int)
#  cuda: 0 use GPU0 else cpu(Type:int)
#  data_type: 'CITE', 'RNA_ATAC', or 'multiple RNA'
#  envPath: The address for environment to use if use.env is TRUE(Type:str)
# output:
#  HGT_result: a list containing requried results of HGT model as follows:
#  parameters: given parameters from user --epoch, lr, n_hid, n_heads, cuda
#  cell_hgt_matrix: cell embedding matrix
#  feature_hgt_matrix : gene embedding matrix and protein embedding matrix when data_type is 'CITE';
#  attention: attention meassage for features and cells
#  data_type: 'CITE', 'RNA_ATAC', or 'multiple RNA'
#  result_dir: The address for storing the models and optimization results
#  GAS: the spliced and normalized matrix obtained from CLR function

run_HGT <- function(GAS, result_dir, data_type, envPath = NULL, lr = NULL, epoch = NULL, n_hid = NULL, n_heads = NULL, cuda = 0) {
  if (data_type == "scRNA_scATAC") {
    if (is.null(lr)) {
      lr <- 0.1
    }
    if (is.null(epoch)) {
      epoch <- 100
    }
    if (is.null(n_hid)) {
      n_hid <- 128
    }
    if (is.null(n_heads)) {
      n_heads <- 16
    }
  }
  print(epoch)
  cat(lr, epoch, n_hid, n_heads, cuda)
  if (!is.null(envPath)) {
    use_condaenv(envPath)
  }
  list_in <- assign("list_in", list(lr = lr, epoch = epoch, n_hid = n_hid, n_heads = n_heads, result_dir = result_dir, cuda = cuda, data_type = data_type, cell_gene = GAS, gene_name = rownames(GAS), cell_name = colnames(GAS)), envir = .GlobalEnv)
  source_python("./arg.py")
  cell_hgt_matrix <- py$cell_matrix
  gene_hgt_matrix <- py$gene_matrix
  attention <- py$df2
  rownames(cell_hgt_matrix) <- list_in$cell_name
  rownames(gene_hgt_matrix) <- list_in$gene_name
  HGT_result <- list()
  HGT_result[["parameters"]] <- data.frame(lr, epoch, n_hid, n_heads, cuda)
  HGT_result[["GAS"]] <- GAS
  HGT_result[["cell_hgt_matrix"]] <- cell_hgt_matrix
  HGT_result[["feature_hgt_matrix"]] <- gene_hgt_matrix
  HGT_result[["attention"]] <- attention
  HGT_result[["result_dir"]] <- result_dir
  HGT_result[["data_type"]] <- data_type
  return(HGT_result)
}


umap_HGT<-function(HGT_result,save_path){
  cell_hgt_matrix <- HGT_result[['cell_hgt_matrix']]
  rownames(cell_hgt_matrix) <- colnames(GAS)
  lymph_obj <- lymph_obj[, colnames(GAS)]
  cell_hgt_matrix <- cell_hgt_matrix[colnames(GAS),]
  colnames(cell_hgt_matrix)<- paste("HGT", 1:128, sep = "_")
  write.csv(cell_hgt_matrix,file.path(result_path,"cell_hgt_matrix.csv"))
  
  HGT_embedding <-
    CreateDimReducObject(embeddings = cell_hgt_matrix,
                         key = "HGT_",
                         assay = "RNA")
  lymph_obj@reductions[['HGT']] <- HGT_embedding
  lymph_obj <- FindVariableFeatures(lymph_obj, selection.method = "vst", nfeatures = 2000)
  lymph_obj <- ScaleData(lymph_obj, features = VariableFeatures(lymph_obj))
  lymph_obj <-RunUMAP(
    lymph_obj,
    reduction = 'HGT',
    dims = 1:ncol(cell_hgt_matrix),
    reduction.name = "umap.rna",
    reduction.key = "rnaUMAP_"
  )
  lymph_obj <-FindNeighbors(lymph_obj,reduction = "HGT",dims = 1:ncol(cell_hgt_matrix))
  lymph_obj <- FindClusters(lymph_obj, resolution = 1)
  DefaultAssay(lymph_obj) <- "RNA"
  
  pdf(path, width = 8, height = 6)
  DimPlot(lymph_obj, reduction = 'umap.rna',label=TRUE,repel=TRUE)
  # dev.off()
}

get_gene_module <- function(obj, GAS, att, cutoff = 1.6) {
  # 在R中，Idents(obj) 函数通常与单细胞基因组学数据分析库 Seurat 相关联。
  # Seurat 是一个流行的R包，用于单细胞RNA测序数据的分析和探索。
  # Idents 函数用于获取或设置 Seurat 对象中细胞的活性身份，这些身份通常用于区分不同的细胞群或条件。
  # obj=lymph_obj
  # att=HGT_result[["attention"]]
  graph.out <- Idents(obj)
  nhead <- ncol(att)
  gene_name <- rownames(GAS)[att$gene + 1]
  cell_name <- colnames(GAS)[att$cell + 1]

  att$ct <- graph.out[cell_name]
  att$gene_name <- gene_name
  att$cell_name <- cell_name
  mod <- function(x) {
    return(sqrt(sum(c(x^2))))
  }
  nor <- function(x) {
    return((x - min(x)) / (max(x) - min(x)))
  }

  att[, 3:nhead] <- nor(att[, 3:nhead])
  attention <-
    aggregate(
      x = as.list(att[, 3:nhead]),
      by = list(att$ct, att$gene_name),
      mean
    )
  weight <- apply(att[, 3:nhead], 1, mod)
  df <-
    data.frame(
      "node1" = att$gene_name,
      "node2" = att$cell_name,
      "weight" = weight,
      "ct" = att$ct
    )
  return(df)
}

get_cell_time <- function(lymph_obj) {
  # sce <- SingleCellExperiment(assays = list(counts = GAS))
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(as.matrix(GetAssayData(lymph_obj, assay = "RNA", layer = "data")))))
  sce <- logNormCounts(sce)
  sce <- runPCA(sce)
  clus <- kmeans(reducedDim(sce, "PCA"), centers = 5)$cluster
  # colLabels(sce) <- factor(clus)
  # sce <- slingshot(sce, clusterLabels = "label", reducedDim = "PCA")
  sce <- slingshot(sce, clusterLabels = lymph_obj@active.ident, reducedDim = "PCA")
  pseudotime <- slingPseudotime(SlingshotDataSet(sce))
  pseudotime_data <- as.data.frame(pseudotime)

  if (ncol(pseudotime_data) == 1) {
    # 只有一列，直接使用其值
    print("只有一条轨迹")
    pseudotime_processed <- pseudotime_data
  }
  # } else {
  #   # 多列时处理
  #   pseudotime_processed <- apply(pseudotime_data, 1, function(row) {
  #     if (all(!is.na(row))) {
  #       # 如果每一行都有值，取平均值
  #       return(mean(row, na.rm = TRUE))
  #     } else {
  #       # 用0填充缺失值，其他值保持不变
  #       row[is.na(row)] <- 0
  #       return(mean(row))
  #     }
  #   })
  # }
  #
  # cell_names <- Cells(gas_object)
  # pseudotime_processed <- pseudotime_processed[cell_names, ]

  return(pseudotime_processed) 
}

get_cell_class <- function(pseudotime_data, df) {
  # 稀疏矩阵转为密集矩阵
  gas_data <- acast(df, node1 ~ node2, value.var = "weight", fill = 0)
  gas_object <- CreateSeuratObject(count = as(gas_data, "dgCMatrix"))
  # gas_object <- AddMetaData(gas_object, metadata = pseudotime_data$Lineage1, col.name = "pseudotime")
  # 将 pseudotime_data 的每一列分别添加到 gas_object 中
  for (col_name in colnames(pseudotime_data)) {
    # print(paste0("pseudotime_", col_name))
    gas_object <- AddMetaData(gas_object, metadata = pseudotime_data[[col_name]], col.name = paste0("pseudotime_", col_name))
  }
  # 将伪时间分为三类
  pseudotime_data$class <- cut(pseudotime_data$Lineage1, breaks = 3, labels = c("early", "mid", "late"))

  pseudotime_data <- pseudotime_data[Cells(gas_object), ]

  # 将分类信息添加到 Seurat 对象
  gas_object <- AddMetaData(gas_object, metadata = pseudotime_data$class, col.name = "pseudotime_class")

  # 提取伪时间分类为 'mid' 的细胞名称
  # early_cells <- rownames(gas_object@meta.data)[gas_object@meta.data$pseudotime_class == "early"]
  mid_cells <- rownames(gas_object@meta.data)[gas_object@meta.data$pseudotime_class == "mid"]
  # late_cells <- rownames(gas_object@meta.data)[gas_object@meta.data$pseudotime_class == "late"]

  mid_seurat_object <- subset(gas_object, cells = mid_cells)
  
  return(mid_seurat_object)
}

getTopGenesIntersection <- function(seurat_object, num = 500, n = 0.8) {
  # seurat_object: Seurat 对象或矩阵
  # num: 每个细胞中取表达量前 num 的基因
  # n: 交集的阈值比例（例如 0.8 表示基因需要在至少 80% 的细胞中出现）
  
  # 获取数据矩阵
  if (inherits(seurat_object, "Seurat")) { 
    data_matrix <- GetAssayData(seurat_object, slot = "count")
  }else{
    data_matrix<-as.matrix(seurat_object)
  }
  
  # 获取每个细胞中得分前num的基因
  top_genes_list <- apply(data_matrix, 2, function(x) {
    # 获取前num个基因的名称
    top_genes <- rownames(data_matrix)[order(x, decreasing = TRUE)[1:num]]
    return(top_genes)
  })

  # 将列表转换为向量，统计每个基因出现的频率
  all_top_genes <- unlist(top_genes_list)
  gene_counts <- table(all_top_genes)
  
  # 计算交集基因的阈值
  threshold <- n * ncol(top_genes_list)
  gene_counts[gene_counts >= threshold]
  
  # 筛选出频率大于等于阈值的基因
  selected_genes <- names(gene_counts[gene_counts >= threshold])

  return(list(top_genes_list = top_genes_list, gene_intersection = selected_genes))
}
