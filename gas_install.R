# 1. 安装和加载Bioconductor管理器
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. 安装和加载基因组和注释数据库
if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
  BiocManager::install("GenomeInfoDb")
}
if (!requireNamespace("IRanges", quietly = TRUE)) {
  BiocManager::install("IRanges")
}
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  BiocManager::install("rtracklayer")
}
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
  BiocManager::install("EnsDb.Hsapiens.v86")
}
if (!requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE)) {
  BiocManager::install("EnsDb.Mmusculus.v79")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

# 3. 安装和加载单细胞分析工具
if (!requireNamespace("scater", quietly = TRUE)) {
  BiocManager::install("scater")
}
if (!requireNamespace("bluster", quietly = TRUE)) {
  BiocManager::install("bluster")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
if (!requireNamespace("Signac", quietly = TRUE)) {
  install.packages("Signac")
}
if (!requireNamespace("rhdf5", quietly = TRUE)) {
  BiocManager::install("rhdf5")
}
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  BiocManager::install("glmGamPoi")
}
if (!requireNamespace("AUCell", quietly = TRUE)) {
  BiocManager::install("AUCell")
}

BiocManager::install(c("doMC", "doRNG", "doSNOW"))
BiocManager::install(c("mixtools", "SummarizedExperiment"))
BiocManager::install(c("DT", "plotly", "NMF", "d3heatmap", "dynamicTreeCut", "R2HTML", "Rtsne", "zoo"))
BiocManager::install("topGO")

# 4. 安装和加载可视化和绘图工具
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("patchwork", quietly = TRUE)) {
  devtools::install_github("thomasp85/patchwork")
}

# 5. 安装和加载数据处理和统计分析工具
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
if (!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel")
}

# 6. 安装和加载其他辅助工具
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages('RColorBrewer')
}
if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}
if (!requireNamespace("plyr", quietly = TRUE)) {
  install.packages('plyr')
}
if (!requireNamespace("dsb", quietly = TRUE)) {
  install.packages('dsb')
}
if (!requireNamespace("Hmisc", quietly = TRUE)) {
  install.packages("Hmisc")
}
if (!requireNamespace("CellChat", quietly = TRUE)) {
  devtools::install_github("sqjin/CellChat")
}
if (!requireNamespace("MAESTRO", quietly = TRUE)) {
  install_github("liulab-dfci/MAESTRO")
}
if (!requireNamespace("hdf5r", quietly = TRUE)) {
  install.packages("hdf5r")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("enrichR", quietly = TRUE)) {
  install.packages("enrichR")
}
if (!requireNamespace("limma", quietly = TRUE)) {
  install.packages("limma")
}
devtools::install_github("mojaveazure/seurat-disk")
devtools::install_github("10XGenomics/loupeR")
devtools::install_github('xzhoulab/SPARK')

