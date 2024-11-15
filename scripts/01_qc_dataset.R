# Load required libraries
library(scuttle)
library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(glmGamPoi)
library(ggplot2)
library(SeuratData)
library(SeuratDisk)
library(Matrix)

# Load expression matrix
expressionMatrix <- ReadMtx(
  mtx = 'Exp_data_UMIcounts.mtx', 
  cells = 'Cells.csv', 
  features = 'Genes.txt', 
  feature.column = 1, 
  skip.cell = 1, 
  cell.sep = ','
)

# Create Seurat object
obj <- CreateSeuratObject(
  counts = expressionMatrix, 
  min.cells = 3, 
  project = "sc"
)

# Add mitochondrial percentage metadata
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Set object identity
Idents(obj) <- obj@meta.data$orig.ident

# Subset data based on feature counts and mitochondrial content
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 25)

# Extract RNA matrix and filter based on expression
RNAmatrix <- GetAssayData(obj, assay = "RNA", slot = "counts")
subRNAmatrix <- RNAmatrix[rowMeans(RNAmatrix > 0) > 0.01, ]

# Save filtered RNA matrix to file
write.table(subRNAmatrix, 'RNAmatrix.tsv', sep = '\t', quote = FALSE, col.names = NA)
