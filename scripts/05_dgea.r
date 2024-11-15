# Load required libraries
library(ggtree)
library(Seurat)
library(data.table)
library(dplyr)
library(tidytree)
library(scMayoMap)

# Define paths
cell_path <- '/home/celia/Documents/VisualStudioCode/WERBA/Werba/Werba.Cells.Annotation.V1.tsv'
expression_mtx_path <- '/home/celia/Documents/VisualStudioCode/WERBA/Werba/Werba.rawcount.V1.mtx'
features_path <- '/home/celia/Documents/VisualStudioCode/WERBA/Werba/Werba.geneAnnotation.V1.tsv'
subset_base_path <- '/home/celia/Documents/VisualStudioCode/WERBA/cellstatesappli'
output_path <- "/home/celia/Documents/annotation/findmarkers/werba"

# Define subsets to process
sample_file <- "/home/celia/Documents/VisualStudioCode/samples/werba_samples.txt"
subsets <- readLines(sample_file)

# Load cell annotation data
cell <- qtl2::fread_csv(cell_path, sep = "\t")

# Load expression matrix
expressionMatrix <- ReadMtx(
  mtx = expression_mtx_path, 
  cells = cell_path, 
  features = features_path,
  feature.column = 2, skip.cell = 1, cell.sep = '\t', feature.sep = '\t', skip.feature = 1
)

# Create Seurat object and normalize data
obj <- CreateSeuratObject(counts = expressionMatrix)
obj$cell <- colnames(obj)
obj <- NormalizeData(object = obj)

# Add metadata and process each subset
for (subset in subsets) {
  cat("Starting processing subset:", subset, "\n")
  
  # Define paths for the current subset
  RNA_matrix_path <- file.path(subset_base_path, paste0('RNAmatrix_', subset, '.tsv'))
  clusters_path <- file.path(subset_base_path, paste0('RNAmatrix_', subset, '/optimized_clusters.txt'))
  
  # Load cluster data and RNA matrix header
  clusters <- fread(clusters_path)
  data <- read.table(file = RNA_matrix_path, header = TRUE, nrows = 1, sep = '\t')
  
  # Update column names
  col_names <- gsub("\\.", "-", colnames(data))
  colnames(data) <- col_names
  
  # Subset cell annotation data
  cells_subset <- subset(cell, rownames(cell) %in% colnames(data))
  cells_subset$state <- clusters$V1
  cells_subset$sample <- subset
  
  # Add metadata to the Seurat object
  obj <- AddMetaData(object = obj, metadata = cells_subset)
  cat("Finished processing subset:", subset, "\n\n")
}

# Process nodes for each subset
for (subset in subsets) {
  cat("Starting processing subset:", subset, "\n")
  
  # Load subset Seurat object and tree
  obj_subset <- subset(obj, sample == subset)
  sub_path <- file.path(subset_base_path, paste0('RNAmatrix_', subset))
  tree_path <- file.path(sub_path, 'newick_tree.txt')
  t <- read.tree(file = tree_path)
  nodes <- setdiff(c(t$tip.label, t$node.label), "I0")
  tree_tibble <- as_tibble(t)
  
  # Process each node
  for (node in nodes) {
    cat("Processing node:", node, "\n")
    
    tryCatch({
      # Determine clusters and update metadata
      if (startsWith(node, "I")) {
        children <- offspring(tree_tibble, node)
        clusters <- children$label[grep("^C", children$label)]
        clust_num <- as.numeric(gsub("[^0-9]", "", clusters))
        
        obj_subset@meta.data$cluster <- ifelse(obj_subset@meta.data$state %in% clust_num, clust_num[1], clust_num[1] + 1)
        Idents(object = obj_subset) <- "cluster"
        
        markers_subset <- FindMarkers(obj_subset, ident.1 = clust_num[1], min.cells.group = 10)
      } else {
        clust_num <- as.numeric(gsub("[^0-9]", "", node))
        obj_subset@meta.data$cluster <- ifelse(obj_subset@meta.data$state %in% clust_num, clust_num, clust_num + 1)
        Idents(object = obj_subset) <- "cluster"
        
        markers_subset <- FindMarkers(obj_subset, ident.1 = clust_num, min.cells.group = 10)
      }
      
      # Save markers to file if found
      if (!is.null(markers_subset) && nrow(markers_subset) > 0) {
        subset_dir <- file.path(output_path, subset)
        if (!dir.exists(subset_dir)) dir.create(subset_dir, recursive = TRUE)
        csv_file <- file.path(subset_dir, paste0(node, ".csv"))
        write.csv(markers_subset, file = csv_file, row.names = TRUE)
      }
      
    }, error = function(e) {
      cat("Error processing node:", node, "-", e$message, "\n")
    })
  }
}
