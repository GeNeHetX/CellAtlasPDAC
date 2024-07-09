# Introduction
This document outlines the process of creating a cell atlas through a series of scripts. Each script performs a specific step in the workflow, and the inputs and outputs of each step are clearly defined.

Each dataset must follow the first 6 scripts sequentially. The next 6 scripts should only be executed after executing the first 6 scripts for each dataset.

When adding a new dataset, we have to execute the first 6 scripts for this new dataset then execute the next 6 scripts. Also, the output directory of the 7th script should be the same as the one used when creating the atlas.

# Step 1: Data Loading and Preprocessing

## Summary

This script loads a single-cell RNA-seq dataset using Seurat, performs basic quality control (QC) steps to filter out low-quality cells based on feature counts and mitochondrial content, and saves the filtered RNA matrix for downstream analysis.

## Inputs

- **Expression Matrix**: Exp_data_UMIcounts.mtx
- **Cells metadata**: Cells.csv
- **Genes**: Genes.txt
- **Metadata**: Meta-data.csv

## Outputs

- **Filtered RNA Matrix**: Saved as `RNAmatrix.tsv` in the current working directory. This file contains the expression data of filtered cells and genes after preprocessing.

## Script
```
Rscript 01_qc_dataset.R
```

# Step 2: Sample-based Filtering and Directory Creation

## Summary

This script filters the RNA matrix based on sample IDs extracted from column names, and saves the filtered data for each sample into separate files (`RNAmatrix_{sample_id}.tsv`). Additionally, it creates directories named after each sample ID to store eventual *cellstates* outputs.

## Inputs

- **RNA Matrix**: Path to a RNA matrix of a dataset, e.g., `/cellstatesappli/RNAmatrix.tsv`.
- **User Choice**: Users are prompted to choose how cell barcodes are separated per samples (e.g., `_`, `-`, etc.).

## Outputs

- **Filtered RNA Matrices**: Saved as `RNAmatrix_{sample_id}.tsv` for each sample ID in the directory `/cellstatesappli/`.
- **Directories**: Directories named `RNAmatrix_{sample_id}` are created in `/cellstatesappli/` to store eventual *cellstates* outputs.

## Script: 02_sample_matrix.py

```
python3 02_sample_matrix.py
```

# Step 3: Running Cell States Analysis

## Summary

This script runs a Python script (`run_cellstates.py`) for each sample listed in `ITEMS`. It specifies input file paths (`INPUT_FILE`) and output directories (`OUTPUT_DIR`) for each sample.

## Inputs

- **Base Directory**: Base path that contains where all outputs will be generated, e.g.,`/cellstatesappli/`
- **Script Path**: Path to the script that runs cellstates `/run_cellstates.py`
- **List of Samples**: `ITEMS`, containing sample identifiers.

## Outputs

- **Cell States Analysis Results**: Generated in each output directory (`OUTPUT_DIR`) for every sample processed. The directories are generated in the script `02_sample_matrix.py`.

## Script:

```
bash 03_cellstates_samples.sh
```

# Step 4: Generate Newick Tree Strings for Clusters

## Summary

This Python script processes RNA matrices, cluster assignments, Dirichlet pseudocounts, and cluster hierarchies for multiple subsets within a specified dataset. It generates Newick tree strings for each subset and writes them to output files.

## Inputs

- **Base Path**: Base path containing all cellstates outputs
- **Dataset Name**: Name of the folder containing cellstates outputs of a specific dataset
- **Subsets**: A list of subset names to process, e.g., `['P01', 'P02', 'P04', 'P05', 'P06', 'P07', 'P08', 'P09', 'P10', 'P11', 'P12', 'P13']`

## Outputs

- **Newick Tree Strings**: Generated and saved to `newick_tree.txt` files within respective subset directories.

## Script

```
python 04_generate_newick.py
```
## Dependencies

- **CellStates Environment**: Ensure the `cellstates-env` Conda environment is activated before running the script.

```bash
conda activate cellstates-env
```

# Step 5: Identify Markers with DGEA

## Summary

This R script processes RNA matrices, cell annotations, and cluster data for multiple subsets within a specified dataset. It generates a Seurat object, normalizes the data, adds metadata, processes nodes for each subset, and identifies marker genes. The results are saved to CSV files.

## Inputs

- **Cell Annotation Path**: Path to file containing the cell annotation data or cell barcodes of a dataset, e.g., `/Zhang.Cells.Annotation.V1.tsv`
- **Expression Matrix Path**: Path to file containing the expression matrix of a dataset in 'mtx' format, e.g.,`/Zhang.rawcount.V1.mtx`
- **Features Path**: Path to file containing gene annotation data or gene list of a dataset, e.g., `/Zhang.geneAnnotation.V1.tsv`
- **Subset Matrix Base Path**: Path to folder containing matrices of all subsets of a dataset
- **Output Path**: Path to output folder
- **Subsets to Process**: A list of subsets to process, e.g., `['Case1.YF', 'Case1.ZY', 'Case2.YF', 'Case2.ZC', 'Case2.ZY', 'Case3.YF', 'Case3.ZY', 'Case4.ZY']`

## Outputs

- **Markers CSV Files**: Marker genes for each node in each subset, saved in the output directory.

## Dependencies

- **Newick Tree File**: Ensure that all Newick Tree files are generated by *Script 4* for all samples before running the script.

## Script

```
Rscript 05_dgea.r
```

# Step 6: Create Node-Cell Matrix

## Summary

This script generates a node-cell matrix using node IDs from Newick tree files and RNA matrix data for multiple subsets. It utilizes multiprocessing to efficiently process nodes and outputs a sparse matrix in .npz format.

## Inputs

- **Base Path**: Base path containing all cellstates outputs
- **Dataset Name**: Name of the folder containing cellstates outputs of a specific dataset
- **Subsets to Process**: A list of subsets to process, e.g., `['Case1.YF', 'Case1.ZY', 'Case2.YF', 'Case2.ZC', 'Case2.ZY', 'Case3.YF', 'Case3.ZY', 'Case4.ZY']`
- **Output File Path**: Path to output file in .npz format, e.g., `/sparse_mtx/cell_mtx_zhang.npz`
- **Number of Processes**: Number of processes for multiprocessing

## Outputs

- **Sparse Node-Cell Matrix**: Saved as .npz file at the specified output path

## Dependencies

- **Newick Tree File**: Ensure that all Newick Tree files are generated by *Script 4* for all samples before running the script.

## Script

```
python 06_sparse_mtx.py
```

# Step 7: Compute Similarity and Reciprocal Best Hits (RBH) Between Node Pairs

## Summary

This script calculates the Jaccard similarity between gene expression markers of nodes across different datasets and subsets, identifies reciprocal best hits (RBH), and outputs the results to CSV files. It uses multiprocessing to handle the large number of node comparisons efficiently.

## Inputs

- **Output Path**: Folder path for output files, e.g., `/nodes/`
- **Cellstates Base Path**: Base path containing all cellstates outputs (including Newick tree files), e.g., `/VisualStudioCode/`
- **Sparse Matrix Directory**: Directory containing the sparse matrix files, e.g., `/sparse_mtx/`
- **Markers Base Path**: Base path for marker files, e.g., `/markers/`
- **Number of Processes**: Number of processes for multiprocessing

## Outputs

- **Graph Nodes Data CSV Files**: CSV files containing node pairs, Jaccard scores, RBH status, and dataset information.

## Dependencies

- **Newick Tree File**: Ensure that all Newick Tree files are generated by *Script 4* for all samples before running the script.
- **Sparse Matrices**: Ensure that all node-cell matrices are generated by *Script 6* and located inside a single directory.

## Script

```
python3 07_similarity_RBH.py --output_path /nodes/ --sparse_matrix_dir /sparse_mtx/ --cellstates_base_path /home/celia/Documents/VisualStudioCode --markers_base_path /markers/ --num_processes 2
```

# Step 8: Combine Pairwise Similarity Data

## Summary

This script combines multiple CSV files containing pairwise similarity data into a single CSV file. It reads all CSV files in the specified input directory, combines their content, and writes the combined data to the specified output file.

## Inputs

- **Input Directory**: Path to the folder containing all samples' pairwise comparisons, e.g., `/data/nodes/`
- **Output File Path**: Path to save the combined output file in CSV format, e.g., `/combined_data.csv`

## Outputs

- **Combined CSV File**: A single CSV file that combines data from all input CSV files, saved at the specified output path

## Dependencies

- **CSV Files**: Ensure that all CSV files generated by *Script 7* are present in the input directory before running the script.

## Script

```
python3 08_combined_pairwise_similarity.py /data/nodes/ /data/combined_data.csv
```

# Step 9: Filter Node Data

## Summary

This script filters nodes based on Reciprocal Best Hits (RBH) and a specified Jaccard score threshold. It reads an input CSV file, applies the filters, and saves the filtered data to a new CSV file.

## Inputs

- **Input File Path**: Path to the input CSV file containing node data, e.g., `/data.csv`
- **Output File Path**: Path to save the filtered output CSV file, e.g., `data_filtered.csv`
- **RBH Filter**: Boolean value to filter by RBH (Reciprocal Best Hits), e.g., `True`
- **Jaccard Score Threshold**: Minimum Jaccard score to keep, e.g., `0.20`

## Outputs

- **Filtered CSV File**: A CSV file containing nodes that meet the specified filters, saved at the specified output path

## Dependencies

- **Data File**: Ensure that the input data file to be filtered is the combined data file generated by *Script 8*.

## Script

```
python3 09_filter_similarities.py /path/to/input/full_data.csv /path/to/output/full_data_filtered.csv --rbh True --threshold 0.20
```

# Step 10: Perform Hierarchical Clustering on Node Data

## Summary

This script performs hierarchical clustering on node data based on a distance metric (Jaccard score) and saves the filtered results. It reads an input CSV file containing similarity data, creates a distance matrix, performs clustering, filters clusters based on size, and saves the results as a condensed distance matrix and node labels.

## Inputs

- **Input File Path**: Path to the input CSV file containing similarity data, e.g., `/path/to/input/full_data_filtered.csv`
- **Output Directory**: Directory to save the output files, e.g., `/path/to/output_dir`
- **Distance Column**: Column name containing the distance metric in the input file, default is `Jaccard Score`
- **Linkage Method**: Linkage method for hierarchical clustering, default is `average`
- **Threshold**: Threshold for clusters in the dendrogram, default is `0.8`
- **Minimum Cluster Size**: Minimum size of clusters to keep, default is `3`

## Outputs

- **Condensed Matrix File**: A CSV file containing the condensed distance matrix, saved in the specified output directory.
- **Labels File**: A text file containing node labels, saved in the specified output directory.

## Dependencies

- **Data File**: Ensure that the input data file to be filtered is the combined data file generated by *Script 8* or filtered by *Script 9*.

## Script

```
python 10_nodes_clustering.py /path/to/input/full_data.csv /path/to/output_dir --distance_column 'Jaccard Score' --linkage_method 'average' --threshold 0.8 --min_cluster_size 3
```
# Script 11: Generate Pairwise Matrix Comparing Clusters

## Summary

This script  generates a pairwise matrix comparing clusters based on their cell compositions. It reads input files including a condensed matrix, labels, and annotation data, and calculates pairwise similarities between clusters based on cell overlap. The results are saved as a CSV file.

## Inputs

- **Condensed Matrix Path**: Path to the condensed matrix CSV file, e.g., `/condensed_matrix.csv`
- **Labels Path**: Path to the labels CSV file, e.g., `/labels.csv`
- **Output Path**: Path to save the output pairwise matrix CSV file, e.g., `/clusters_comparisons.csv`
- **Sparse Matrix Directory**: Directory containing the sparse matrix files, e.g., `/sparse_mtx/`

## Outputs

- **Pairwise Matrix**: CSV file containing pairwise comparisons of clusters based on cell overlap.

## Dependencies

- **Condensed Matrix & labels**: These two files are generated by *Script 10*
- **Sparse Matrices**: Ensure that all node-cell matrices are generated by *Script 6* and located inside a single directory.

## Script

```
python3 11_comparing_clusters.py /condensed_matrix.csv /labels.csv /clusters_comparisons.csv /sparse_mtx/
```

# Script 12: Interactive Graph for Clustering Analysis

## Summary

This script generates an interactive sunburst chart. It reads input files including a condensed matrix, labels, and annotation data, and visualizes clusters with hierarchical relationships. Users can specify an annotation to color the sunburst chart and visualize different attributes of the clusters.

## Inputs

- **Condensed Matrix Path**: Path to the condensed matrix CSV file, e.g., `/condensed_matrix.csv`
- **Labels Path**: Path to the labels CSV file, e.g., `/labels.csv`
- **Node Annotation Path**: Optional. Path to the annotation CSV file containing node annotations, e.g., `/annotations.csv`
- **Cluster Comparison Path**: Path to the pairwise matrix CSV file comparing clusters, e.g., `/cluster_comparison.csv`
- **Sparse Matrix Directory**: Directory containing the sparse matrix files, e.g., `/sparse_mtx/`
- **Annotation**: Optional. Annotation column name to color the sunburst chart based on cluster attributes.

## Outputs

- **Interactive Sunburst Chart**: Visualization of hierarchical clustering results with interactive features for exploring cluster relationships and attributes.

## Dependencies

- **Condensed Matrix & labels**: These two files are generated by *Script 10*
- **Sparse Matrices**: Ensure that all node-cell matrices are generated by *Script 6* and located inside a single directory.


## Script

```
python3 12_interactive_graph.py /condensed_matrices/condensed_matrix.csv /condensed_matrices/labels.csv /annotation.csv /clusters_comparisons.csv --annotation Epithelial
```
## Authors

- [Usama AKHTAR](https://github.com/usama-ak)
