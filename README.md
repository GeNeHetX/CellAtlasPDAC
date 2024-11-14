## Authors
- [Usama Akhtar](https://github.com/usama-ak)
- 
# Introduction

This document outlines the process of creating a comprehensive cell atlas using a series of scripts. Each script performs a specific step in the workflow, with clearly defined inputs and outputs. The scripts are designed to process single-cell RNA-seq data and build a hierarchical cell state atlas.

## Workflow Overview

The workflow is divided into two main phases:

1. **Dataset Preparation (Scripts 1 to 6-2)**:
    - These scripts must be executed **independently for each dataset**. 
    - This phase involves data preprocessing, quality control, clustering with the `cellstates` package, and marker identification.
    - It is mandatory to complete **all 6 scripts** for every new dataset before proceeding to the next phase.

2. **Atlas Creation (Scripts 7 to 11)**:
    - Once the initial preparation is done for all datasets using scripts 1 to 6-2, you can proceed with scripts 7 to 11 to create a unified cell atlas.
    - **Important**: The output directory specified in **Script 7** should be the same as the one used when creating the cell atlas to ensure consistency across datasets.

## Guidelines for Adding a New Dataset

- **Step 1**: For any new dataset, start by executing scripts 1 through 6-2 sequentially to process and analyze the dataset independently.
- **Step 2**: After completing scripts 1 to 6-2 for all your datasets, move on to scripts 7 through 11 to add the new generated data with the existing data to generate the final cell atlas.

## Repository Information

- **GitHub Repository**: https://github.com/GeNeHetX/CellAtlasPDAC
- **Version**: v1.2.0

Please refer to the detailed instructions provided in each script's section for further guidance on how to set up the inputs, run the analysis, and interpret the outputs.


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

# Step 2: Split RNA Matrix For Each Sample and Directory Creation

## Summary

This script filters the RNA matrix based on sample IDs extracted from column names, and saves the filtered data **for each sample** into separate files (`RNAmatrix_{sample_id}.tsv`), and it creates directories named after each sample ID to store eventual *cellstates* outputs. Additionally, it generates a .txt file that lists all sample IDs.

## Inputs

- **RNA Matrix**: Path to a RNA matrix of a dataset, e.g., `/cellstatesappli/RNAmatrix.tsv`.
- **Sample Separator**: Character that separates sample identifiers in the column names of the RNA matrix (e.g., `_`, `-`)
- **Sample Position**: Position of the sample identifier in the column name (e.g., 0 for `sample1_cellbarcode1` if the separator is `_`).
- **Samples File Directory**: Directory where the output file containing the list of sample IDs will be saved, e.g., `/home/celia/Documents/VisualStudioCode/samples/`.

## Outputs

- **Filtered RNA Matrices**: Saved as `RNAmatrix_{sample_id}.tsv` **for each sample ID** in the directory `/cellstatesappli/`. The filtering process removes genes with zero expression across all cells and cells with zero expression across all genes.
- **Directories**: Directories named `RNAmatrix_{sample_id}` are created in `/cellstatesappli/` to store eventual *cellstates* outputs.
- **Sample ID List File**: A .txt file containing all extracted sample IDs is created and saved as {dataset}_samples.txt in the sample output directory, e.g., `/home/celia/Documents/VisualStudioCode/samples/`.

## Script: 02_sample_matrix.py

```
python3 02_sample_matrix.py
```

# Step 3: Running Cell States Analysis

## Summary

This script runs a Python script (`run_cellstates.py`) for each sample listed in `ITEMS`. It specifies input file paths (`INPUT_FILE`) and output directories (`OUTPUT_DIR`) for each sample.  

Cellstates is a python package for analysis of UMI-based single-cell RNA-seq data. It describes the higher-order relationship of these cell-states in a hierarchical tree and provide scores for marker-genes within this tree. For more information, visit the CellStates GitHub repository : https://github.com/nimwegenLab/cellstates.git.

## Dependencies

- **Cellstates Python Package**: Ensure that you install the `cellstates` package before running the script.
    - **Installation Instructions**:  
        ```  
        cd /Documents/GitHub/
        git clone https://github.com/nimwegenLab/cellstates.git  
        ```

- **CellStates Environment**: Ensure the `cellstates-env` Conda environment is activated before running the script.
    - **Conda Environment Set Up Instructions**:  
        ```  
        # Create Environment
        conda create -n cellstates-env python=3.11  
        conda activate cellstates-env  

        # Instal python packages
        conda install numpy scipy cython pandas matplotlib  
        pip install scanpy  
        pip install --upgrade ete3  
        conda install -c conda-forge pypickle  

        # Build and install the CellStates package
        cd /Documents/GitHub/cellstates/  
        python setup.py build_ext --inplace  
        python setup.py install  


        ```

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

This Python script processes RNA matrices, and Cellstates outputs (cluster assignments, Dirichlet pseudocounts, and cluster hierarchies) for multiple subsets within a specified dataset. It generates Newick tree strings for each subset and writes them to output files.

## Dependencies

- **CellStates Environment**: Ensure the `cellstates-env` Conda environment is activated before running the script.

```bash
conda activate cellstates-env
```

## Inputs

- **Base Path**: Base path containing all cellstates outputs
- **Dataset Name**: Name of the folder containing cellstates outputs of a specific dataset
- **Sample Directory**: Path to a directory containing sample lists, e.g., `/home/celia/Documents/VisualStudioCode/samples/`
    - For each dataset, there should be a corresponding text file with sample names, e.g., `peng_samples.txt`
        - Content of `peng_samples.txt`:

            ```  
            T2  
            T3  
            T6  
            T7  
            ```  

## Outputs

- **Newick Tree Strings**: Generated and saved to `newick_tree.txt` files within respective subset directories.

## Script

```
python 04_generate_newick.py
```


# Step 5: Identify Markers with DGEA

## Summary

This R script processes RNA matrices, cell annotations, and cluster data for multiple subsets within a specified dataset. It processes nodes for each subset, and identifies marker genes. The results are saved to CSV files.

## Inputs

- **Cell Annotation Path**: Path to file containing the cell annotation data or cell barcodes of a dataset, e.g., `/Zhang.Cells.Annotation.V1.tsv`
- **Expression Matrix Path**: Path to file containing the expression matrix of a dataset in 'mtx' format, e.g.,`/Zhang.rawcount.V1.mtx`
- **Features Path**: Path to file containing gene annotation data or gene list of a dataset, e.g., `/Zhang.geneAnnotation.V1.tsv`
- **Subset Matrix Base Path**: Path to folder containing matrices of all subsets of a dataset
- **Output Path**: Path to output folder
- **Sample File**: Path to a file containing samples of the dataset, e.g., `/home/celia/Documents/VisualStudioCode/samples/werba_samples.txt`

## Outputs

- **Markers CSV Files**: Marker genes for each node in each subset, saved in the output directory.

### Example of Output CSV file**:  

|        | p_val             | avg_log2FC | pct.1 | pct.2 | p_val_adj        |
|--------|-------------------|------------|-------|-------|------------------|
| VENTX  | 3.38E-35         | 9.07       | 0.25  | 0.001 | 6.18E-31         |
| ZNF843 | 3.38E-35         | 9.20       | 0.25  | 0.001 | 6.18E-31         |
| DDIT4L | 4.73E-35         | 7.26       | 0.50  | 0.005 | 8.66E-31         |
| ST8SIA1| 3.11E-32         | 6.53       | 0.50  | 0.006 | 5.68E-28         |
| LBP    | 1.08E-26         | 8.42       | 0.25  | 0.002 | 1.97E-22         |
| GPC6   | 2.29E-22         | 5.69       | 0.50  | 0.009 | 4.19E-18         |

### Explanation of Columns:
- **Gene**: The name of the gene identified as a marker.
- **p_val**: The p-value of the differential expression test.
- **avg_log2FC**: Average log2 fold-change between groups.
- **pct.1**: Percentage of cells in the first group expressing the gene.
- **pct.2**: Percentage of cells in the second group expressing the gene.
- **p_val_adj**: Adjusted p-value for multiple testing.

## Dependencies

- **Newick Tree File**: Ensure that all Newick Tree files are generated by *Script 4* for all samples before running the script.
- **Presto library**: To optimize performance, it is highly recommended to install the `presto` library. This library accelerates the computation, significantly reducing processing time, especially when handling large datasets. You can install it in R with the following command:  
```
# install.packages("devtools")  
devtools::install_github("immunogenomics/presto")
```

## Script

```
Rscript 05_dgea.r
```

# Step 6-1: Create Node-Cell Matrix

## Summary

This script generates a node-cell matrix using node IDs from Newick tree files and RNA matrix data for multiple subsets. It utilizes multiprocessing to efficiently process nodes and outputs a sparse matrix in .npz format.

## Inputs

- **Base Path**: Base path containing all cellstates outputs, e.g., `/home/celia/Documents/VisualStudioCode/`
- **Sample Directory**: Path to a directory containing sample lists, e.g., `/home/celia/Documents/VisualStudioCode/samples/`
    - For each dataset, there should be a corresponding text file with sample names, e.g., `peng_samples.txt`
- **Dataset Name**: Name of the dataset for which to generate the matrix e.g., `peng`
- **Output File Path**: Path where the resulting .npz sparse matrix will be saved, e.g., `/media/celia/data/cell_mtx/cell_mtx_peng_data.npz`
- **Number of Processes**: Number of processes for multiprocessing

## Outputs

- **Sparse Node-Cell Matrix**: Saved as .npz file at the specified output path

## Dependencies

- **Newick Tree File**: Ensure that all Newick Tree files are generated by *Script 4* for all samples before running the script.

## Script

```
python 06_1_cell_mtx.py
```

# Step 6-2: Create Node-Gene Matrix

## Summary

This script generates a Node-Gene Matrix by processing differential gene expression analysis (DGEA) output files for multiple datasets. The matrix captures over-expressed and under-expressed genes for each node, creating sparse matrices stored in compressed .npz format. The script leverages multiprocessing to parallelize processing across multiple samples for efficiency.

## Inputs

- **DGEA Output Path**: Base path containing differential gene expression analysis output files for various datasets, e.g., `/home/celia/Documents/annotation/findmarkers/`
- **Sample Directory**: Path to a directory containing sample lists, e.g., `/home/celia/Documents/VisualStudioCode/samples/`
    - For each dataset, there should be a corresponding text file with sample names, e.g., `peng_samples.txt`
- **Datasets Name**: List of dataset names to process, e.g., `['peng', 'steele', 'lin', 'zhang', 'hwang', 'werba']`
- **Output Directory**: Path to save the resulting sparse matrices, e.g., `/media/celia/data/gene_mtx/`
- **Number of Processes**: Number of processes for multiprocessing (default set to 18 in the script)

## Outputs

- **Sparse Node-Gene Matrix**:
    - For each dataset and sample, two sparse matrices (over-expressed and under-expressed genes) are generated and saved in a compressed .npz format. The output includes:
        - `sparse_matrix_over`: Sparse matrix of over-expressed genes.
        - `sparse_matrix_under`: Sparse matrix of under-expressed genes.
        - `gene_index_map_over`: Gene index mapping for over-expressed genes.
        - `gene_index_map_under`: Gene index mapping for under-expressed genes.
        - `node_index_map`: Mapping of nodes to matrix rows.
    
    - **Example Output Files**
        ```
        /media/celia/data/gene_mtx/peng/
        ├── peng_Sample1_data.npz
        └── peng_Sample2_data.npz
        ```

## Dependencies

- **DGEA CSV Files**: Ensure that each dataset's differential gene expression analysis outputs are organized in folders by sample names. 
    - **Example Directory Structure**

        ```
        /home/celia/Documents/annotation/findmarkers/  
        ├── peng/  
        │   ├── Sample1/  
        │   │   ├── NodeA.csv  
        │   │   ├── NodeB.csv  
        │   │   └── NodeC.csv  
        │   └── Sample2/  
        │       ├── NodeA.csv  
        │       ├── NodeB.csv  
        │       └── NodeC.csv  
        └── steele/  
            ├── Sample1/  
            ├── Sample2/  
        ```

## Script

```
python 06_2_gene_mtx.py
```

# Step 7: Compute Similarity and Reciprocal Best Hits (RBH) Between Node Pairs

## Summary

This script calculates the Jaccard similarity between gene expression markers of nodes across different datasets and subsets, identifies reciprocal best hits (RBH), and outputs the results to CSV files. It uses multiprocessing to handle the large number of node comparisons efficiently.

## Inputs

- **Output Path**: Path to the directory where the output files will be saved, e.g., `/media/celia/data/nodes/`
- **Node Cell Matrix Directory**: Directory containing the sparse cell matrix files for each dataset., e.g., `/media/celia/data/cell_mtx/`
- **Node Gene Matrix Directory**: Directory containing the sparse gene matrix files for each dataset, e.g., `/media/celia/data/gene_mtx/`
- **Sample Directory**: Path to a directory containing sample lists, e.g., `/home/celia/Documents/VisualStudioCode/samples/`
- **Number of Processes**: Number of processes for multiprocessing

## Outputs

- **Graph Nodes Data CSV Files**: CSV files containing node pairs, Jaccard scores, RBH status, and dataset information.

## Dependencies

- **Node Cell Sparse Matrices**: Ensure that all node-cell matrices are generated by *Script 6-1*.
- **Node Gene Sparse Matrices**: Ensure that all node-genes matrices are generated by *Script 6-2*.

## Script

```
python3 07_similarity_rbh.py \
  --output_path <output_folder_path> \
  --cell_mtx_dir <cell_matrix_directory> \
  --gene_mtx_dir <gene_matrix_directory> \
  --samples_dir <samples_directory> \
  --num_processes <num_processes>
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
python3 08_combined_pairwise_similarity.py \
    --input_dir /data/nodes/ \
    --output_path /data/combined_data.csv
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
# Step 11: Establishing hierarchy between the clusters

## Summary 

This Jupyter notebook focuses on refining and organizing the clusters obtained from hierarchical clustering (*Step 10*). The goal is to establish a hierarchy among clusters, identifying parent-child relationships based on node connectivity and similarity.

## Key Steps

1) Cluster Cleaning:

- For each cluster, nodes are analyzed to identify the most representative (highest) nodes. This is achieved by filtering out redundant nodes within the same branch, ensuring only the top nodes are kept for further analysis.

2) Cluster Comparison:

- A comprehensive comparison is performed between all clusters to evaluate their hierarchical relationships. This involves counting how frequently nodes from one cluster appear above or below nodes from another cluster. This comparison helps identify potential parent-child relationships between clusters.

3) Flow-Based Hierarchical Assignment:

- To establish a clear hierarchy, a flow-based approach is used to determine parent clusters. By modeling the clusters as a directed graph, we compute the maximum flow between clusters to assign parent-child relationships. This method ensures that each cluster is linked to its most relevant parent, creating a robust hierarchy.

4) Hierarchy Visualization:

- The final step involves visualizing the established hierarchy using a graph with a hierarchical layout. This visualization helps in understanding the structure and relationships between clusters, making it easier to interpret the clustering results.

## Inputs

- Clustering Output: Path to the CSV file containing the clustering results obtained from the previous script (*Script 10*). This file serves as the input for establishing hierarchical relationships between clusters.

## Jupyter Notebook

```
evaluate_clusters.ipynb
```

