# conda activate cellstates-env

import cellstates as cs
import pandas as pd
import numpy as np
from pathlib import Path

base_path = '/home/celia/Documents/VisualStudioCode/'
sample_dir = Path('/home/celia/Documents/VisualStudioCode/samples/')
dataset = 'STEELE'

file_path = sample_dir / f'{dataset.lower()}_samples.txt'
if file_path.exists():
    with open(file_path, 'r') as f:
        subsets = [line.strip() for line in f.readlines()]
else:
    print(f"Warning: File {file_path} not found!")

for subset in subsets:
    clusters_path = base_path + f"{dataset}/cellstatesappli/RNAmatrix_{subset}/optimized_clusters.txt"
    lmbd_path = base_path + f"{dataset}/cellstatesappli/RNAmatrix_{subset}/dirichlet_pseudocounts.txt"
    hierarchy_df_path = base_path + f"{dataset}/cellstatesappli/RNAmatrix_{subset}/cluster_hierarchy.tsv"
    RNA_matrix_path = base_path + f"{dataset}/cellstatesappli/RNAmatrix_{subset}.tsv"

    data_ = pd.read_csv(RNA_matrix_path, sep="\t", index_col=0)
    clusters = np.loadtxt(clusters_path, dtype=np.int64)
    lmbd = np.loadtxt(lmbd_path)
    hierarchy_df = pd.read_csv(hierarchy_df_path, sep='\t')
    data = data_.values.astype(np.int64)
    clst = cs.Cluster(data, lmbd, clusters)
    newick_string = cs.hierarchy_to_newick(hierarchy_df, clst.clusters, cell_leaves=False)

    file_path = base_path + f"{dataset}/cellstatesappli/RNAmatrix_{subset}/newick_tree.txt"

    with open(file_path, "w") as file:
        file.write(newick_string)
    print(f"Newick tree string written to: {file_path}")
