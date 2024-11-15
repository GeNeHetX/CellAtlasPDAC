from multiprocessing import Pool
from functools import partial
from scipy.sparse import lil_matrix
import numpy as np
import pandas as pd
import os
from pathlib import Path
from tqdm import tqdm

def process_sample(datas_key, sample, dgea_output_path, output_dir):
    print(f"Processing sample: {sample}")
    output_path = f'{output_dir}{datas_key}/'
    os.makedirs(output_path, exist_ok=True)
    output_file = os.path.join(output_path, f'{datas_key}_{sample}_data.npz')
       
    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"Output file already exists for sample: {sample}. Skipping processing.")
        return sample

    folder_path = f'{dgea_output_path}{datas_key}/{sample}/'
    
    all_genes_over = set()
    all_genes_under = set()
    node_gene_map = {}

    try :
        # Process each file in the folder
        for filename in os.listdir(folder_path):
            if filename.endswith('.csv'):
                node_name = os.path.splitext(filename)[0]  
                filepath = os.path.join(folder_path, filename)
                
                df = pd.read_csv(filepath)

                if datas_key == 'hwang':
                    df['Unnamed: 0'] = df['Unnamed: 0'].str.split('\t').str[0]
                
                # Filter for over-expressed and under-expressed genes
                over_expressed_genes = df[(df['avg_log2FC'] > 0) & (df['p_val_adj'] <= 0.05)]['Unnamed: 0'].tolist()
                under_expressed_genes = df[(df['avg_log2FC'] < 0) & (df['p_val_adj'] <= 0.05)]['Unnamed: 0'].tolist()

                all_genes_over.update(over_expressed_genes)
                all_genes_under.update(under_expressed_genes)
                
                node_gene_map[node_name] = (over_expressed_genes, under_expressed_genes)

    except Exception as e:
        print(f"Error processing sample {sample}: {e}")
        return None

    # Sorted gene lists and index mappings
    all_genes_over = sorted(list(all_genes_over))
    all_genes_under = sorted(list(all_genes_under))

    gene_index_map_over = {gene: idx for idx, gene in enumerate(all_genes_over)}
    gene_index_map_under = {gene: idx for idx, gene in enumerate(all_genes_under)}

    num_nodes = len(node_gene_map)
    num_genes_over = len(gene_index_map_over)
    num_genes_under = len(gene_index_map_under)

    # Create sparse matrices for over/under expressed genes
    sparse_matrix_over = lil_matrix((num_nodes, num_genes_over), dtype=np.int8)
    sparse_matrix_under = lil_matrix((num_nodes, num_genes_under), dtype=np.int8)

    node_index_map = {node: idx for idx, node in enumerate(node_gene_map.keys())}

    # Fill sparse matrices
    for node, (over_genes, under_genes) in node_gene_map.items():
        row_idx = node_index_map[node]
        for gene in over_genes:
            if gene in gene_index_map_over:
                col_idx = gene_index_map_over[gene]
                sparse_matrix_over[row_idx, col_idx] = 1
        for gene in under_genes:
            if gene in gene_index_map_under:
                col_idx = gene_index_map_under[gene]
                sparse_matrix_under[row_idx, col_idx] = 1

    sparse_matrix_over = sparse_matrix_over.tocsr()
    sparse_matrix_under = sparse_matrix_under.tocsr()
        
    np.savez_compressed(output_file,
                        sparse_matrix_over=sparse_matrix_over,
                        sparse_matrix_under=sparse_matrix_under,
                        gene_index_map_over=np.array(list(gene_index_map_over.items()), dtype=object),
                        gene_index_map_under=np.array(list(gene_index_map_under.items()), dtype=object),
                        node_index_map=np.array(list(node_index_map.items()), dtype=object))

    print(f"Finished processing sample: {sample}")
    return sample


def process_dataset(datas_key, samples, dgea_output_path, output_dir):
    with Pool(18) as pool:
        # Use tqdm for progress tracking
        list(tqdm(pool.starmap(process_sample, [(datas_key, sample, dgea_output_path, output_dir) for sample in samples]), 
                  total=len(samples), desc=f'Processing {datas_key}'))

if __name__ == "__main__":
    dgea_output_path = '/home/celia/Documents/annotation/findmarkers/'
    sample_dir = '/home/celia/Documents/VisualStudioCode/samples/'
    output_dir = '/media/celia/data/gene_mtx/'
    datasets_name = ['peng' , 'steele', 'lin' , 'hwang', 'zhang', 'werba']

    datasets = {}
    for dataset in datasets_name:
        sample_dir = Path(sample_dir)
        file_path = sample_dir / f'{dataset}_samples.txt'
        
        if file_path.exists():
            with open(file_path, 'r') as f:
                samples = [line.strip() for line in f.readlines()]
                datasets[dataset] = samples
        else:
            print(f"Warning: File {file_path} not found!")


    # Process each dataset in sequence
    for datas_key, samples in datasets.items():
        print(f"Processing dataset: {datas_key}")
        process_dataset(datas_key, samples, dgea_output_path, output_dir)
        print(f"Finished processing dataset: {datas_key}")

    print("All datasets processed.")
