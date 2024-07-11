import pandas as pd
import plotly.graph_objects as go
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, set_link_color_palette
from ete3 import Tree
from pathlib import Path
import scipy as sp
import argparse

''' INITIALIZE DATA '''

def initialize_data(annotation_path, condensed_matrix_path, labels_path):

    
    condensed_matrix = np.loadtxt(condensed_matrix_path, delimiter=',')
    with open(labels_path, 'r') as f:
        labels = [line.strip() for line in f]

    Z = linkage(condensed_matrix, 'average')
    set_link_color_palette([f"C{i}" for i in range(1, 100000)])
    d = dendrogram(Z, labels=labels, color_threshold=0.8, above_threshold_color='black', no_plot=True)

    trtm_status_hwang = {'U': 'untreated', 'T': 'treated'}
    trtm_status = {'peng': 'untreated', 'steele': 'untreated', 'lin': 'untreated', 'zhang': 'untreated'}

    nodes = []
    clusters = []
    annot1 = []
    annot2 = []
    treatments = []
    datasets = []
    if annotation_path != False:
        colors = pd.read_csv(annotation_path)
        for i, node in enumerate(d['ivl']):
            nodes.append(node)
            clusters.append(d['leaves_color_list'][i])
            annot1.append(colors['Annotation1'][colors['Node'] == node].values[0])
            annot2.append(colors['Annotation2'][colors['Node'] == node].values[0])
            if node.split('-')[-1] == 'hwang':
                trtm = trtm_status_hwang.get(node.split('_')[2][0], 'black')
            else:
                trtm = trtm_status.get(node.split('-')[-1])
            treatments.append(trtm)
            datasets.append(node.split('-')[-1])
    else :
        for i, node in enumerate(d['ivl']):
            nodes.append(node)
            clusters.append(d['leaves_color_list'][i])
            if node.split('-')[-1] == 'hwang':
                trtm = trtm_status_hwang.get(node.split('_')[2][0], 'black')
            else:
                trtm = trtm_status.get(node.split('-')[-1])
            treatments.append(trtm)
            datasets.append(node.split('-')[-1])
    
    nodes_description = pd.DataFrame({
        'Node': nodes,
        'Cluster': clusters,
        'Dataset': datasets,
        'Annotation1': annot1,
        'Annotation2': annot2,
        'Treatment': treatments
    })

    return nodes_description


''' CELL BARCODES'''
datasets = {
    'peng': ['T2', 'T3', 'T6', 'T7', 'T8', 'T9', 'T11', 'T13', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'T21', 'T22'],
    'steele': ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11A", "11B", "12", "13", "15", "16"],
    'lin': ["P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "MET01", "MET02", "MET03", "MET04", "MET05", "MET06"],
    'hwang': ['003_U15', '004_U18', '007_U2', '008_T22', '009_T20', '010T_U12', '011_U9', '2083_T4', '2100_T2', '2223_T17',
              '2229_T18', '2276_U4', '2364_U5', '2376_U1', '2443_U10', '2462_T16', '2490_U11', '2498_U7', '2507_T25', '2523_U3',
              '2540_T6', '2591_U13', '2603_U6', '2626_U16', '2634_T5', '2664_U14', '2667_T23', '2668_T7', '2675_T13', 'MGH2076_T9',
              'MGH2101_T8', 'MGH2381_T19', 'MGH2675_T13', 'MGHR1_T14', 'MGHR2_T15', 'MGHR3_T3', 'MGHR5_T1', 'MGHR6_T10', 'MGHR7_T12', 'MGHR8_T24',
              'MGHR9_T21', 'MGHR11_T11', 'MGHR16_U17', 'MGHR17_U8'],
    'zhang': ['Case1.YF', 'Case1.ZY', 'Case2.YF', 'Case2.ZC', 'Case2.ZY', 'Case3.YF', 'Case3.ZY', 'Case4.ZY']
}

def create_node_ids(subset, dataset):
    tree_path = f"/home/celia/Documents/VisualStudioCode/{dataset.upper()}/cellstatesappli/RNAmatrix_{subset}/newick_tree.txt"
    t = Tree(tree_path, format=1)
    return [f"{node.name}_{subset}-{dataset}" for node in t.traverse() if not node.is_root()]

def create_dataset_node_list(dataset):
    dataset_node_list = []
    for subset in datasets[dataset]:
        dataset_node_list.extend(create_node_ids(subset, dataset))
    return dataset_node_list

def create_dataset_cell_list(dataset):
    data_path = Path(f'/home/celia/Documents/VisualStudioCode/{dataset.upper()}/cellstatesappli/RNAmatrix.tsv')
    data = pd.read_csv(data_path, delimiter='\t', index_col=0, nrows=1)
    return data.columns.tolist()

def get_cells_barcodes(node, sparse_matrices):
    dataset = node.split('-')[-1].lower()
    nodes = node_lists[dataset]
    cells = cell_lists[dataset]
    node_index = nodes.index(node)
    matching_indices = sparse_matrices[dataset].row == node_index
    columns = sparse_matrices[dataset].col[matching_indices]
    corresponding_cells = [cells[col_index] for col_index in columns]
    return corresponding_cells




''' GET PARENTS'''

def get_parents(pairwise_matrix, nodes_description, sparse_matrices) :
    labels = pairwise_matrix.index.tolist()
    parents = []
    for row_index, row in pairwise_matrix.iterrows():
        one_indices = [idx for idx, val in enumerate(row) if val == 1]
        if len(one_indices) == 1:
            parent = row.index[one_indices[0]]
        elif len(one_indices) > 1:
            max_index = pairwise_matrix[row_index].idxmax()
            if row[max_index] == 1:
                parent = max_index
            else:
                parent = row.index[one_indices[0]]
        elif len(one_indices) == 0:
            highest = pairwise_matrix.loc[[row_index], :].idxmax(axis=1).values[0]
            if (pairwise_matrix.loc[[row_index], [highest]].max(axis=1).values[0] > pairwise_matrix.loc[[highest], [row_index]].max(axis=1).values[0]):
                parent = highest

            elif (pairwise_matrix.loc[[row_index], [highest]].max(axis=1).values[0] == pairwise_matrix.loc[[highest], [row_index]].max(axis=1).values[0]):
                clst1 = nodes_description[nodes_description['Cluster'] == row_index]
                cluster1_cells = set()
                for _, row in clst1.iterrows():
                    cell_barcodes = get_cells_barcodes(row['Node'], sparse_matrices)
                    cluster1_cells.update(cell_barcodes)
                
                clst2 = nodes_description[nodes_description['Cluster'] == highest]
                cluster2_cells = set()
                for _, row in clst2.iterrows():
                    cell_barcodes = get_cells_barcodes(row['Node'], sparse_matrices)
                    cluster2_cells.update(cell_barcodes)
                
                intersect = cluster1_cells.intersection(cluster2_cells)
                if len(cluster1_cells) != len(cluster2_cells) :   
                    if (len(intersect)/len(cluster1_cells)) < (len(intersect)/len(cluster2_cells)):
                        parent = row_index
                    else :
                        parent = highest
                elif len(cluster1_cells) == len(cluster2_cells) :
                    second_highest = pairwise_matrix.loc[[row_index], :].drop(columns=pairwise_matrix.loc[[row_index], :].idxmax(axis=1)).idxmax(axis=1).values[0]
                    if (pairwise_matrix.loc[[row_index], [second_highest]].max(axis=1).values[0] > pairwise_matrix.loc[[second_highest], [row_index]].max(axis=1).values[0]):
                        parent = second_highest
                    else : 
                        parent = ''
            else :
                parent = ''
        
        parents.append(parent)
    return labels, parents

def prepare_annotation(nodes_description):
    counts1 = nodes_description.groupby(['Cluster', 'Annotation1']).size().unstack(fill_value=0)
    ordered_index1 = sorted(counts1.index, key=lambda x: int(x[1:]))
    counts1 = counts1.loc[ordered_index1]
    annotation1 = counts1.div(counts1.sum(axis=1), axis=0)

    annotation1['Epithelial'] = annotation1['Epithelial cells'] + annotation1['Epithelial (malignant)'] + annotation1['Ductal cells'] + annotation1['Acinar cell'] + annotation1['Endocrine cells']

    counts2 = nodes_description.groupby(['Cluster', 'Annotation2']).size().unstack(fill_value=0)
    ordered_index2 = sorted(counts2.index, key=lambda x: int(x[1:]))
    counts2 = counts2.loc[ordered_index2]
    annotation2 = counts2.div(counts2.sum(axis=1), axis=0)

    counts3 = nodes_description.groupby(['Cluster', 'Treatment']).size().unstack(fill_value=0)
    ordered_index3 = sorted(counts3.index, key=lambda x: int(x[1:]))
    counts3 = counts3.loc[ordered_index3]
    treatment_annotation = counts3.div(counts3.sum(axis=1), axis=0)

    counts4 = nodes_description.groupby(['Cluster', 'Dataset']).size().unstack(fill_value=0)
    ordered_index4 = sorted(counts4.index, key=lambda x: int(x[1:]))
    counts4 = counts4.loc[ordered_index4]
    dataset_annotation = counts4.div(counts4.sum(axis=1), axis=0)

    return annotation1, annotation2, treatment_annotation, dataset_annotation



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process clustering and generate a sunburst chart.')
    parser.add_argument('condensed_matrix_path', type=str, help='Path to the condensed matrix CSV file.')
    parser.add_argument('labels_path', type=str, help='Path to the labels CSV file.')
    parser.add_argument('node_annotation_path', type=str, help='Path to the node annotation CSV file.')
    parser.add_argument('cluster_comparison_path', type=str, help='Path to the pairwise matrix CSV file comparing clusters.')
    parser.add_argument("--sparse_matrix_dir", type=str, required=True, help="Directory containing the sparse matrix files")
    parser.add_argument('--annotation', type=str, help='Annotation to color the sunburst chart.')
    args = parser.parse_args()

    # Initialize node lists
    node_lists = {}
    for dataset in datasets:
        node_lists[dataset] = create_dataset_node_list(dataset)

    # Initialize cell lists
    data = {'peng', 'steele', 'lin', 'hwang', 'zhang'}
    cell_lists = {}
    for dataset in data:
        cell_lists[dataset] = create_dataset_cell_list(dataset)

    if args.node_annotation_path is not None:
        nodes_description = initialize_data(annotation_path=args.node_annotation_path, condensed_matrix_path=args.condensed_matrix_path, labels_path=args.labels_path)
    else : 
        nodes_description = initialize_data(annotation_path=False, condensed_matrix_path=args.condensed_matrix_path, labels_path=args.labels_path)

    annotation1, annotation2, treatment_annotation, dataset_annotation = prepare_annotation(nodes_description)

    if args.annotation in annotation1.columns:
        colors = annotation1[args.annotation]
    elif args.annotation in annotation2.columns:
        colors = annotation2[args.annotation]
    elif args.annotation in treatment_annotation.columns:
        colors = treatment_annotation[args.annotation]
    elif args.annotation in dataset_annotation.columns:
        colors = dataset_annotation[args.annotation]
    else:
        colors = None
        print(f"Annotation '{args.annotation}' not found in any of the dataframes. Default to None")

    sparse_matrix_dir = Path(args.sparse_matrix_dir)
    datasets = ['peng', 'steele', 'lin', 'hwang', 'zhang']
    
    # Load sparse matrices
    sparse_matrices = {}
    for dataset in datasets:
        sparse_matrix_path = sparse_matrix_dir / f'cell_mtx_{dataset}.npz'
        if sparse_matrix_path.exists():
            sparse_matrices[dataset] = sp.load_npz(sparse_matrix_path)
        else:
            print(f"Warning: Sparse matrix file not found for dataset {dataset} at {sparse_matrix_path}")

    pairwise_matrix = pd.read_csv(args.cluster_comparison_path, index_col=0)
    labels, parents = get_parents(pairwise_matrix, nodes_description, sparse_matrices)

    fig = go.Figure(go.Sunburst(
        labels=labels,
        parents=parents,
    ))

    if colors is not None:
        fig.update_traces(marker=dict(
            colorscale='Bluered',
            colors=colors,
            colorbar=dict(title=''),
        ))

    fig.update_layout(margin=dict(t=0, l=0, r=0, b=0), height=1000)

    fig.show()


# annotation_path = '/media/celia/data/colors.csv'
# condensed_matrix_path = '/media/celia/data/condensed_matrices/condensed_matrix_zhang_RBH02_cut08.csv'
# labels_path = '/media/celia/data/condensed_matrices/labels_zhang_RBH02_cut08.csv'
# pairwise_matrix = pd.read_csv('/media/celia/data/condensed_matrices/clusters_comparisons2.csv', index_col=0)

# python3 12_interactive_graph.py /media/celia/data/condensed_matrices/condensed_matrix_zhang_RBH02_cut08.csv /media/celia/data/condensed_matrices/labels_zhang_RBH02_cut08.csv /media/celia/data/colors.csv /media/celia/data/condensed_matrices/clusters_comparisons2.csv --annotation Epithelial
