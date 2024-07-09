import pandas as pd
import argparse

def filter_nodes(input_path, rbh, threshold):
    """
    Filter nodes based on RBH and Jaccard score threshold.

    Parameters:
    input_path (str): Path to the input CSV file.
    rbh (bool): Filter by RBH (Reciprocal Best Hits) if True.
    threshold (float): Minimum Jaccard score to keep.

    Returns:
    filtered_nodes (set): Set of filtered nodes.
    """
    filtered_nodes = set()
    chunksize = 100000000

    for i, chunk in enumerate(pd.read_csv(input_path, sep='\t', chunksize=chunksize)):
        df_filtered = chunk[(chunk['RBH'] == rbh) & (chunk['Jaccard Score'] >= threshold)]
        filtered_nodes.update(df_filtered['Node1'].tolist())
        filtered_nodes.update(df_filtered['Node2'].tolist())
        print(f'1/2 \t chunk {i} done !')

    return filtered_nodes

def filter_and_save_data(input_path, output_path, filtered_nodes):
    """
    Filter the data to include only rows with filtered nodes and save to a new CSV file.

    Parameters:
    input_path (str): Path to the input CSV file.
    output_path (str): Path to the output CSV file.
    filtered_nodes (set): Set of filtered nodes.
    """
    filtered_nodes_list = list(filtered_nodes)
    chunksize = 100000000
    header_written = False

    for i, chunk in enumerate(pd.read_csv(input_path, sep='\t', chunksize=chunksize)):
        chunk_filtered = chunk[chunk['Node1'].isin(filtered_nodes_list) & chunk['Node2'].isin(filtered_nodes_list)]
        mode = 'w' if not header_written else 'a'
        header = not header_written
        chunk_filtered.to_csv(output_path, sep='\t', index=False, mode=mode, header=header)
        if not header_written:
            header_written = True

        print(f'2/2 \t chunk {i} done !')

def main(input_path, output_path, rbh, threshold):
    """
    Main function to filter and save node data based on RBH and Jaccard score threshold.

    Parameters:
    input_path (str): Path to the input CSV file.
    output_path (str): Path to the output CSV file.
    rbh (bool): Filter by RBH (Reciprocal Best Hits) if True.
    threshold (float): Minimum Jaccard score to keep.
    """
    filtered_nodes = filter_nodes(input_path, rbh, threshold)
    filter_and_save_data(input_path, output_path, filtered_nodes)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter and save node data based on RBH and Jaccard score threshold.')
    parser.add_argument('input_path', type=str, help='Path to the input CSV file.')
    parser.add_argument('output_path', type=str, help='Path to the output CSV file.')
    parser.add_argument('--rbh', type=bool, default=True, help='Filter by RBH (Reciprocal Best Hits) if True.')
    parser.add_argument('--threshold', type=float, default=0.20, help='Minimum Jaccard score to keep.')

    args = parser.parse_args()

    main(args.input_path, args.output_path, args.rbh, args.threshold)

# python script_name.py /path/to/input/full_data.csv /path/to/output/full_data_filtered.csv --rbh True --threshold 0.20
