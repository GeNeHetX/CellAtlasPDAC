import glob
from tqdm import tqdm
import argparse

def combine_csv_files(base_path, output_file):
    with open(output_file, 'w', newline='') as fout:
        header_written = False
        csv_files = glob.glob(str(base_path + '/' + 'graph_nodes_data_*.csv'))
        for csv_file in tqdm(csv_files, desc=f"Processing files in {base_path}", unit="file"):
            with open(csv_file, 'r') as fin:
                if not header_written:
                    header = fin.readline()
                    fout.write(header)
                    header_written = True
                else:
                    fin.readline()
                for line in fin:
                    fout.write(line)

def main(input_dir, output_path):
    combine_csv_files(input_dir, output_path)
    print(f'Combined data saved to {output_path}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform hierarchical clustering on node data and save filtered results.')
    parser.add_argument('input_dir', type=str, help='Path to the folder containing all samples pairwise comparisons')
    parser.add_argument('output_path', type=str, help='Path to save the output file.')
    args = parser.parse_args()
    main(args.input_dir, args.output_path)

# 08_combined_pairwise_similarity.py /media/celia/data/nodes/ /media/celia/data/combined_data.csv
