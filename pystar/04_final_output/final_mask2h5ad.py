import pandas as pd
import numpy as np
import anndata as ad
import copy
from scipy.sparse import csr_matrix
import os
import time
import tifffile as tiff
from skimage.measure import regionprops
from skimage import measure
import argparse

start_time = time.time()

def parse_args():
    parser = argparse.ArgumentParser(description="Convert stitched mask/read outputs to h5ad and summary CSVs.")
    parser.add_argument("--sample", required=True, help="Sample name used in output h5ad filename.")
    parser.add_argument("--stitch-dir", required=True, help="Directory containing cell_centerouter.csv and remain_readsouter.csv.")
    parser.add_argument("--results-dir", required=True, help="Directory for output CSV and h5ad files.")
    parser.add_argument("--dapi-mask", required=True, help="Path to maskdapiuse_segmentation.tif.")
    parser.add_argument("--output-tag", required=True, help="Tag used in output filenames, e.g. 20251224.")
    return parser.parse_args()

args = parse_args()
sample = args.sample
stitch_dir = args.stitch_dir
results_dir = args.results_dir
dapi_mask_path = args.dapi_mask
output_tag = args.output_tag

# Create the directory if it doesn't exist
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Read the cell center and remaining read data (with somata and process)
cell_center = pd.read_csv(os.path.join(stitch_dir, "cell_centerouter.csv"), index_col=0)
remain_reads = pd.read_csv(os.path.join(stitch_dir, "remain_readsouter.csv"), index_col=0, na_filter=False)

# Process reads: set cell_barcode to -1 for process reads (raw_cell_barcode == -1)
# In new_stitch0220_outerpoints.py, process reads have cell_barcode=0 and raw_cell_barcode=-1
if 'raw_cell_barcode' in remain_reads.columns:
    process_mask = remain_reads['raw_cell_barcode'] == -1
    remain_reads.loc[process_mask, 'cell_barcode'] = -1
    print(f"Process reads (raw_cell_barcode=-1): {process_mask.sum()}")
    print(f"Somata reads: {(~process_mask).sum()}")

# Handle gene column name (could be 'gene' or 'gene_name')
if 'gene_name' not in remain_reads.columns and 'gene' in remain_reads.columns:
    remain_reads['gene_name'] = remain_reads['gene']
    print("Using 'gene' column as 'gene_name'")

dup_cell_barcodes = cell_center['cell_barcode'].duplicated().sum()
if dup_cell_barcodes > 0:
    print(f"Duplicated cell_barcode rows in cell_center: {dup_cell_barcodes}; keeping first occurrence.")
    cell_center = cell_center.drop_duplicates(subset='cell_barcode', keep='first').copy()

cell_center_index = copy.deepcopy(cell_center)
cell_center_index.set_index('cell_barcode', inplace=True, drop=True) 

# Get valid cell_barcodes that exist in cell_center_index (only somata cells)
valid_cell_barcodes = set(cell_center_index.index)

# Create expression matrix: only include somata reads (cell_barcode != -1) 
# AND only include cell_barcodes that exist in cell_center_index
remain_reads_somata = remain_reads[
    (remain_reads['cell_barcode'] != -1) & 
    (remain_reads['cell_barcode'].isin(valid_cell_barcodes))
]
remain_reads_t = remain_reads_somata.loc[:,['cell_barcode','gene_name']]
remain_reads_t['value'] = 1

# Create cell-by-gene expression matrix (only for somata cells that exist in cell_center)
exp_matrix = pd.pivot_table(remain_reads_t, index='cell_barcode', columns='gene_name', aggfunc='count', fill_value=0)
var_raw = [str(s2) for (s1, s2) in exp_matrix.columns.tolist()]
exp_matrix = exp_matrix.set_axis(var_raw, axis=1)

# Create cell metadata (obs) - only for somata cells that exist in cell_center
obs = cell_center_index.loc[exp_matrix.index.values, ['column', 'row', 'z']]

# Create gene metadata (var)
var = pd.DataFrame(index=var_raw)  ## Using gene names as index

# Read DAPI mask (nuclei mask) TIF file
dapi_mask = tiff.imread(dapi_mask_path)

# Function to assign nuclei labels to molecular points (only need the nuclei mask)
def assign_labels_to_points(nuclei_mask, mol):
    """
    Assign labels to molecular points based on the provided nuclei mask
    """
    H, W = nuclei_mask.shape
    # Get the position of molecular points
    x_int = np.clip(mol['column'].values - 1, 0, W - 1).astype(int)  # Adjust based on actual column names
    y_int = np.clip(mol['row'].values - 1, 0, H - 1).astype(int)  # Adjust based on actual row names

    # Check if the molecular points are within the nuclei
    n_vals = nuclei_mask[y_int, x_int]
    mol['nuclei'] = np.where(n_vals > 0, 1, 0)  # Set to 1 if in the nuclei, otherwise set to 0

    return mol

# Use the assign_labels_to_points function to assign nuclei labels to all reads
remain_reads = assign_labels_to_points(dapi_mask, remain_reads)

# For process reads (cell_barcode == -1), set nuclei to -1 (they are not in nuclei)
if 'raw_cell_barcode' in remain_reads.columns:
    process_mask = remain_reads['raw_cell_barcode'] == -1
    remain_reads.loc[process_mask, 'nuclei'] = -1

# Select columns for remain_reads_info, include raw_cell_barcode if available
cols_to_keep = ['gene_name', 'column', 'row', 'z', 'cell_barcode', 'is_noise', 'nuclei']
if 'raw_cell_barcode' in remain_reads.columns:
    cols_to_keep.append('raw_cell_barcode')
remain_reads_selected = remain_reads[cols_to_keep]

# Save the expression matrix (cell-by-gene count matrix) as CSV files
exp_matrix.to_csv(os.path.join(results_dir, f"cell_by_gene_{output_tag}.csv"))
exp_matrix.transpose().to_csv(os.path.join(results_dir, f"gene_by_cell_{output_tag}.csv"))

# Save cell metadata as cell_metadata.csv
obs.to_csv(os.path.join(results_dir, f"cell_metadata_{output_tag}.csv"))

# Save gene metadata as gene_metadata.csv
var.to_csv(os.path.join(results_dir, f"gene_metadata_{output_tag}.csv"))

# Create AnnData object and add the 'nuclei' information to uns
adata = ad.AnnData(
    X=csr_matrix(np.array(exp_matrix)),
    var=var,
    obs=obs
)

# Add the 'nuclei' column to AnnData's uns
adata.uns['remain_reads_info'] = remain_reads_selected
# Save the AnnData object as an h5ad file
adata.write_h5ad(os.path.join(results_dir, f"{sample}_{output_tag}.h5ad"))

end_time = time.time()

execution_time_seconds = end_time - start_time
execution_time_minutes = execution_time_seconds / 60

print(f"Execution time: {execution_time_minutes} minutes")
