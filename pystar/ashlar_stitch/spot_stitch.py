import os
import numpy as np
import pandas as pd
from scipy.spatial import distance

def rotate_coords(x, y, w, h, rotation):
    if rotation == 0:
        return x, y
    elif rotation == 90:
        return h - y, x
    elif rotation == 180:
        return w - x, h - y
    elif rotation == 270:
        return y, w - x
    return x, y

def empty_reads():
    return pd.DataFrame({
        'spot_location_1': pd.Series(dtype=float),
        'spot_location_2': pd.Series(dtype=float),
        'spot_location_3': pd.Series(dtype=float),
        'cell_barcode': pd.Series(dtype=int),
        'is_noise': pd.Series(dtype=int),
    })

def empty_cells():
    return pd.DataFrame({
        'cell_barcode': pd.Series(dtype=int),
        'column': pd.Series(dtype=float),
        'row': pd.Series(dtype=float),
        'z': pd.Series(dtype=float),
    })

def read_reads_csv(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return empty_reads()
    reads = pd.read_csv(path)
    if reads.empty:
        return empty_reads()
    if 'clustermap' in reads.columns:
        reads.rename(columns={'clustermap': 'cell_barcode'}, inplace=True)
    if 'x' in reads.columns:
        reads.rename(columns={'x': 'spot_location_1', 'y': 'spot_location_2', 'z': 'spot_location_3'}, inplace=True)
    for col in ('spot_location_1', 'spot_location_2', 'spot_location_3'):
        if col not in reads.columns:
            reads[col] = pd.Series(dtype=float)
    if 'cell_barcode' not in reads.columns:
        reads['cell_barcode'] = -1
    if 'is_noise' not in reads.columns:
        reads['is_noise'] = 0
    return reads

def read_cells_csv(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return empty_cells()
    cells = pd.read_csv(path)
    if cells.empty:
        return empty_cells()
    if 'x' in cells.columns:
        cells.rename(columns={'x': 'column', 'y': 'row'}, inplace=True)
    if 'z_axis' in cells.columns:
        cells.rename(columns={'z_axis': 'z'}, inplace=True)
    for col in ('column', 'row'):
        if col not in cells.columns:
            cells[col] = pd.Series(dtype=float)
    if 'cell_barcode' not in cells.columns:
        cells['cell_barcode'] = pd.Series(dtype=int)
    if 'z' not in cells.columns:
        cells['z'] = pd.Series(dtype=float)
    return cells

def stitch_spots(positions, tile_size, output_dir, pos_coords, spot_rotation):
    H, W = tile_size
    
    if spot_rotation in [90, 270]:
        orig_w, orig_h = H, W
    else:
        orig_w, orig_h = W, H

    centers_y = positions[:, 0] + (H / 2.0)
    centers_x = positions[:, 1] + (W / 2.0)
    fov_centers = np.column_stack((centers_x, centers_y))
    
    all_reads = []
    all_cells = []
    cell_barcode_offset = 0
    
    for i, entry in enumerate(pos_coords):
        pos_path = entry[0]
        pos_name = os.path.basename(pos_path)
        
        # Directly construct paths based on the input directory structure
        reads_file = os.path.join(pos_path, 'remain_reads_raw.csv')
        cells_file = os.path.join(pos_path, 'cell_center.csv')
        
        reads = read_reads_csv(reads_file)
        cells = read_cells_csv(cells_file)
            
        reads['raw_cell_barcode'] = reads['cell_barcode']
        cells['raw_cell_barcode'] = cells['cell_barcode']
        
        rx, ry = rotate_coords(
            reads['spot_location_1'].values,
            reads['spot_location_2'].values,
            orig_w,
            orig_h,
            spot_rotation,
        )
        reads['spot_location_1'], reads['spot_location_2'] = rx, ry
        
        cx, cy = rotate_coords(
            cells['column'].values,
            cells['row'].values,
            orig_w,
            orig_h,
            spot_rotation,
        )
        cells['column'], cells['row'] = cx, cy
        
        global_top_y, global_top_x = positions[i]
        reads['spot_location_1'] += global_top_x
        reads['spot_location_2'] += global_top_y
        cells['column'] += global_top_x
        cells['row'] += global_top_y
        
        cell_coords = cells[['column', 'row']].values
        if len(cell_coords) > 0:
            dists = distance.cdist(cell_coords, fov_centers)
            closest_fov_idx = np.argmin(dists, axis=1)
            cells = cells[closest_fov_idx == i]
        
        valid_cell_barcodes = set(cells['raw_cell_barcode'].values)
        
        if 'is_noise' in reads.columns:
            reads = reads[reads['is_noise'] != -1]
            
        process_mask = reads['raw_cell_barcode'] == -1
        cell_mask = reads['raw_cell_barcode'].isin(valid_cell_barcodes)
        
        reads_process = reads[process_mask]
        if len(reads_process) > 0:
            rp_coords = reads_process[['spot_location_1', 'spot_location_2']].values
            rp_dists = distance.cdist(rp_coords, fov_centers)
            rp_closest = np.argmin(rp_dists, axis=1)
            reads_process = reads_process[rp_closest == i]
            
        reads_within_cells = reads[cell_mask & (~process_mask)]
        reads = pd.concat([reads_within_cells, reads_process])
        
        process_mask_final = reads['raw_cell_barcode'] == -1
        reads.loc[~process_mask_final, 'cell_barcode'] += cell_barcode_offset
        reads.loc[process_mask_final, 'cell_barcode'] = 0
        cells['cell_barcode'] += cell_barcode_offset
        
        if len(cells) > 0:
            cell_barcode_offset = cells['cell_barcode'].max()
            
        reads['FOV'] = pos_name
        cells['FOV'] = pos_name
        all_reads.append(reads)
        all_cells.append(cells)
    
    final_reads = pd.concat(all_reads, ignore_index=True) if all_reads else empty_reads()
    final_cells = pd.concat(all_cells, ignore_index=True) if all_cells else empty_cells()
    
    final_reads.rename(columns={'spot_location_1': 'column', 'spot_location_2': 'row', 'spot_location_3': 'z'}, inplace=True)
    final_reads.drop_duplicates(inplace=True)
    
    final_cells.to_csv(os.path.join(output_dir, 'cell_centerouter.csv'))
    final_reads.to_csv(os.path.join(output_dir, 'remain_readsouter.csv'))
    print(f"    Saved {len(final_cells)} cells and {len(final_reads)} reads to CSVs.")
