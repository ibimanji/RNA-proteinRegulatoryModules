import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors

def get_figsize(max_col, max_row, scale=5):
    # Match the original FIJI script's dynamic sizing logic
    width = (max_col / 100 / scale) * 2
    height = (max_row / 100 / scale) * 2
    
    # Add a safety cap to prevent MemoryError on extremely large datasets
    # while ensuring extremely high resolution
    width = min(width, 150)
    height = min(height, 150)
    return (width, height)

def visualize_spots(output_dir):
    reads_file = os.path.join(output_dir, 'remain_readsouter.csv')
    cells_file = os.path.join(output_dir, 'cell_centerouter.csv')
    coords_file = os.path.join(output_dir, 'registered_tilecoord.csv')
    
    reads = pd.read_csv(reads_file, index_col=0)
    cells = pd.read_csv(cells_file, index_col=0)
    coords = pd.read_csv(coords_file)
    
    # Calculate canvas boundaries based on actual pixel coordinates
    # Adding 3000 buffer to approximate the width/height of the edge FOVs
    max_col = coords['Global_X_px'].max() + 3000  
    max_row = coords['Global_Y_px'].max() + 3000
    shape_row = max_row
    
    # Generate exactly the same random colormap as the FIJI script
    cmap_random = matplotlib.colors.ListedColormap(np.random.rand(256, 3))
    
    # Extract numeric FOV IDs for scatter coloring and labeling
    fov_numbers = []
    for fov_str in coords['FOV']:
        match = re.search(r'\d+', str(fov_str))
        fov_numbers.append(int(match.group()) if match else 0)
        
    plt.style.use('dark_background')
    
    # Use dynamic figsize to ensure dot sizes (s=0.1, s=10) render proportionally
    fig_size = get_figsize(max_col, max_row, scale=5)
    plt.figure(figsize=fig_size)
    
    # Reads with cell centers
    plt.subplot(2, 2, 1)
    plt.title('Reads with cell centers', fontsize=30)
    plt.scatter(reads['column'], shape_row - reads['row'], s=0.1, alpha=0.8, 
                c=pd.Categorical(np.array(reads['raw_cell_barcode'])).codes, cmap=cmap_random)
    plt.scatter(cells['column'], shape_row - cells['row'], s=1, c='red', alpha=1)
    plt.axis('off')
    
    # Reads with cell centers and tile order
    plt.subplot(2, 2, 2)
    plt.title('Reads with cell centers and tile order', fontsize=30)
    plt.scatter(reads['column'], shape_row - reads['row'], s=0.1, alpha=0.8, 
                c=pd.Categorical(np.array(reads['raw_cell_barcode'])).codes, cmap=cmap_random)
    plt.scatter(cells['column'], shape_row - cells['row'], s=1, c='red', alpha=1)
    
    # Plot FOV anchor points and text directly at the top-left coordinates without offset
    plt.scatter(x=coords['Global_X_px'], y=shape_row - coords['Global_Y_px'], c=fov_numbers, s=100)
    for i, row in coords.iterrows():
        plt.text(x=row['Global_X_px'], y=shape_row - row['Global_Y_px'], s=str(fov_numbers[i]), 
                 fontdict=dict(fontsize=20, color='white'))
    plt.axis('off')
    
    # Cell centers
    plt.subplot(2, 2, 3)
    plt.title('Cell centers', fontsize=30)
    plt.scatter(cells['column'], shape_row - cells['row'], s=10, c='red', alpha=0.8)
    plt.axis('off')
    
    # Cell centers and tile order
    plt.subplot(2, 2, 4)
    plt.title('Cell centers and tile order', fontsize=30)
    plt.scatter(cells['column'], shape_row - cells['row'], s=10, c='red', alpha=0.8)
    
    # Plot FOV anchor points and text directly at the top-left coordinates without offset
    plt.scatter(x=coords['Global_X_px'], y=shape_row - coords['Global_Y_px'], c=fov_numbers, s=100)
    for i, row in coords.iterrows():
        plt.text(x=row['Global_X_px'], y=shape_row - row['Global_Y_px'], s=str(fov_numbers[i]), 
                 fontdict=dict(fontsize=20, color='white'))
    plt.axis('off')
    
    plt.tight_layout()
    out_img = os.path.join(output_dir, 'cell_reads_profile.png')
    
    plt.savefig(out_img)
    print("    Saved diagnostic dashboard -> cell_reads_profile.png")