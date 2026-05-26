#!/usr/bin/env python3
import argparse
import re
import os
import pandas as pd
from pathlib import Path
from read_maf import get_theoretical_coords
from align import align_reference
from apply_stitch import stitch_all_channels
from local_overlay import write_local_overlays
from spot_stitch import stitch_spots
from visualize import visualize_spots

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', str(s))]

def position_name(pos_id):
    return f"Position{int(pos_id):03d}"

def resolve_position_dir(input_dir, pos_id):
    root = Path(input_dir).resolve()
    padded = root / position_name(pos_id)
    if padded.is_dir():
        return padded
    unpadded = root / f"Position{int(pos_id)}"
    if unpadded.is_dir():
        return unpadded
    return padded

def resolve_image_path(image_dir, spot_dir, pos_id, ref_channel):
    if image_dir:
        root = Path(image_dir).resolve()
        candidates = [
            root / f"{position_name(pos_id)}.tif",
            root / f"Position{int(pos_id)}.tif",
        ]
        for p in candidates:
            if p.is_file():
                return p
        return candidates[0]

    if ref_channel:
        candidates = sorted(spot_dir.glob(f"*{ref_channel}*.tif"))
    else:
        candidates = sorted(spot_dir.glob("*.tif"))
    return candidates[0] if candidates else None

def choose_position_ids(input_dir, maf_coords_dict, maf_start, maf_end):
    if maf_start is not None or maf_end is not None:
        if maf_start is None or maf_end is None:
            raise ValueError("--maf_start and --maf_end must be provided together.")
        if maf_start > maf_end:
            raise ValueError(f"--maf_start must be <= --maf_end, got {maf_start}>{maf_end}")
        ids = list(range(int(maf_start), int(maf_end) + 1))
        missing = [pid for pid in ids if pid not in maf_coords_dict]
        if missing:
            raise ValueError(f"MAF missing PositionID values in requested range: {missing[:20]}")
        return ids

    root = Path(input_dir).resolve()
    pos_dirs = [d for d in root.iterdir() if d.is_dir() and re.fullmatch(r"Position\d+", d.name, re.IGNORECASE)]
    pos_dirs.sort(key=lambda x: natural_sort_key(x.name))
    if not pos_dirs:
        raise ValueError(f"No Position* directories found in {input_dir}!")
    ids = [int(re.search(r"\d+", d.name).group()) for d in pos_dirs]
    missing = [pid for pid in ids if pid not in maf_coords_dict]
    if missing:
        raise ValueError(f"MAF missing PositionID values for input folders: {missing[:20]}")
    return ids

def main():
    parser = argparse.ArgumentParser(description="Ashlar Integrated Image & Spot Stitching Pipeline")
    parser.add_argument('-i', '--input_dir', required=True, help="Parent folder containing Position subfolders")
    parser.add_argument('-o', '--output_dir', required=True, help="Output save directory")
    parser.add_argument('-m', '--maf_file', required=True, help="MAF file path")
    parser.add_argument('-p', '--pixelsize', type=float, required=True, help="Physical pixel size (XY), um")
    parser.add_argument('-r', '--ref_channel', required=True, help="Reference registration channel (ignored if single channel)")
    parser.add_argument('-rot', '--rotation', type=int, choices=[0, 90, 180, 270], default=0, 
                        help="FOV clockwise rotation angle (0, 90, 180, 270)")
    parser.add_argument('--spot_rotation', type=int, choices=[0, 90, 180, 270], default=None,
                        help="Clockwise rotation applied to remain_reads/cell_center coordinates. Default: same as --rotation.")
    parser.add_argument('--local_qc_count', type=int, default=8,
                        help="Number of FOV-level image/point overlay QC PNGs to write. Use 0 to disable.")
    parser.add_argument('--image_dir', default=None,
                        help="Optional folder containing external FOV images, e.g. 02_registration/IF/PI/PositionXXX.tif")
    parser.add_argument('--maf_start', type=int, default=None,
                        help="First continuous MAF PositionID to stitch.")
    parser.add_argument('--maf_end', type=int, default=None,
                        help="Last continuous MAF PositionID to stitch.")
    
    args = parser.parse_args()
    
    print("\n" + "="*50)
    print("Integrated Image & Spot Stitching Pipeline Start!")
    spot_rotation = args.rotation if args.spot_rotation is None else args.spot_rotation
    print(f"(image rotate: clockwise {args.rotation} deg)")
    print(f"(spot/cell rotate: clockwise {spot_rotation} deg)")
    print("="*50)

    # Read MAF and map coordinates by explicit PositionID.
    print("\nStep 1: Read MAF and map coordinates by PositionID...")
    maf_coords_dict = get_theoretical_coords(args.maf_file)
    position_ids = choose_position_ids(args.input_dir, maf_coords_dict, args.maf_start, args.maf_end)
    channels = [""]
    actual_ref_channel = ""
    mode_text = "External single-channel PI/DAPI" if args.image_dir else "Single-channel"
    print(f"Detected mode: {mode_text}")
    local_offset_base = position_ids[0] if args.maf_start is not None else 1
    print(f"Requested MAF PositionID range/list: {position_ids[0]}-{position_ids[-1]} (n={len(position_ids)})")
    if args.maf_start is not None:
        print(f"Local Position mapping: MAF Position{position_ids[0]:03d} -> local Position001")
    
    pos_coords = []
    spot_matched_count = 0
    img_matched_count = 0
    
    maf_ids_for_output = []
    local_ids_for_output = []

    for pos_id in position_ids:
        local_pos_id = pos_id - local_offset_base + 1
        pos_dir = resolve_position_dir(args.input_dir, local_pos_id)
        image_path = resolve_image_path(args.image_dir, pos_dir, local_pos_id, args.ref_channel)
        current_coord = maf_coords_dict[pos_id]
        
        # Check Spots
        reads_file = pos_dir / 'remain_reads_raw.csv'
        cells_file = pos_dir / 'cell_center.csv'
        has_spot = reads_file.exists() and cells_file.exists()
        if has_spot:
            spot_matched_count += 1
            
        # Check Images
        has_img = image_path is not None and Path(image_path).is_file()
        if has_img:
            img_matched_count += 1
            
        if has_img:
            pos_coords.append([str(pos_dir), current_coord, str(image_path)])
            maf_ids_for_output.append(pos_id)
            local_ids_for_output.append(local_pos_id)
        else:
            print(
                f"    [WARN] skip MAF PositionID={pos_id} "
                f"(local PositionID={local_pos_id}): image missing: {image_path}"
            )

    print("    Data completeness diagnostics:")
    print(f"    - Total FOVs requested: {len(position_ids)}")
    print(f"    - FOVs with Spot data: {spot_matched_count}")
    print(f"    - FOVs with Ref Images: {img_matched_count}")
    print(f"    - Fully matched valid FOVs: {len(pos_coords)}")

    if not pos_coords:
        raise ValueError("No valid FOVs found with both Image and Spot data!")

    # Perform registration
    print(f"\nStep 2: Calculate feature drift matrix...")
    aligner = align_reference(pos_coords, actual_ref_channel, args.pixelsize, args.rotation)
    
    # Save registered FOV coordinates (Names derived purely from input_dir)
    print("\nStep 3: Save registered FOV coordinates...")
    fov_names = [Path(p[0]).name for p in pos_coords] 
    reg_coords_px = aligner.positions
    reg_coords_um = reg_coords_px * args.pixelsize
    
    df_coords = pd.DataFrame({
        'FOV': fov_names,
        'MAF_PositionID': maf_ids_for_output,
        'Local_PositionID': local_ids_for_output,
        'Global_Y_px': reg_coords_px[:, 0],
        'Global_X_px': reg_coords_px[:, 1],
        'Global_Y_um': reg_coords_um[:, 0],
        'Global_X_um': reg_coords_um[:, 1]
    })
    os.makedirs(args.output_dir, exist_ok=True)
    df_coords.to_csv(os.path.join(args.output_dir, 'registered_tilecoord.csv'), index=False)
    print("    Saved -> registered_tilecoord.csv")

    # Apply stitching for all image channels
    print("\nStep 4: Apply stitching coordinates for all image channels...")
    stitch_all_channels(aligner, pos_coords, channels, args.output_dir, args.pixelsize, args.rotation)
    
    # Spot stitching
    print("\nStep 5: Stitch spots and cell centers (Voronoi deduplication)...")
    stitch_spots(aligner.positions, aligner.metadata.size, args.output_dir, pos_coords, spot_rotation)
    write_local_overlays(
        pos_coords,
        args.output_dir,
        image_rotation=args.rotation,
        spot_rotation=spot_rotation,
        max_fovs=args.local_qc_count,
    )
    
    # Visualization
    print("\nStep 6: Generate Spot Visualizations...")
    visualize_spots(args.output_dir)
    
    print("\nAll tasks successfully completed!")

if __name__ == "__main__":
    main()
