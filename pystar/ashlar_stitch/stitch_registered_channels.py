#!/usr/bin/env python3
import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile as tif


def position_id_from_name(name):
    match = re.search(r"\d+", str(name))
    if not match:
        raise ValueError(f"Cannot parse Position ID from {name!r}")
    return int(match.group())


def position_name(pos_id):
    return f"Position{int(pos_id):03d}"


def resolve_image_path(channel_dir, fov_name):
    root = Path(channel_dir).resolve()
    pos_id = position_id_from_name(fov_name)
    candidates = [
        root / f"{position_name(pos_id)}.tif",
        root / f"Position{pos_id}.tif",
    ]
    for path in candidates:
        if path.is_file():
            return path
    return candidates[0]


def read_rotated_projection(path, rotation):
    img = tif.imread(path)
    if img.ndim == 3:
        img = np.max(img, axis=0)
    elif img.ndim > 3:
        raise ValueError(f"Expected 2D or 3D TIFF, got shape={img.shape}: {path}")
    return np.rot90(img, k=-(rotation // 90))


def stitch_channel(coords, channel_dir, output_path, rotation, blend):
    first_img = None
    matched = []

    for _, row in coords.iterrows():
        image_path = resolve_image_path(channel_dir, row["FOV"])
        if not image_path.is_file():
            print(f"[WARN] missing image for {row['FOV']}: {image_path}", flush=True)
            continue
        if first_img is None:
            first_img = read_rotated_projection(image_path, rotation)
        matched.append((row, image_path))

    if first_img is None or not matched:
        raise ValueError(f"No images found under {channel_dir}")

    tile_h, tile_w = first_img.shape
    xs = np.rint([row["Global_X_px"] for row, _ in matched]).astype(int)
    ys = np.rint([row["Global_Y_px"] for row, _ in matched]).astype(int)
    shift_x = max(0, -int(xs.min()))
    shift_y = max(0, -int(ys.min()))
    canvas_w = int(xs.max()) + shift_x + tile_w
    canvas_h = int(ys.max()) + shift_y + tile_h

    canvas = np.zeros((canvas_h, canvas_w), dtype=first_img.dtype)
    count = None
    if blend == "mean":
        canvas = np.zeros((canvas_h, canvas_w), dtype=np.float64)
        count = np.zeros((canvas_h, canvas_w), dtype=np.uint16)

    for idx, (row, image_path) in enumerate(matched, start=1):
        img = first_img if idx == 1 else read_rotated_projection(image_path, rotation)
        if img.shape != (tile_h, tile_w):
            raise ValueError(
                f"Tile shape mismatch for {image_path}: got {img.shape}, expected {(tile_h, tile_w)}"
            )

        x0 = int(round(row["Global_X_px"])) + shift_x
        y0 = int(round(row["Global_Y_px"])) + shift_y
        y1 = y0 + tile_h
        x1 = x0 + tile_w

        if blend == "max":
            canvas[y0:y1, x0:x1] = np.maximum(canvas[y0:y1, x0:x1], img)
        elif blend == "overwrite":
            canvas[y0:y1, x0:x1] = img
        elif blend == "mean":
            canvas[y0:y1, x0:x1] += img.astype(np.float64, copy=False)
            count[y0:y1, x0:x1] += 1
        else:
            raise ValueError(f"Unknown blend mode: {blend}")

        if idx % 10 == 0 or idx == len(matched):
            print(f"    stitched {idx}/{len(matched)} tiles", flush=True)

    if blend == "mean":
        nonzero = count > 0
        canvas[nonzero] = canvas[nonzero] / count[nonzero]
        canvas = np.clip(canvas, 0, np.iinfo(first_img.dtype).max).astype(first_img.dtype)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    tif.imwrite(output_path, canvas, photometric="minisblack", bigtiff=True)
    print(f"[OK] saved {output_path} shape={canvas.shape} dtype={canvas.dtype}", flush=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Stitch external IF channels with an existing registered_tilecoord.csv."
    )
    parser.add_argument("--registered_tilecoord", required=True,
                        help="registered_tilecoord.csv from the verified ashlar stitch run.")
    parser.add_argument("--if_dir", required=True,
                        help="IF folder containing channel subfolders, for example .../02_registration/IF")
    parser.add_argument("--channels", nargs="+", required=True,
                        help="Channel folder names under --if_dir, for example PI X34")
    parser.add_argument("--output_dir", required=True,
                        help="Output folder for stitched TIFFs.")
    parser.add_argument("--rotation", type=int, choices=[0, 90, 180, 270], default=0,
                        help="Clockwise image rotation before stitching. Use the same image rotation as the verified run.")
    parser.add_argument("--blend", choices=["max", "mean", "overwrite"], default="max",
                        help="How to combine overlaps. max is conservative for fluorescence max projections.")
    return parser.parse_args()


def main():
    args = parse_args()
    coords = pd.read_csv(args.registered_tilecoord)
    required = {"FOV", "Global_Y_px", "Global_X_px"}
    missing = sorted(required - set(coords.columns))
    if missing:
        raise ValueError(f"registered_tilecoord missing columns: {missing}")

    coords = coords.sort_values(
        by="FOV",
        key=lambda s: s.map(position_id_from_name),
    ).reset_index(drop=True)

    for channel in args.channels:
        channel_dir = os.path.join(args.if_dir, channel)
        out_path = os.path.join(args.output_dir, f"{channel}_registered_stitched.tif")
        print(f"\n[{channel}] channel_dir={channel_dir}", flush=True)
        stitch_channel(coords, channel_dir, out_path, args.rotation, args.blend)


if __name__ == "__main__":
    main()
