#!/usr/bin/env python3

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile as tif

## /gpfs/share/home/2301920002/software/miniconda3/envs/clustermaptest/bin/python plot_clustermap_2d.py

# Edit these defaults when you want to run the script without command-line args.
DEFAULT_OUT_DIR = "/gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mWT/03_segmentation/Position078"
DEFAULT_OUTPUT = None
DEFAULT_PI_TIF = "/gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mWT/02_registration/IF/PI/Position078.tif"
DEFAULT_PI_ROTATE = 270
DEFAULT_PLOT_RAW = False
DEFAULT_NO_DAPI = False
DEFAULT_POINT_SIZE = 1.0
DEFAULT_RAW_POINT_SIZE = 0.4
DEFAULT_DPI = 300
DEFAULT_LABEL_CENTERS = False
DEFAULT_CENTER_NEIGHBORS = 8


def _hex(rgb):
    vals = np.clip(np.asarray(rgb) * 255, 0, 255).astype(int)
    return "#{:02x}{:02x}{:02x}".format(vals[0], vals[1], vals[2])


def _normalize_image(img):
    img = np.asarray(img, dtype=np.float32)
    lo, hi = np.percentile(img[np.isfinite(img)], [1, 99.8])
    if hi <= lo:
        return img
    return np.clip((img - lo) / (hi - lo), 0, 1)


def _base_palette(n=256):
    names = ["tab20", "tab20b", "tab20c", "Set1", "Set2", "Set3", "Dark2", "Paired"]
    colors = []
    for name in names:
        cmap = plt.get_cmap(name)
        count = getattr(cmap, "N", 20)
        colors.extend([cmap(i)[:3] for i in range(count)])

    # Add many HSV colors in an order that jumps around hue space.
    for i in range(max(n, 256)):
        h = (i * 0.618033988749895) % 1.0
        colors.append(plt.cm.hsv(h)[:3])

    # Drop very pale colors; white background and gray PI need saturated points.
    arr = np.asarray(colors, dtype=float)
    keep = (arr.max(axis=1) - arr.min(axis=1) > 0.25) & (arr.mean(axis=1) < 0.82)
    arr = arr[keep]

    # Deduplicate roughly.
    rounded = np.round(arr, 3)
    _, idx = np.unique(rounded, axis=0, return_index=True)
    return arr[np.sort(idx)]


def _cell_centers(reads, centers):
    if centers is not None and centers.shape[0] > 0 and {"cell_barcode", "column", "row"}.issubset(centers.columns):
        c = centers.loc[:, ["cell_barcode", "column", "row"]].copy()
        c["cell_barcode"] = pd.to_numeric(c["cell_barcode"], errors="coerce").astype("Int64")
        c = c.dropna(subset=["cell_barcode", "column", "row"])
        c["cell_barcode"] = c["cell_barcode"].astype(int)
        return c.rename(columns={"cell_barcode": "clustermap"})

    grouped = reads.groupby("clustermap", as_index=False).agg(
        column=("spot_location_1", "mean"),
        row=("spot_location_2", "mean"),
    )
    grouped["clustermap"] = pd.to_numeric(grouped["clustermap"], errors="coerce").astype(int)
    return grouped


def _assign_spatial_colors(center_df, n_neighbors):
    cell_ids = center_df["clustermap"].astype(int).to_numpy()
    xy = center_df[["column", "row"]].to_numpy(dtype=float)
    n = len(cell_ids)
    palette = _base_palette(max(n * 3, 256))

    if n == 0:
        return {}
    if n == 1:
        return {int(cell_ids[0]): palette[0]}

    dist = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(dist, np.inf)
    k = min(max(1, n_neighbors), n - 1)
    neighbors = np.argsort(dist, axis=1)[:, :k]
    degree_dist = np.take_along_axis(dist, neighbors, axis=1).sum(axis=1)
    order = np.lexsort((degree_dist, -np.isfinite(dist).sum(axis=1)))

    assigned = {}
    used = set()
    for idx in order:
        neighbor_color_idx = [assigned[j] for j in neighbors[idx] if j in assigned]
        if not neighbor_color_idx:
            pick = next((p for p in range(len(palette)) if p not in used), 0)
        else:
            nc = palette[np.asarray(neighbor_color_idx, dtype=int)]
            color_dist = np.sqrt(((palette[:, None, :] - nc[None, :, :]) ** 2).sum(axis=2))
            score = color_dist.min(axis=1)
            for p in used:
                if p < len(score):
                    score[p] -= 0.15
            pick = int(np.argmax(score))
        assigned[idx] = pick
        used.add(pick)

    return {int(cell_ids[i]): palette[p] for i, p in assigned.items()}


def parse_args():
    p = argparse.ArgumentParser(description="Plot ClusterMap 2D visualization for one Position output.")
    p.add_argument("out_dir", nargs="?", default=DEFAULT_OUT_DIR,
                   help="Output directory, e.g. .../03_segmentation/Position078")
    p.add_argument("-o", "--output", default=DEFAULT_OUTPUT,
                   help="Output PNG path. Default: OUT_DIR/clustermap_2d_visualization.png")
    p.add_argument("--pi-tif", default=DEFAULT_PI_TIF,
                   help="PI/DAPI image to overlay. Default: OUT_DIR/max_rotated_dapi.tif")
    p.add_argument("--pi-rotate", type=int, default=DEFAULT_PI_ROTATE,
                   help="Rotate PI max projection by this many degrees before overlay. Use 0 for already-rotated images.")
    p.add_argument("--raw", action="store_true", default=DEFAULT_PLOT_RAW,
                   help="Also plot unassigned reads from remain_reads_raw.csv in gray.")
    p.add_argument("--no-dapi", action="store_true", default=DEFAULT_NO_DAPI,
                   help="Do not plot max_rotated_dapi.tif as background.")
    p.add_argument("--point-size", type=float, default=DEFAULT_POINT_SIZE, help="Assigned read point size.")
    p.add_argument("--raw-point-size", type=float, default=DEFAULT_RAW_POINT_SIZE, help="Unassigned read point size.")
    p.add_argument("--dpi", type=int, default=DEFAULT_DPI, help="Output image DPI.")
    p.add_argument("--label-centers", action="store_true", default=DEFAULT_LABEL_CENTERS,
                   help="Write cell IDs next to center markers.")
    p.add_argument("--center-neighbors", type=int, default=DEFAULT_CENTER_NEIGHBORS,
                   help="Nearby cells considered when choosing contrasting colors.")
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = os.path.abspath(args.out_dir)
    output = args.output or os.path.join(out_dir, "clustermap_2d_visualization.png")

    reads_csv = os.path.join(out_dir, "remain_reads.csv")
    raw_csv = os.path.join(out_dir, "remain_reads_raw.csv")
    center_csv = os.path.join(out_dir, "cell_center.csv")
    dapi_tif = args.pi_tif or os.path.join(out_dir, "max_rotated_dapi.tif")

    if not os.path.isfile(reads_csv):
        raise FileNotFoundError(f"Missing remain_reads.csv: {reads_csv}")

    reads = pd.read_csv(reads_csv)
    if reads.shape[0] == 0:
        raise ValueError(f"remain_reads.csv is empty: {reads_csv}")

    required = {"spot_location_1", "spot_location_2", "clustermap"}
    missing = sorted(required - set(reads.columns))
    if missing:
        raise ValueError(f"remain_reads.csv missing columns: {missing}")

    reads = reads.copy()
    reads["clustermap"] = pd.to_numeric(reads["clustermap"], errors="coerce").fillna(-1).astype(int)
    reads = reads.loc[reads["clustermap"] >= 0].copy()

    centers = pd.read_csv(center_csv) if os.path.isfile(center_csv) else None
    center_df = _cell_centers(reads, centers)
    color_map = _assign_spatial_colors(center_df, args.center_neighbors)
    fallback = np.array([1.0, 0.0, 0.0])
    point_colors = np.vstack([color_map.get(int(cid), fallback) for cid in reads["clustermap"]])

    key = pd.DataFrame({
        "cell_barcode": sorted(color_map),
        "color_hex": [_hex(color_map[cid]) for cid in sorted(color_map)],
    })
    key_path = os.path.join(out_dir, "cell_color_key.csv")
    key.to_csv(key_path, index=False)

    fig, ax = plt.subplots(figsize=(12, 12))

    if not args.no_dapi and os.path.isfile(dapi_tif):
        dapi = tif.imread(dapi_tif)
        if dapi.ndim == 3:
            dapi = dapi.max(axis=0)
        if args.pi_rotate % 360 != 0:
            dapi = np.rot90(dapi, k=(args.pi_rotate // 90) % 4)
        ax.imshow(_normalize_image(dapi), cmap="gray", alpha=0.9)

    if args.raw and os.path.isfile(raw_csv):
        raw = pd.read_csv(raw_csv)
        if {"spot_location_1", "spot_location_2", "clustermap"}.issubset(raw.columns):
            raw_cm = pd.to_numeric(raw["clustermap"], errors="coerce").fillna(-1)
            unassigned = raw.loc[raw_cm < 0]
            if unassigned.shape[0] > 0:
                ax.scatter(
                    unassigned["spot_location_1"],
                    unassigned["spot_location_2"],
                    s=args.raw_point_size,
                    c="lightgray",
                    alpha=0.25,
                    linewidths=0,
                    label="unassigned",
                )

    ax.scatter(
        reads["spot_location_1"],
        reads["spot_location_2"],
        c=point_colors,
        s=args.point_size,
        alpha=0.82,
        linewidths=0,
        label="assigned",
    )

    if center_df.shape[0] > 0:
        center_colors = np.vstack([color_map.get(int(cid), fallback) for cid in center_df["clustermap"]])
        ax.scatter(
            center_df["column"],
            center_df["row"],
            s=62,
            c="white",
            marker="x",
            linewidths=2.4,
            label="cell center",
            zorder=5,
        )
        ax.scatter(
            center_df["column"],
            center_df["row"],
            s=42,
            c=center_colors,
            marker="x",
            linewidths=1.6,
            zorder=6,
        )
        if args.label_centers:
            for _, row in center_df.iterrows():
                ax.text(
                    row["column"] + 8,
                    row["row"] + 8,
                    str(int(row["clustermap"])),
                    fontsize=5,
                    color="white",
                    path_effects=[],
                    zorder=7,
                )

    title = os.path.basename(os.path.normpath(out_dir))
    ax.set_title(f"{title} reads colored by stitched cell ID")
    ax.set_xlabel("column / spot_location_1")
    ax.set_ylabel("row / spot_location_2")
    ax.set_aspect("equal")
    ax.legend(loc="upper right", markerscale=4)

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    plt.tight_layout()
    plt.savefig(output, dpi=args.dpi)
    print(f"saved: {output}")
    print(f"saved: {key_path}")


if __name__ == "__main__":
    main()
