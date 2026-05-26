import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import skimage.io

from spot_stitch import read_cells_csv, read_reads_csv, rotate_coords


def _normalize_image(img):
    img = np.asarray(img, dtype=np.float32)
    finite = img[np.isfinite(img)]
    if finite.size == 0:
        return img
    lo, hi = np.percentile(finite, [1, 99.8])
    if hi <= lo:
        return img
    return np.clip((img - lo) / (hi - lo), 0, 1)


def _read_rotated_image(path, image_rotation):
    img = skimage.io.imread(path)
    if img.ndim == 3:
        img = np.max(img, axis=0)
    return np.rot90(img, k=-(image_rotation // 90))


def _transform_local_xy(df, x_col, y_col, image_shape, spot_rotation):
    h, w = image_shape
    if spot_rotation in (90, 270):
        orig_w, orig_h = h, w
    else:
        orig_w, orig_h = w, h

    x, y = rotate_coords(
        df[x_col].to_numpy(dtype=float),
        df[y_col].to_numpy(dtype=float),
        orig_w,
        orig_h,
        spot_rotation,
    )
    return x, y


def write_local_overlays(
    pos_coords,
    output_dir,
    image_rotation,
    spot_rotation,
    max_fovs=8,
    max_points=50000,
):
    if max_fovs <= 0:
        return

    qc_dir = os.path.join(output_dir, "local_overlay_qc")
    os.makedirs(qc_dir, exist_ok=True)

    made = 0
    for entry in pos_coords[:max_fovs]:
        pos_path = entry[0]
        pos_name = os.path.basename(pos_path)
        image_path = entry[2] if len(entry) > 2 else None
        if not image_path or not os.path.isfile(image_path):
            continue

        reads = read_reads_csv(os.path.join(pos_path, "remain_reads_raw.csv"))
        cells = read_cells_csv(os.path.join(pos_path, "cell_center.csv"))
        if reads.empty and cells.empty:
            continue

        img = _read_rotated_image(image_path, image_rotation)
        img_h, img_w = img.shape

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.imshow(_normalize_image(img), cmap="gray", alpha=0.92)

        if not reads.empty:
            x, y = _transform_local_xy(
                reads,
                "spot_location_1",
                "spot_location_2",
                img.shape,
                spot_rotation,
            )
            if len(x) > max_points:
                rng = np.random.default_rng(0)
                keep = rng.choice(len(x), size=max_points, replace=False)
                x = x[keep]
                y = y[keep]
            ax.scatter(x, y, s=0.35, c="#00d5ff", alpha=0.55, linewidths=0)

        if not cells.empty:
            cx, cy = _transform_local_xy(cells, "column", "row", img.shape, spot_rotation)
            ax.scatter(cx, cy, s=10, marker="x", c="#ff3b30", alpha=0.9, linewidths=0.8)

        ax.set_xlim(0, img_w)
        ax.set_ylim(img_h, 0)
        ax.set_aspect("equal")
        ax.set_title(
            f"{pos_name}  image_rotation={image_rotation}, spot_rotation={spot_rotation}",
            fontsize=9,
        )
        ax.set_xlabel("column")
        ax.set_ylabel("row")
        plt.tight_layout()

        out_path = os.path.join(qc_dir, f"{pos_name}_overlay.png")
        plt.savefig(out_path, dpi=220)
        plt.close(fig)
        made += 1

    print(f"    Saved {made} local overlay QC image(s) -> local_overlay_qc/")
