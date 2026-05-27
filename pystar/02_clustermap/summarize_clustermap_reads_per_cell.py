#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_SEGMENTATION_DIR = Path(
    "/gpfs/share/home/2301920002/labShare/2301920002/ADdecon/14mWT_14mFAD/14mWT/03_segmentation"
)
DEFAULT_OUTPUT_DIR = None


def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize reads-per-cell distribution from ClusterMap remain_reads.csv files."
    )
    p.add_argument("--segmentation-dir", type=Path, default=DEFAULT_SEGMENTATION_DIR,
                   help="Directory containing Position*/remain_reads.csv.")
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
                   help="Output directory. Default: SEGMENTATION_DIR/clustermap_reads_summary.")
    p.add_argument("--min-reads", type=int, default=1,
                   help="Only include cells with at least this many assigned reads.")
    p.add_argument("--log-x", action="store_true",
                   help="Use log10 x-axis for histogram.")
    p.add_argument("--bins", type=int, default=80,
                   help="Number of histogram bins.")
    return p.parse_args()


def position_sort_key(path):
    name = path.name
    digits = "".join(ch for ch in name if ch.isdigit())
    return int(digits) if digits else name


def _rna_mask(df, rna_type):
    if rna_type == "totalRNA":
        return pd.Series(True, index=df.index)

    text_cols = [c for c in ("gene_name", "gene") if c in df.columns]
    if not text_cols:
        return pd.Series(False, index=df.index)

    text = pd.Series("", index=df.index, dtype="object")
    for col in text_cols:
        text = text.str.cat(df[col].astype(str), sep=" ")
    return text.str.contains(rna_type, case=False, na=False, regex=False)


def load_counts(segmentation_dir, min_reads):
    rows_by_type = {rna_type: [] for rna_type in ("totalRNA", "ntRNA", "rbRNA")}
    missing = []
    for pos_dir in sorted(segmentation_dir.glob("Position*"), key=position_sort_key):
        if not pos_dir.is_dir():
            continue
        reads_csv = pos_dir / "remain_reads.csv"
        if not reads_csv.exists():
            missing.append(pos_dir.name)
            continue

        df = pd.read_csv(reads_csv, usecols=lambda c: c in {"clustermap", "gene_name", "gene"})
        if "clustermap" not in df.columns or df.empty:
            continue

        df = df.copy()
        df["clustermap"] = pd.to_numeric(df["clustermap"], errors="coerce")
        df = df.dropna(subset=["clustermap"])
        df["clustermap"] = df["clustermap"].astype(int)
        df = df.loc[df["clustermap"] >= 0]

        for rna_type in rows_by_type:
            sub = df.loc[_rna_mask(df, rna_type)]
            counts = sub["clustermap"].value_counts().sort_index()
            counts = counts.loc[counts >= min_reads]
            for cell_id, n_reads in counts.items():
                rows_by_type[rna_type].append({
                    "rna_type": rna_type,
                    "position": pos_dir.name,
                    "cell_barcode": int(cell_id),
                    "n_reads": int(n_reads),
                })

    return {k: pd.DataFrame(v) for k, v in rows_by_type.items()}, missing


def summarize_counts(cell_counts):
    if cell_counts.empty:
        return pd.DataFrame([{
            "rna_type": "",
            "n_positions": 0,
            "n_cells": 0,
            "total_reads": 0,
            "mean_reads_per_cell": np.nan,
            "median_reads_per_cell": np.nan,
            "std_reads_per_cell": np.nan,
            "min_reads_per_cell": np.nan,
            "max_reads_per_cell": np.nan,
        }])

    rna_type = str(cell_counts["rna_type"].iloc[0]) if "rna_type" in cell_counts.columns else ""
    return pd.DataFrame([{
        "rna_type": rna_type,
        "n_positions": int(cell_counts["position"].nunique()),
        "n_cells": int(cell_counts.shape[0]),
        "total_reads": int(cell_counts["n_reads"].sum()),
        "mean_reads_per_cell": float(cell_counts["n_reads"].mean()),
        "median_reads_per_cell": float(cell_counts["n_reads"].median()),
        "std_reads_per_cell": float(cell_counts["n_reads"].std(ddof=1)),
        "min_reads_per_cell": int(cell_counts["n_reads"].min()),
        "max_reads_per_cell": int(cell_counts["n_reads"].max()),
    }])


def write_one_rna_outputs(rna_type, cell_counts, output_dir, bins, log_x):
    output_dir.mkdir(parents=True, exist_ok=True)

    cell_csv = output_dir / f"{rna_type}_reads_per_cell.csv"
    position_csv = output_dir / f"{rna_type}_reads_per_cell_by_position_summary.csv"
    overall_csv = output_dir / f"{rna_type}_reads_per_cell_overall_summary.csv"
    png = output_dir / f"{rna_type}_reads_per_cell_distribution.png"

    cell_counts.to_csv(cell_csv, index=False)

    if cell_counts.empty:
        overall = summarize_counts(cell_counts)
        overall.loc[0, "rna_type"] = rna_type
        overall.to_csv(overall_csv, index=False)
        return cell_csv, position_csv, overall_csv, png, overall

    by_pos = (
        cell_counts
        .groupby("position")["n_reads"]
        .agg(
            n_cells="size",
            mean_reads_per_cell="mean",
            median_reads_per_cell="median",
            min_reads_per_cell="min",
            max_reads_per_cell="max",
        )
        .reset_index()
    )
    by_pos.to_csv(position_csv, index=False)

    overall = summarize_counts(cell_counts)
    overall.to_csv(overall_csv, index=False)

    values = cell_counts["n_reads"].to_numpy()
    fig, ax = plt.subplots(figsize=(9, 6))
    if log_x:
        positive = values[values > 0]
        edges = np.logspace(np.log10(positive.min()), np.log10(positive.max()), bins + 1)
        ax.hist(positive, bins=edges, color="#4c78a8", alpha=0.85)
        ax.set_xscale("log")
    else:
        ax.hist(values, bins=bins, color="#4c78a8", alpha=0.85)

    mean_v = overall.loc[0, "mean_reads_per_cell"]
    median_v = overall.loc[0, "median_reads_per_cell"]
    ax.axvline(mean_v, color="#e45756", linewidth=2, label=f"mean = {mean_v:.2f}")
    ax.axvline(median_v, color="#54a24b", linewidth=2, label=f"median = {median_v:.2f}")
    ax.set_title(f"ClusterMap {rna_type} reads per cell distribution")
    ax.set_xlabel("reads per cell")
    ax.set_ylabel("number of cells")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(png, dpi=300)
    plt.close(fig)

    return cell_csv, position_csv, overall_csv, png, overall


def write_outputs(counts_by_type, missing, output_dir, bins, log_x):
    output_dir.mkdir(parents=True, exist_ok=True)
    missing_txt = output_dir / "positions_missing_remain_reads.txt"
    missing_txt.write_text("\n".join(missing) + ("\n" if missing else ""))

    outputs = {}
    overall_rows = []
    for rna_type in ("totalRNA", "ntRNA", "rbRNA"):
        cell_csv, position_csv, overall_csv, png, overall = write_one_rna_outputs(
            rna_type, counts_by_type[rna_type], output_dir, bins, log_x
        )
        outputs[rna_type] = {
            "cell_csv": cell_csv,
            "position_csv": position_csv,
            "overall_csv": overall_csv,
            "png": png,
            "overall": overall,
        }
        overall_rows.append(overall)

    combined_overall = pd.concat(overall_rows, ignore_index=True)
    combined_csv = output_dir / "totalRNA_ntRNA_rbRNA_reads_per_cell_overall_summary.csv"
    combined_overall.to_csv(combined_csv, index=False)
    return outputs, missing_txt, combined_csv, combined_overall


def main():
    args = parse_args()
    segmentation_dir = args.segmentation_dir
    output_dir = args.output_dir or (segmentation_dir / "clustermap_reads_summary")

    counts_by_type, missing = load_counts(segmentation_dir, args.min_reads)
    outputs, missing_txt, combined_csv, combined_overall = write_outputs(
        counts_by_type, missing, output_dir, args.bins, args.log_x
    )

    for _, row in combined_overall.iterrows():
        print(
            f"[DONE] {row['rna_type']}: cells={int(row['n_cells'])} "
            f"positions={int(row['n_positions'])} "
            f"mean={row['mean_reads_per_cell']:.3f} "
            f"median={row['median_reads_per_cell']:.3f}"
        )
    for rna_type, paths in outputs.items():
        print(f"[DONE] {rna_type} cell counts: {paths['cell_csv']}")
        print(f"[DONE] {rna_type} position summary: {paths['position_csv']}")
        print(f"[DONE] {rna_type} overall summary: {paths['overall_csv']}")
        print(f"[DONE] {rna_type} plot: {paths['png']}")
    print(f"[DONE] combined overall summary: {combined_csv}")
    print(f"[DONE] missing list: {missing_txt}")


if __name__ == "__main__":
    main()
