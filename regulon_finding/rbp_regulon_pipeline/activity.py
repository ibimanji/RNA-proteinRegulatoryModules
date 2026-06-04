from __future__ import annotations
import numpy as np, pandas as pd

def _weighted_auc_for_one_cell(ranked_genes: pd.Index, regulon: dict[str,float], top_k: int, use_weights=True) -> float:
    top_genes = ranked_genes[:top_k]; hits = []; weights = []
    for i, g in enumerate(top_genes, start=1):
        if g in regulon: hits.append(i); weights.append(float(regulon[g]) if use_weights else 1.0)
    if not hits: return 0.0
    hit_mask = np.zeros(top_k, dtype=float)
    for pos, w in zip(hits, weights): hit_mask[pos-1] = w
    cum = np.cumsum(hit_mask); denom = hit_mask.sum() * top_k
    return 0.0 if denom <= 0 else float(cum.sum() / denom)

def score_regulon_activity(expr_df: pd.DataFrame, regulons: dict[str,dict[str,float]], top_fraction=0.05, use_weights=True) -> pd.DataFrame:
    if not regulons: return pd.DataFrame(index=expr_df.index)
    top_k = max(1, int(round(expr_df.shape[1] * top_fraction)))
    ranked_gene_lists = {cell: expr_df.loc[cell].sort_values(ascending=False).index for cell in expr_df.index}
    activity = pd.DataFrame(index=expr_df.index, columns=list(regulons.keys()), dtype=float)
    for cell in expr_df.index:
        ranked_genes = ranked_gene_lists[cell]
        for rbp, gene2weight in regulons.items():
            activity.at[cell, rbp] = _weighted_auc_for_one_cell(ranked_genes, gene2weight, top_k, use_weights)
    return activity
