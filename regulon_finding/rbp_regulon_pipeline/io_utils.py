from __future__ import annotations
import os, re
import numpy as np, pandas as pd, scanpy as sc

def sanitize_gene_names(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"^\ufeff", "", regex=True).str.strip()

def load_expression_from_h5ad(h5ad_path: str, layer: str, min_cells_expressed: int = 10):
    adata = sc.read_h5ad(h5ad_path)
    if layer not in adata.layers:
        raise KeyError(f"Layer '{layer}' not found. Available layers: {list(adata.layers.keys())}")
    X = adata.layers[layer]
    if hasattr(X, 'toarray'): X = X.toarray()
    expr = pd.DataFrame(X, index=adata.obs_names.astype(str), columns=adata.var_names.astype(str))
    keep = (expr > 0).sum(axis=0) >= int(min_cells_expressed)
    return adata, expr.loc[:, keep]

def ensure_dir(path: str): os.makedirs(path, exist_ok=True)

def minmax_scale(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors='coerce').fillna(0.0)
    vmin, vmax = float(s.min()), float(s.max())
    if vmax == vmin: return pd.Series(np.zeros(len(s)), index=s.index, dtype=float)
    return (s - vmin) / (vmax - vmin)

def read_hydra_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path); df['gene'] = sanitize_gene_names(df['gene']); return df

def read_gene_scores_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path); df['gene'] = sanitize_gene_names(df['gene']); return df

def read_rapid_table(path: str, rbp_col: str='RBP', target_col: str='TG') -> pd.DataFrame:
    try: df = pd.read_csv(path)
    except Exception: df = pd.read_csv(path, sep=None, engine='python')
    df.columns = [re.sub(r'^\ufeff', '', str(c)).strip() for c in df.columns]
    df[rbp_col] = sanitize_gene_names(df[rbp_col]); df[target_col] = sanitize_gene_names(df[target_col])
    return df
