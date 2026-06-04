from __future__ import annotations
import pandas as pd
from .io_utils import minmax_scale

def select_candidate_rbps(hydra_df: pd.DataFrame, gene_scores_df: pd.DataFrame, hydra_gene_col='gene', hydra_score_col='seqSVM_seqDNN_ProteinBERT_score', hydra_rbp_flag_col='HydRa_RBPs', hydra_score_threshold=0.5, require_hydra_flag=True, combined_gene_col='gene', combined_score_col='combined_score_01', combined_score_threshold=0.3, selection_logic='or', expressed_genes=None) -> pd.DataFrame:
    h = hydra_df[[hydra_gene_col, hydra_score_col] + ([hydra_rbp_flag_col] if hydra_rbp_flag_col in hydra_df.columns else [])].copy().rename(columns={hydra_gene_col:'gene', hydra_score_col:'hydra_score'})
    if hydra_rbp_flag_col in h.columns:
        h['hydra_flag'] = h[hydra_rbp_flag_col].astype(bool); h = h.drop(columns=[hydra_rbp_flag_col])
    else: h['hydra_flag'] = True
    g = gene_scores_df[[combined_gene_col, combined_score_col]].copy().rename(columns={combined_gene_col:'gene', combined_score_col:'combined_score_01'})
    merged = pd.merge(h, g, on='gene', how='outer')
    merged['hydra_score'] = pd.to_numeric(merged['hydra_score'], errors='coerce').fillna(0.0)
    merged['combined_score_01'] = pd.to_numeric(merged['combined_score_01'], errors='coerce').fillna(0.0)
    merged['hydra_flag'] = merged['hydra_flag'].fillna(False)
    cond_hydra = merged['hydra_score'] >= hydra_score_threshold
    if require_hydra_flag: cond_hydra = cond_hydra & merged['hydra_flag']
    cond_combined = merged['combined_score_01'] >= combined_score_threshold
    keep = (cond_hydra & cond_combined) if selection_logic.lower() == 'and' else (cond_hydra | cond_combined)
    if expressed_genes is not None: keep = keep & merged['gene'].isin(set(expressed_genes))
    merged['is_candidate_rbp'] = keep
    merged['hydra_score_norm'] = minmax_scale(merged['hydra_score'])
    merged['combined_score_01_norm'] = minmax_scale(merged['combined_score_01'])
    merged['candidate_prior_score'] = 0.5 * merged['hydra_score_norm'] + 0.5 * merged['combined_score_01_norm']
    return merged.sort_values(['is_candidate_rbp','candidate_prior_score','hydra_score','combined_score_01'], ascending=[False,False,False,False]).reset_index(drop=True)
