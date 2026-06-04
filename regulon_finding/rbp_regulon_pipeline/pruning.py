from __future__ import annotations
import pandas as pd
from .io_utils import minmax_scale

def prune_edges(adj_df: pd.DataFrame, candidate_rbp_df: pd.DataFrame, rapid_df: pd.DataFrame, weight_importance=0.5, weight_rapid=0.5, min_importance=None, min_rapid_score=None, top_targets_per_rbp=100, min_targets_per_regulon=10, rapid_rbp_col='RBP', rapid_target_col='TG', rapid_score_col='score'):
    rapid_use = rapid_df[[rapid_rbp_col, rapid_target_col, rapid_score_col]].copy().rename(columns={rapid_rbp_col:'RBP', rapid_target_col:'TG', rapid_score_col:'rapid_score'})
    rbp_prior = candidate_rbp_df[['gene','candidate_prior_score','hydra_score','combined_score_01']].copy().rename(columns={'gene':'RBP'})
    merged = pd.merge(adj_df, rapid_use, on=['RBP','TG'], how='inner')
    merged = pd.merge(merged, rbp_prior, on='RBP', how='left')
    merged['importance'] = pd.to_numeric(merged['importance'], errors='coerce').fillna(0.0)
    merged['rapid_score'] = pd.to_numeric(merged['rapid_score'], errors='coerce').fillna(0.0)
    if min_importance is not None: merged = merged.loc[merged['importance'] >= float(min_importance)].copy()
    if min_rapid_score is not None: merged = merged.loc[merged['rapid_score'] >= float(min_rapid_score)].copy()
    merged['importance_norm'] = minmax_scale(merged['importance'])
    merged['rapid_score_norm'] = minmax_scale(merged['rapid_score'])
    total = weight_importance + weight_rapid
    wi, wr = weight_importance / total, weight_rapid / total
    merged['final_score'] = wi * merged['importance_norm'] + wr * merged['rapid_score_norm']
    merged = merged.sort_values(['RBP','final_score','importance','rapid_score'], ascending=[True,False,False,False]).groupby('RBP', as_index=False, group_keys=False).head(int(top_targets_per_rbp)).reset_index(drop=True)
    counts = merged.groupby('RBP')['TG'].nunique(); good_rbps = counts[counts >= int(min_targets_per_regulon)].index
    merged = merged.loc[merged['RBP'].isin(good_rbps)].copy()
    regulons = {str(rbp): dict(zip(sub['TG'].astype(str), sub['final_score'].astype(float))) for rbp, sub in merged.groupby('RBP')}
    return merged, regulons
