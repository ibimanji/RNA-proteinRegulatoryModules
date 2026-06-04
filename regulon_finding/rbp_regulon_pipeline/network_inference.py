from __future__ import annotations
import pandas as pd

def run_arboreto_inference(expr_df: pd.DataFrame, candidate_rbps: list[str], method='grnboost2', seed=777, verbose=True) -> pd.DataFrame:
    if len(candidate_rbps) == 0: raise ValueError('No candidate RBPs available for arboreto inference.')
    from arboreto.algo import grnboost2, genie3
    if verbose:
        print(f'[network_inference] expr shape = {expr_df.shape}')
        print(f'[network_inference] candidate RBPs = {len(candidate_rbps)}')
    if method.lower() == 'grnboost2':
        adj = grnboost2(expression_data=expr_df, tf_names=candidate_rbps, seed=seed, verbose=verbose)
    elif method.lower() == 'genie3':
        adj = genie3(expression_data=expr_df, tf_names=candidate_rbps, seed=seed, verbose=verbose)
    else:
        raise ValueError("method must be 'grnboost2' or 'genie3'.")
    return adj.rename(columns={'TF':'RBP','target':'TG','importance':'importance'})
