from __future__ import annotations
import os, json, pandas as pd
from .config import PipelineConfig
from .io_utils import ensure_dir, load_expression_from_h5ad, read_hydra_table, read_gene_scores_table, read_rapid_table
from .candidate_rbps import select_candidate_rbps
from .network_inference import run_arboreto_inference
from .pruning import prune_edges
from .activity import score_regulon_activity

def run_pipeline(config: PipelineConfig):
    ensure_dir(config.output_dir)
    _, expr_df = load_expression_from_h5ad(config.expression_h5ad, layer=config.expression_layer, min_cells_expressed=config.min_cells_expressed)
    hydra_df = read_hydra_table(config.hydra_csv)
    gene_scores_df = read_gene_scores_table(config.gene_scores_csv)
    rapid_df = read_rapid_table(config.rapid_csv, rbp_col=config.rapid_rbp_col, target_col=config.rapid_target_col)
    candidate_rbp_df = select_candidate_rbps(hydra_df, gene_scores_df, hydra_gene_col=config.hydra_gene_col, hydra_score_col=config.hydra_score_col, hydra_rbp_flag_col=config.hydra_rbp_flag_col, hydra_score_threshold=config.hydra_score_threshold, require_hydra_flag=config.require_hydra_flag, combined_gene_col=config.combined_gene_col, combined_score_col=config.combined_score_col, combined_score_threshold=config.combined_score_threshold, selection_logic=config.selection_logic, expressed_genes=expr_df.columns.tolist())
    candidate_rbps = candidate_rbp_df.loc[candidate_rbp_df['is_candidate_rbp'], 'gene'].tolist()

    print("[debug] expr_df.shape:", expr_df.shape)
    print("[debug] expr_df first 20 columns:", expr_df.columns[:20].tolist())

    print("[debug] candidate_rbps raw len:", len(candidate_rbps))
    print("[debug] candidate_rbps first 20:", list(candidate_rbps[:20]))

    expr_genes = pd.Index(expr_df.columns).astype(str)
    cand_genes = pd.Index(candidate_rbps).astype(str)

    print("[debug] unique candidate_rbps:", cand_genes.nunique())
    print("[debug] overlap len:", len(expr_genes.intersection(cand_genes)))
    print("[debug] first 20 overlap:", expr_genes.intersection(cand_genes)[:20].tolist())
    print("[debug] first 20 candidate not in expr:", cand_genes.difference(expr_genes)[:20].tolist())

    adj_df = run_arboreto_inference(expr_df, candidate_rbps, method=config.arboreto_method, seed=config.seed, verbose=config.verbose,num_workers=config.arboreto_num_workers,)
    pruned_edges_df, regulons = prune_edges(adj_df, candidate_rbp_df, rapid_df, weight_importance=config.weight_importance, weight_rapid=config.weight_rapid, min_importance=config.min_importance, min_rapid_score=config.min_rapid_score, top_targets_per_rbp=config.top_targets_per_rbp, min_targets_per_regulon=config.min_targets_per_regulon, rapid_rbp_col=config.rapid_rbp_col, rapid_target_col=config.rapid_target_col, rapid_score_col=config.rapid_score_col)
    activity_df = score_regulon_activity(expr_df, regulons, top_fraction=config.activity_top_fraction, use_weights=config.activity_use_weights)
    if config.save_candidate_rbps: candidate_rbp_df.to_csv(os.path.join(config.output_dir, 'candidate_rbps.csv'), index=False)
    if config.save_adjacencies: adj_df.to_csv(os.path.join(config.output_dir, 'arboreto_adjacencies.csv'), index=False)
    if config.save_pruned_edges: pruned_edges_df.to_csv(os.path.join(config.output_dir, 'pruned_edges.csv'), index=False)
    activity_df.to_csv(os.path.join(config.output_dir, 'regulon_activity.csv'))
    summary_df = pd.DataFrame({'RBP': list(regulons.keys()), 'n_targets': [len(v) for v in regulons.values()]}).sort_values('n_targets', ascending=False)
    summary_df.to_csv(os.path.join(config.output_dir, 'regulon_summary.csv'), index=False)
    with open(os.path.join(config.output_dir, 'regulons.json'), 'w') as f: json.dump(regulons, f, indent=2)
    return {'candidate_rbp_df': candidate_rbp_df, 'adj_df': adj_df, 'pruned_edges_df': pruned_edges_df, 'regulons': regulons, 'activity_df': activity_df, 'summary_df': summary_df}
