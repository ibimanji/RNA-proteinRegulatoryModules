from __future__ import annotations
import argparse
from .config import PipelineConfig
from .pipeline import run_pipeline

def main():
    p = argparse.ArgumentParser(description='RBP regulon pipeline')
    p.add_argument('--expression_h5ad', required=True)
    p.add_argument('--expression_layer', default='totalRNA_gaussian_kernel_k3_adaptive')
    p.add_argument('--hydra_csv', required=True)
    p.add_argument('--gene_scores_csv', required=True)
    p.add_argument('--rapid_csv', required=True)
    p.add_argument('--output_dir', required=True)
    p.add_argument('--hydra_score_threshold', type=float, default=0.5)
    p.add_argument('--combined_score_threshold', type=float, default=0.3)
    p.add_argument('--selection_logic', default='or', choices=['or','and'])
    p.add_argument('--method', default='grnboost2', choices=['grnboost2','genie3'])
    p.add_argument('--top_targets_per_rbp', type=int, default=100)
    p.add_argument('--min_targets_per_regulon', type=int, default=10)
    p.add_argument('--weight_importance', type=float, default=0.5)
    p.add_argument('--weight_rapid', type=float, default=0.5)
    p.add_argument('--min_cells_expressed', type=int, default=10)
    p.add_argument('--activity_top_fraction', type=float, default=0.05)
    a = p.parse_args()
    config = PipelineConfig(expression_h5ad=a.expression_h5ad, expression_layer=a.expression_layer, hydra_csv=a.hydra_csv, gene_scores_csv=a.gene_scores_csv, rapid_csv=a.rapid_csv, output_dir=a.output_dir, hydra_score_threshold=a.hydra_score_threshold, combined_score_threshold=a.combined_score_threshold, selection_logic=a.selection_logic, arboreto_method=a.method, top_targets_per_rbp=a.top_targets_per_rbp, min_targets_per_regulon=a.min_targets_per_regulon, weight_importance=a.weight_importance, weight_rapid=a.weight_rapid, min_cells_expressed=a.min_cells_expressed, activity_top_fraction=a.activity_top_fraction)
    run_pipeline(config)

if __name__ == '__main__': main()
