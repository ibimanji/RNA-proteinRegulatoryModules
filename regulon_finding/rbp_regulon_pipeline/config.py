from dataclasses import dataclass
from typing import Optional

@dataclass
class PipelineConfig:
    expression_h5ad: str
    expression_layer: str = "totalRNA_gaussian_kernel_k3_adaptive"
    hydra_csv: str = ""
    gene_scores_csv: str = ""
    rapid_csv: str = ""
    output_dir: str = "rbp_regulon_output"
    hydra_score_col: str = "seqSVM_seqDNN_ProteinBERT_score"
    hydra_gene_col: str = "gene"
    hydra_rbp_flag_col: str = "HydRa_RBPs"
    hydra_score_threshold: float = 0.827
    require_hydra_flag: bool = True
    combined_gene_col: str = "gene"
    combined_score_col: str = "combined_score_01"
    combined_score_threshold: float = 0.1
    selection_logic: str = "or"
    min_cells_expressed: int = 10
    arboreto_method: str = "grnboost2"
    seed: int = 777
    verbose: bool = True
    rapid_rbp_col: str = "RBP"
    rapid_target_col: str = "TG"
    rapid_score_col: str = "score"
    min_rapid_score: Optional[float] = None
    min_importance: Optional[float] = None
    weight_importance: float = 0.5
    weight_rapid: float = 0.5
    top_targets_per_rbp: int = 100
    min_targets_per_regulon: int = 10
    activity_top_fraction: float = 0.05
    activity_use_weights: bool = True
    save_candidate_rbps: bool = True
    save_adjacencies: bool = True
    save_pruned_edges: bool = True
    arboreto_num_workers: int = 8
