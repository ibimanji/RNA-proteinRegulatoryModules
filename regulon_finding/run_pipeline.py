from rbp_regulon_pipeline import PipelineConfig, run_pipeline


def main():
    config = PipelineConfig(
        expression_h5ad='/media/zenglab/result/lingyuan/STEM/aggregate/GaussianKernel_results/k_3_mult1.5/processed_data_k3_adaptive.h5ad',
        expression_layer='totalRNA_gaussian_kernel_k3_adaptive',
        hydra_csv='/media/zenglab/result/lingyuan/STEM/mechanistic_plausibility/hydra/mouse_res/predict_mouse_2322genes_HydRa_predictions_with_gene.csv',
        gene_scores_csv='/media/zenglab/result/lingyuan/STEM/mechanistic_plausibility/gene_integrated_scores_weighted.csv',
        rapid_csv='/media/zenglab/result/lingyuan/STEM/mechanistic_plausibility/scRAPID/catrapid_RBP_TG_both_in_detected_genes.csv',
        output_dir='rbp_regulon_output',
        hydra_score_threshold=0.827,
        combined_score_threshold=0.2,
        selection_logic='or',
        arboreto_method='grnboost2',
        weight_importance=0.5,
        weight_rapid=0.5,
        top_targets_per_rbp=100,
        min_targets_per_regulon=10,
        activity_top_fraction=0.05,
        activity_use_weights=True,
        arboreto_num_workers=8,
    )

    results = run_pipeline(config)
    print(results['summary_df'].head())
    print(results['activity_df'].shape)


if __name__ == "__main__":
    main()