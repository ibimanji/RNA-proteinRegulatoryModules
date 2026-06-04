# RBP regulon pipeline

Workflow:

```text
expression matrix
   ↓
select candidate RBPs by HydRa / combined_score_01
   ↓
run arboreto using these candidate RBPs
   ↓
prune edges by RAPID + edge importance score
   ↓
build regulons
   ↓
score regulon activity
```

Outputs:
- candidate_rbps.csv
- arboreto_adjacencies.csv
- pruned_edges.csv
- regulons.json
- regulon_summary.csv
- regulon_activity.csv

Run:
```bash
python -m rbp_regulon_pipeline.cli   --expression_h5ad /media/zenglab/result/lingyuan/STEM/aggregate/GaussianKernel_results/k_3_mult1.5/processed_data_k3_adaptive.h5ad   --expression_layer totalRNA_gaussian_kernel_k3_adaptive   --hydra_csv /media/zenglab/result/lingyuan/STEM/mechanistic_plausibility/hydra/mouse_res/predict_mouse_2322genes_HydRa_predictions_with_gene.csv   --gene_scores_csv /media/zenglab/result/lingyuan/STEM/mechanistic_plausibility/gene_integrated_scores_weighted.csv   --rapid_csv /media/zenglab/result/lingyuan/STEM/mechanistic_plausibility/scRAPID/catrapid_RBP_TG_both_in_detected_genes.csv   --output_dir rbp_regulon_output
```
