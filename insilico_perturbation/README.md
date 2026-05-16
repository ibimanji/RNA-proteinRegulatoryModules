# Virtual RBP Knockout Simulation Pipeline

This repository contains two complementary computational pipelines for performing **in silico (virtual) RNA-Binding Protein (RBP) knockout experiments** on single-cell RNA-seq data. Both pipelines operate on the same input gene regulatory network but employ fundamentally different modeling approaches, making their results mutually validating.

---

## Repository Structure

```
.
├── Sctenifoldpy.ipynb       # Pipeline 1: scTenifold manifold alignment-based KO
└── CellOracle.ipynb      # Pipeline 2: CellOracle GRN simulation-based KO
```

---

## Pipeline 1 — scTenifold-based Virtual Knockout (`Sctenifoldpy.ipynb`)

### Overview

This pipeline uses the **scTenifold** framework to perform virtual RBP knockouts via **manifold alignment**. For each RBP, the pipeline constructs a "knockout" gene regulatory network by zeroing out all outgoing edges from that RBP, then compares the wild-type (WT) and knockout (KO) networks in a shared low-dimensional manifold to identify significantly perturbed downstream target genes.

### Input

| File | Description |
|---|---|
| `pruned_edges.csv` | Pre-computed RBP-Target regulatory network with columns: `RBP` (source RNA-binding protein), `TG` (target gene), `final_score` (edge weight representing regulatory strength) |

### Workflow

1. **Environment Setup** — Imports core libraries (`scTenifold`, `pandas`, `numpy`, `networkx`) and sets global parameters including manifold alignment dimensions (`MANIFOLD_DIM=15`) and random seed for reproducibility.

2. **WT Network Construction** — Reads `pruned_edges.csv` and builds a global wild-type (WT) weighted adjacency matrix. All unique RBPs and target genes form the node set; `final_score` values populate the directed edges. Self-regulatory loops (diagonal) are set to zero.

3. **Batch Virtual Knockout Loop** — Iterates over every RBP in the network:
   - Creates a KO adjacency matrix by zeroing all outgoing edges from the target RBP.
   - Runs `manifold_alignment()` to co-embed WT and KO networks into a shared `MANIFOLD_DIM`-dimensional space.
   - Runs `d_regulation()` (differential regulation test) to compute perturbation distance, fold change, p-value, and FDR-adjusted p-value for each downstream target gene.
   - Flags target genes as significant at FDR < 0.05.

4. **Results Aggregation** — Concatenates all per-RBP results into a master dataframe, sorted by RBP and adjusted p-value.

5. **Visualization** — For every RBP with at least one significant target gene, generates two plots:
   - **Ego-Network** (`EgoNetwork_{RBP}.png`): A directed network graph where node size encodes perturbation distance, edge width encodes original regulatory weight, and node color indicates significance (red = significant, grey = not significant).
   - **Lollipop Chart** (`Lollipop_{RBP}.png`): A horizontal ranked chart of the top 30 most-affected target genes by perturbation distance.

### Output

```
Sctenifoldpy_RBP_KO_Analysis_Results/
├── All_Tested_RBP_Targets_Report.csv    # Full results for all RBP-Target pairs tested
├── Significant_Only_Targets_Report.csv  # Filtered results: only FDR < 0.05 pairs
└── KO_{RBP}/                            # One folder per RBP with significant targets
    ├── EgoNetwork_{RBP}.png
    └── Lollipop_{RBP}.png
```

**Key columns in output CSVs:**

| Column | Description |
|---|---|
| `KO_RBP` | The knocked-out RNA-binding protein |
| `Target_Gene` | The downstream target gene being evaluated |
| `original_edge_weight` | Regulatory strength from the input network |
| `Distance` | Perturbation distance in the aligned manifold (larger = more affected) |
| `FC` | Fold change after KO |
| `p-value` | Raw p-value from differential regulation test |
| `adjusted p-value` | FDR-corrected p-value (Benjamini-Hochberg) |
| `is_significant` | Boolean flag: `True` if adjusted p-value < 0.05 |

---

## Pipeline 2 — CellOracle-based Virtual Knockout (`CellOracle.ipynb`)

### Overview

This pipeline uses **CellOracle** to simulate RBP knockouts within the context of single-cell transcriptomic data. Unlike scTenifold which operates purely on the network topology, CellOracle integrates the regulatory network with actual single-cell gene expression, fits a linear regression GRN model per cell cluster, and then propagates the knockout perturbation signal through the network to predict cell state transitions and downstream expression changes.

### Input

| File | Description |
|---|---|
| `processed_data_k3_adaptive.h5ad` | Processed single-cell RNA-seq AnnData object. Must contain: raw count matrix, `X_pca_harmony` (batch-corrected PCA embeddings), and cluster labels in `obs["rna_nn_alg1_label3"]` |
| `pruned_edges.csv` | Same RBP-Target regulatory network as Pipeline 1 (columns: `RBP`, `TG`, `final_score`) |

### Workflow

1. **Environment Setup & Data Loading** — Loads the `.h5ad` single-cell dataset and initializes a `co.Oracle` object. Recomputes the cell neighborhood graph using `X_pca_harmony` (Harmony batch-corrected PCA) and recalculates UMAP coordinates.

2. **KNN Imputation** — Performs balanced K-Nearest Neighbor imputation (`oracle.knn_imputation`) on the raw count matrix to smooth out technical noise. The value of `k` is dynamically set to 2.5% of total cell count (clamped between 20 and 100).

3. **GRN Integration** — Loads `pruned_edges.csv`, renames columns to CellOracle's expected format (`source`, `target`, `coef_mean`), and assigns the same network to all cell clusters. Builds a unified TF→Target dictionary and fits a cluster-specific linear regression GRN model (`oracle.fit_GRN_for_simulation`, `alpha=10`).

4. **Batch KO — Cell State Transition (Vector Field)** — For each valid RBP (present in both the network and the single-cell expression matrix)(No need, Optional):
   - Sets the RBP's expression to zero and propagates the signal for 3 steps (`oracle.simulate_shift`).
   - Estimates transition probabilities to neighboring cell states in high-dimensional space (`oracle.estimate_transition_prob`).
   - Projects the transition vectors onto the 2D UMAP embedding (`oracle.calculate_embedding_shift`).
   - Saves a dual-panel figure: left panel shows the KO-specific cell flow field; right panel shows a randomized GRN control to confirm biological specificity.

5. **Batch KO — Downstream Target Expression** — Re-runs `simulate_shift` for each RBP to extract the top 8 most strongly regulated downstream targets (by absolute `final_score`). For each target, compares `imputed_count` (WT) vs. `simulated_count` (KO) across all cells and saves:
   - A **bar plot** showing mean WT vs. KO expression per target gene.
   - A **split violin plot** showing the full expression distribution.

6. **Summary Table Generation** — Computes quantitative statistics for every RBP–Target pair: mean WT expression, mean KO expression, absolute delta, and percent change. Saves all records to a single CSV.

### Output

```
Oracle_KO_Results/
└── KO_VectorField_{RBP}.png          # UMAP vector field: cell state transitions after KO (Optional)

cellorcle_KO_Expression_Plots/
└── Expression_KO_{RBP}.png           # Bar + violin plots: target gene expression WT vs KO

KO_summary_expression_change.csv      # Quantitative summary table of all KO effects
```

**Key columns in `KO_summary_expression_change.csv`:**

| Column | Description |
|---|---|
| `KO_RBP` | The knocked-out RNA-binding protein |
| `Target_Gene` | Downstream target gene evaluated |
| `Mean_WT` | Mean imputed expression before knockout |
| `Mean_KO` | Mean simulated expression after knockout |
| `Delta` | Absolute expression change (Mean_KO − Mean_WT) |
| `Delta_pct` | Percent change relative to WT expression |
| `Abs_Delta` | Absolute value of Delta (used for ranking) |

---

## Dependencies

### Pipeline 1 (scTenifold)
```
scTenifold
pandas
numpy
matplotlib
networkx
seaborn
tqdm
```

### Pipeline 2 (CellOracle)
```
celloracle
anndata
scanpy
pandas
numpy
matplotlib
seaborn
tqdm
```

---

## Key Design Decisions

- **Shared input network**: Both pipelines consume the same `pruned_edges.csv`, ensuring results are directly comparable.
- **scTenifold** operates at the **network topology level** — it is fast, does not require single-cell data, and outputs statistically tested differential regulation scores.
- **CellOracle** operates at the **single-cell expression level** — it is computationally heavier but captures cell-type-specific regulatory dynamics and produces interpretable cell state transition visualizations.
- Both pipelines iterate over **all RBPs** in the network automatically, requiring no manual per-gene configuration.
