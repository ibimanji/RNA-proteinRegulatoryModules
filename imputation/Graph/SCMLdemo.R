#0. packages needed and data loading, mention that FNN and Rspectra are also needed for analysis
rm(list = ls()); gc(); options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixStats)
  library(anndata)
  library(reticulate)
  library(ggplot2)
  library(ggrepel)
})
states = read_h5ad("/media/zenglab/data/jinpu/0529statesmouse14month_ADWT/14mWT_14mAD_glanced.h5ad")
rb_raw_mat <- t(states$layers[["rbRNA_raw"]]) 
total_raw_mat <- t(states$layers[["totalRNA_raw"]])
nt_raw_mat <- t(states$layers[["ntRNA_raw"]])
meta = states$obs
states_seu <- CreateSeuratObject(
  counts = total_raw_mat,
  meta.data = meta,
  project = "statesmouse"
)
states_seu[["rbRNA"]] <- CreateAssayObject(
  counts = rb_raw_mat
)
states_seu[["ntRNA"]] <- CreateAssayObject(
  counts = nt_raw_mat
)
all_umap_coords <- as.matrix(states$obsm[['X_umap']])
colnames(all_umap_coords) <- c("UMAP_1", "UMAP_2")
rownames(all_umap_coords) <- rownames(states$obs)
umap_coords <- all_umap_coords[, , drop = FALSE]
umap_reduction <- CreateDimReducObject(
  embeddings = umap_coords,
  key = "UMAP_",       # 坐标轴的前缀，生成 UMAP_1 和 UMAP_2
  assay = "RNA"        # 关联的主 Assay
)
states_seu[["umap"]] <- umap_reduction
DimPlot(states_seu, reduction = "umap", group.by = "celltype")
#修复：1，NA修复，2 Leiden6 CA Leiden7 DG
meta_df <- states_seu@meta.data
na_clusters <- meta_df %>%
  filter(is.na(celltype) | celltype == "NA") %>%
  pull(leiden) %>%
  unique()
print(na_clusters)#只来自于leiden3，这是我当时处理的时候leiden3没细分动都没动直接分为了OLG-TEPN mix，但设置category的时候有一步问题，导致NA的出现，修复即可，也可以画UMAP进行验证
states_seu@meta.data <- states_seu@meta.data %>%
  mutate(celltype = case_when(
    leiden %in% na_clusters & (is.na(celltype) | celltype == "NA") ~ "OLG-TEPN mix", #这三在这里一回事
    TRUE ~ as.character(celltype) # 其余的保持原样
  ))
table(states_seu$celltype, useNA = "ifany")
DimPlot(states_seu, reduction = "umap", group.by = "celltype")

states_seu@meta.data$celltype2 <- states_seu@meta.data$celltype
states_seu@meta.data$celltype2[states_seu@meta.data$leiden == 6] <- "CA"
states_seu@meta.data$celltype2[states_seu@meta.data$leiden == 7] <- "DG"
DimPlot(states_seu, reduction = "umap", group.by = "celltype2")
states_seu@meta.data$celltype <- states_seu@meta.data$celltype2
DimPlot(states_seu, reduction = "umap", group.by = "celltype") 

states_seu <- subset(
  states_seu,
  subset = type == "14mWT"
)
#magic trail via R
#1.1 we focus on the WT here, to see if it will return the linear biological events(cell cycle related ones， especially the UBE2C、TPX2 related ones)
#normalize should be carefully considered，here we shall use totalRNA layer to calculate spike-in-factor, to evaluate the cell sequencing-lib
#With this sequencing-lib, we can use them to normalize rbRNA layer (and ntRNA layer)

#totalRNA layer as spike in factor
rna_counts <- GetAssayData(states_seu, assay = "RNA", layer = "counts")
rna_lib_sizes <- colSums(rna_counts) #as spike-in factor for other layers

#totalRNA layer: which can be done on default seurat function pipeline
DefaultAssay(states_seu) <- "RNA"
states_seu <- NormalizeData(states_seu, normalization.method = "RC", scale.factor = median(rna_lib_sizes))
states_seu <- FindVariableFeatures(states_seu, selection.method = "vst", nfeatures = 2000)
states_seu <- ScaleData(states_seu)
states_seu <- RunPCA(states_seu, npcs = 50, reduction.name = "rna.pca", reduction.key = "rnaPC_")

#rbRNA layer
rb_counts <- GetAssayData(states_seu, assay = "rbRNA", layer = "counts")
scale_factor_ref_rb <- median(colSums(GetAssayData(states_seu, assay = "rbRNA", layer = "counts")))
rb_norm <- sweep(rb_counts, 2, rna_lib_sizes, "/") * scale_factor_ref_rb
states_seu <- SetAssayData(states_seu, assay = "rbRNA", layer = "data", new.data = as(rb_norm, "dgCMatrix"))
states_seu <- FindVariableFeatures(states_seu, selection.method = "vst", nfeatures = 2000, assay = "rbRNA")
states_seu <- ScaleData(states_seu, assay = "rbRNA")
states_seu <- RunPCA(states_seu, npcs = 50, assay = "rbRNA", reduction.name = "rbrna.pca", reduction.key = "rbrnaPC_")

#ntRNA layer
nt_counts <- GetAssayData(states_seu, assay = "ntRNA", layer = "counts")
scale_factor_ref_nt <- median(colSums(GetAssayData(states_seu, assay = "ntRNA", layer = "counts")))
nt_norm <- sweep(nt_counts, 2, rna_lib_sizes, "/") * scale_factor_ref_nt
states_seu <- SetAssayData(states_seu, assay = "ntRNA", layer = "data", new.data = as(nt_norm, "dgCMatrix"))
states_seu <- FindVariableFeatures(states_seu, selection.method = "vst", nfeatures = 2000, assay = "ntRNA")
states_seu <- ScaleData(states_seu, assay = "ntRNA")
states_seu <- RunPCA(states_seu, npcs = 50, assay = "ntRNA", reduction.name = "ntrna.pca", reduction.key = "ntrnaPC_")

Emb_spatial <- states_seu@meta.data[,c("column","row")]
Emb_spatial <- as.matrix(Emb_spatial)
rownames(Emb_spatial) <- colnames(states_seu)

library(FNN)
library(RSpectra)
BuildKNNGraph <- function(
    emb,
    k = 30
){
  knn <- get.knn(
    emb,
    k = k
  ) #FNN:Fast k-nearest neighbor searching algorithms,can return indice and dist-Eu
  idx <- knn$nn.index
  dst <- knn$nn.dist
  n <- nrow(idx)
  sigma <- dst[,k] #an adaptive sigma (cell-wise), can be changed
  sigma[sigma == 0] <- 1e-6
  ii <- rep(
    seq_len(n),
    each = k
  ) 
  jj <- as.vector(
    t(idx)
  ) #flatten idx
  d <- as.vector(
    t(dst)
  )
  sigma_i <- rep(
    sigma,
    each = k
  )
  sigma_j <- sigma[jj]
  w <- exp(
    - d^2 /
      (
        sigma_i *
          sigma_j
      )
  )# adaptive Gaussian Kernel
  W <- sparseMatrix(
    i = ii,
    j = jj,
    x = w,
    dims = c(n,n)
  )
  W <- pmax(
    W,
    t(W)
  )#weight below should be symmetric, that is to say the graph should be no-direction
  return(W)
}
ComputeLaplacian <- function(W){
  deg <- Matrix::rowSums(W)#D matrix
  Dinv <- Diagonal(
    x = 1/sqrt(deg + 1e-10)
  )#D-1/2
  L <- Diagonal(nrow(W)) -
    Dinv %*% W %*% Dinv #this is normalized Laplace 
  return(L)
}
GetSubspace <- function(
    L,
    ndim = 30
){
  eig <- eigs_sym(
    L,
    k = ndim,
    which = "SM"
  )#smallest eigs
  U <- eig$vectors
  U <- qr.Q(qr(U))
  return(U)#for a Grassmann(k, n), k 
}
RunSCML <- function(
    graph_list,
    ndim = 30,
    alpha = 0.5
){
  n <- nrow(graph_list[[1]])
  L_sum <- Matrix(
    0,
    n,
    n,
    sparse = TRUE
  )
  P_sum <- Matrix(
    0,
    n,
    n,
    sparse = TRUE
  )
  
  U_layers <- list()
  
  for(i in seq_along(graph_list)){
    
    cat(
      "Processing layer:",
      i,
      "\n"
    )
    L <- ComputeLaplacian(
      graph_list[[i]]
    )
    U <- GetSubspace(
      L,
      ndim
    )
    
    U_layers[[i]] <- U
    
    L_sum <- L_sum + L
    
    P_sum <- P_sum + tcrossprod(U)#U X UT
  }
  
  L_mod <- L_sum -
    alpha * P_sum
  
  eig <- eigs_sym(
    L_mod,
    k = ndim,
    which = "SM"
  )
  
  U_consensus <- eig$vectors
  
  U_consensus <- sweep(
    U_consensus,
    1,
    sqrt(rowSums(U_consensus^2)),
    "/"
  )
  
  return(
    list(
      U = U_consensus,
      U_layers = U_layers,
      L_mod = L_mod
    )
  )
}
```
batch function sealing
```{r}
RunSCMLPipeline <- function(
    Emb_nt,
    Emb_rb,
    Emb_spatial,
    states_seu,
    k = 30,
    ndim = 30,
    alpha = 1,
    resolution = 0.5
){
  
  cat(
    "\n==============================\n",
    "k =", k,
    " alpha =", alpha,
    "\n==============================\n"
  )
  
  W_nt <- BuildKNNGraph(
    Emb_nt,
    k = k
  )
  
  W_rb <- BuildKNNGraph(
    Emb_rb,
    k = k
  )
  
  W_spatial <- BuildKNNGraph(
    Emb_spatial,
    k = k
  )
  
  scml_res <- RunSCML(
    graph_list = list(
      W_nt,
      W_rb,
      W_spatial
    ),
    ndim = ndim,
    alpha = alpha
  )
  
  U <- scml_res$U
  
  rownames(U) <- colnames(states_seu)
  
  colnames(U) <- paste0(
    "SCML_",
    seq_len(ncol(U))
  )
  
  seu <- states_seu
  
  seu[["scml"]] <- CreateDimReducObject(
    embeddings = U,
    key = "SCML_",
    assay = "RNA"
  )
  
  seu <- FindNeighbors(
    seu,
    reduction = "scml",
    dims = 1:ndim,
    graph.name = "scml.nn"
  )
  
  seu <- FindClusters(
    seu,
    graph.name = "scml.nn",
    resolution = resolution
  )
  
  seu <- RunUMAP(
    seu,
    reduction = "scml",
    dims = 1:ndim,
    reduction.name = "scml.umap"
  )
  
  return(seu)
}
```
plot function seal
```{r}
RunSCMLPipeline <- function(
    Emb_nt,
    Emb_rb,
    Emb_spatial,
    states_seu,
    k = 30,
    ndim = 30,
    alpha = 1,
    resolution = 0.5
){
  
  cat(
    "\n==============================\n",
    "k =", k,
    " alpha =", alpha,
    "\n==============================\n"
  )
  
  ##################################################
  # Build graphs
  ##################################################
  
  W_nt <- BuildKNNGraph(Emb_nt, k = k)
  W_rb <- BuildKNNGraph(Emb_rb, k = k)
  W_spatial <- BuildKNNGraph(Emb_spatial, k = k)
  
  ##################################################
  # SCML
  ##################################################
  
  scml_res <- RunSCML(
    graph_list = list(W_nt, W_rb, W_spatial),
    ndim = ndim,
    alpha = alpha
  )
  
  U <- scml_res$U
  
  rownames(U) <- colnames(states_seu)
  
  colnames(U) <- paste0("SCML_", seq_len(ncol(U)))
  
  seu <- states_seu
  
  seu[["scml"]] <- CreateDimReducObject(
    embeddings = U,
    key = "SCML_",
    assay = "RNA"
  )
  
  seu <- FindNeighbors(
    seu,
    reduction = "scml",
    dims = 1:ndim,
    graph.name = "scml.nn"
  )
  
  seu <- FindClusters(
    seu,
    graph.name = "scml.nn",
    resolution = resolution
  )
  
  seu <- RunUMAP(
    seu,
    reduction = "scml",
    dims = 1:ndim,
    reduction.name = "scml.umap"
  )
  
  ##################################################
  # PLOT 1: celltype UMAP
  ##################################################
  
  p1 <- DimPlot(
    seu,
    reduction = "scml.umap",
    group.by = "celltype"
  ) +
    ggtitle(paste0("celltype | k=", k, " alpha=", alpha))
  
  print(p1)
  
  ##################################################
  # PLOT 2: cluster UMAP
  ##################################################
  
  p2 <- DimPlot(
    seu,
    reduction = "scml.umap",
    group.by = "seurat_clusters",
    label = TRUE
  ) +
    ggtitle(paste0("cluster | k=", k, " alpha=", alpha))
  
  print(p2)
  
  ##################################################
  # PLOT 3: spatial plot
  ##################################################
  
  df <- seu@meta.data
  
  p3 <- ggplot(
    df,
    aes(column, row, color = seurat_clusters)
  ) +
    geom_point(size = 0.4) +
    scale_y_reverse() +
    theme_classic() +
    ggtitle(paste0("spatial | k=", k, " alpha=", alpha))
  
  print(p3)
  
  ##################################################
  # return object
  ##################################################
  
  return(seu)
}
```

```{r}
Emb_nt      <- Embeddings(states_seu[["ntrna.pca"]])[,1:30]
Emb_rb      <- Embeddings(states_seu[["rbrna.pca"]])[,1:30]
Emb_spatial <- states_seu@meta.data[,c("column","row")]
Emb_spatial <- as.matrix(Emb_spatial)
rownames(Emb_spatial) <- colnames(states_seu)

###############################################################################
##parameters scanning: if not wished, change it to the ideal paras with len = 1
###############################################################################
k_grid <- c(
  15,
  30,
  50
)

alpha_grid <- c(
  0,
  0.25,
  0.5,
  1,
  2
)

for(k in k_grid){
  
  for(alpha in alpha_grid){
    
    seu_scml <- RunSCMLPipeline(
      Emb_nt = Emb_nt,
      Emb_rb = Emb_rb,
      Emb_spatial = Emb_spatial,
      states_seu = states_seu,
      k = k,
      ndim = 30,
      alpha = alpha,
      resolution = 0.5
    )
  }
}
