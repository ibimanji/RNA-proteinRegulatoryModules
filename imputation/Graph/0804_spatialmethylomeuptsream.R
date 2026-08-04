source("/home/jinpu/R_project/plot_in_R/R/spaimpute_skill_0718.R")
gc()
p21_rds_path = "/media/zenglab/result/lingyuan2/2025nature/P21_mCG_RNA_spatial/data/GSE270498_spatial.obj_P21.rds"
#3 main assays for the p21rds obj
P21_MCG_ASSAY = "DNAm_frac_CG"
P21_MCA_ASSAY = "DNAm_frac_CA"
P21_RNA_ASSAY = "RNA"

# Ordered gene[i] <-> mCG VMR[i] <-> mCA VMR[i] triplets.
P21_FEATURE_NAMES = "Prox1,Ntrk3,Cux1"
P21_MCG_FEATURES = paste(c(
  "chr1.190146228.190156728",
  "chr7.78477326.78479526",
  "chr5.136361742.136365142"
), collapse = ",")
P21_MCA_FEATURES = paste(c(
  "chr1.190147438.190160338",
  "chr7.78355114.78361914",
  "chr5.136299005.136310005"
), collapse = ",")
P21_RNA_FEATURES = "Prox1,Ntrk3,Cux1"

# 1. preprocessing --------------------------------------------------------

p21 = readRDS(p21_rds_path)
p21
coords <- GetTissueCoordinates(p21)
head(coords)
coords_grid <- coords[, c("imagecol","imagerow")]
xy <- as.matrix(coords_grid)
D <- dist(xy)
D <- as.vector(D)
D <- D[D > 0]
summary(D)
library(FNN)
nn <- get.knn(
  xy,
  k=10
)
head(nn$nn.dist)
nearest_dist <- nn$nn.dist[,1]
summary(nearest_dist)

#so we can see here that for the nearst neighbor, d = 17, which is the grid row-col length for square pattern.
#if we want the nearest neighbor for such lattice cell-cell neighbor, we ought sqrt(2)d = 1.41 x 18 = 25.42,also you can choose 25 or 30
obj = p21
rm(p21)

dims_n <- 30
k <- 50
prune_snn <- 1/15
n_trees <- 100

obj <- FindNeighbors(
  object = obj,
  reduction = "pca",
  dims = 1:dims_n,
  k.param = k,
  compute.SNN = TRUE,
  prune.SNN = prune_snn,
  nn.method = "annoy",
  n.trees = n_trees,
  annoy.metric = "euclidean",
  l2.norm = TRUE,
  graph.name = c(
    "RNA_new_nn",
    "RNA_new_snn"
  )
)
W_RNA_snn <- obj@graphs$RNA_new_snn

obj <- FindNeighbors(
  object = obj,
  reduction = "mcgpca",
  dims = 1:10,
  k.param = k,
  compute.SNN = TRUE,
  prune.SNN = prune_snn,
  nn.method = "annoy",
  n.trees = n_trees,
  annoy.metric = "euclidean",
  l2.norm = TRUE,
  graph.name = c(
    "mCG_new_nn",
    "mCG_new_snn"
  )
)
W_mCG_snn <- obj@graphs$mCG_new_snn


obj <- FindNeighbors(
  object = obj,
  reduction = "mcapca",
  dims = 1:10,
  k.param = k,
  compute.SNN = TRUE,
  prune.SNN = prune_snn,
  nn.method = "annoy",
  n.trees = n_trees,
  annoy.metric = "euclidean",
  l2.norm = TRUE,
  graph.name = c(
    "mCA_new_nn",
    "mCA_new_snn"
  )
)
W_mCA_snn <- obj@graphs$mCA_new_snn


W_spatial_original <- build_radius_graph(
  coords_grid,
  radius=26#[sqrt(2) x mindist.max] + 1
)
rownames(W_spatial_original)<-colnames(obj)
colnames(W_spatial_original)<-colnames(obj)

rna_counts = GetAssayData(obj, assay = 'RNA', layer = 'counts')
neighbor_expr <- NeighborExpressionEmbedding(
  W_spatial_original,
  expr=t(rna_counts),
  alpha=0.0
)

neighbor_expr_data <- t(neighbor_expr)
colnames(neighbor_expr_data) = colnames(obj)

obj[["Neighbor"]] <-
  CreateAssayObject(
    counts=as(
      neighbor_expr_data,
      "dgCMatrix"
    )
  )
DefaultAssay(obj)="Neighbor"


obj <- NormalizeData(
  obj,
  assay="Neighbor",scale.factor = median(colSums(neighbor_expr_data))
)


obj <- FindVariableFeatures(
  obj,
  assay="Neighbor",
  selection.method="vst",
  nfeatures=1500
)


obj <- ScaleData(
  obj,
  assay="Neighbor"
)


obj <- RunPCA(
  obj,
  assay="Neighbor",
  features=VariableFeatures(obj),
  npcs=50,
  reduction.name="neighbor.pca",
  reduction.key="NEIGHBORPC_"
)
p1 = ElbowPlot(
  obj,
  reduction = "pca",
  ndims = 50
)


p2 = ElbowPlot(
  obj,
  reduction = "neighbor.pca",
  ndims = 50
)
p1+p2

obj <- FindNeighbors(
  object=obj,
  reduction="neighbor.pca",
  dims=1:dims_n,
  k.param=k,
  compute.SNN=TRUE,
  prune.SNN=prune_snn,
  nn.method="annoy",
  n.trees=n_trees,
  annoy.metric="euclidean",
  l2.norm=TRUE,
  graph.name=c(
    "nb_new_nn",
    "nb_new_snn"
  )
)


W_neighbor_snn<- obj@graphs$nb_new_snn
W_mCA_snn = obj@graphs$DNAm_CA_snn


# 0. cellwise alpha calculation -------------------------------------------



graph_list_alpha <- list(
  RNA = W_RNA_snn,
  mCG = W_mCG_snn,
  Neighbor = W_neighbor_snn
)
Z_RNA <- Embeddings(
  obj,
  reduction="pca"
)[,1:30]

Z_mCG <- Embeddings(
  obj,
  reduction="mcgpca"
)[,1:10]
Z_neighbor <- Embeddings(
  obj,
  reduction="neighbor.pca"
)[,1:30]
NormalizeEmbedding <- function(Z){
  
  Z <- as.matrix(Z)
  
  Z / sqrt(rowSums(Z^2)+1e-10)
  
}
embedding_list_alpha <- list(
  RNA = Z_RNA,
  mCG = Z_mCG,
  Neighbor = Z_neighbor
)
embedding_list_alpha <- lapply(
  embedding_list_alpha,
  NormalizeEmbedding
)
CheckCellConsistency(
  graph_list_alpha,
  embedding_list_alpha
)
alpha_result <- ComputeCellwiseAlpha(
  graph_list = graph_list_alpha,
  embedding_list = embedding_list_alpha,
  lambda_self = 1,
  beta = 1
)
alpha_cell <- alpha_result$alpha
rownames(alpha_cell) <- rownames(W_RNA_snn)
alpha_result$alpha <- alpha_cell
alpha_cell[, "Neighbor"] <- 0
alpha_cell <-
  alpha_cell /
  rowSums(alpha_cell) *
  length(colnames(alpha_cell)) 

graph_list_scml <- list(
  RNA = W_RNA_snn,
  mCG = W_mCG_snn,
  Neighbor = W_neighbor_snn
)

# 2. RUN SCML CORE！ -------------------------------------------------------

scml_result <- run_scml_core_pan_cell_specificalpha(
  graph_list = graph_list_scml,
  ndim = 20,
  alpha = alpha_cell,
  layer_weight = c(
    RNA = 1,
    mCG = 1,
    Neighbor = 1
  )
)

ord <- order(scml_result$eigenvalues)
eigs <- scml_result$eigenvalues[ord]
eigs
SCML_embeddingraw <-
  scml_result$U
SCML_embedding <-
  scml_result$U[,ord,drop=FALSE]
rownames(SCML_embedding) <- rownames(W_RNA_snn)
colnames(SCML_embedding) <-
  paste0(
    "SCML_",
    seq_len(ncol(SCML_embedding))
  )
DefaultAssay(obj)<-"RNA"


obj[["scml"]] <-
  CreateDimReducObject(
    embeddings = SCML_embedding,
    key="SCML_",
    assay="RNA"
  )

p1 = hist(
  obj[["scml"]]@cell.embeddings[,"SCML_1"],
  breaks=50,
  main="SCML_1 distribution",
  xlab="SCML_1"
)

p2 = hist(
  obj[["scml"]]@cell.embeddings[,"SCML_2"],
  breaks=50,
  main="SCML_2 distribution",
  xlab="SCML_2"
)


obj <- FindNeighbors(
  obj,
  reduction="scml",
  dims=1:20,
  k.param=20,
  annoy.metric="euclidean",
  l2.norm = TRUE, 
  graph.name=c(
    "scml_nn",
    "scml_snn"
  )
)

obj <- FindClusters(
  obj,
  
  graph.name="scml_snn",
  
  resolution=0.5, algorithm = 3, 
  
  cluster.name="scml_cluster",
  
  random.seed=1
)

obj <- RunUMAP(
  obj,
  
  reduction="scml",
  
  dims=1:20,
  
  reduction.name="scml.umap",
  
  reduction.key="SCMLUMAP_",
  
  spread=1,
  
  min.dist=0.2,
  
  seed.use=1
)
p_umap <- DimPlot(
  obj,
  reduction = "scml.umap",
  group.by = "scml_cluster",
  label = TRUE,
  repel = TRUE,
  pt.size = 1
) +
  theme_classic()


spatial_df <- data.frame(
  x = coords$imagecol,
  y = coords$imagerow,
  cluster = obj$scml_cluster
)

rownames(spatial_df) <- rownames(coords)
spatial_df <- spatial_df[colnames(obj),]
library(Polychrome)

ncluster <- length(unique(spatial_df$cluster))

cluster_colors <- createPalette(
  ncluster,
  c("#000000","#FFFFFF")
)

names(cluster_colors) <- sort(unique(spatial_df$cluster))


p_spatial <- ggplot(
  spatial_df,
  aes(
    x=x,
    y=y,
    color=cluster
  )
)+
  geom_point(
    size=2
  )+
  scale_color_manual(
    values=cluster_colors
  )+
  scale_y_reverse()+
  coord_fixed()+
  theme_classic()


p_spatial+p_umap

alpha_RNA_df <- data.frame(
  cell = rownames(alpha_cell),
  alpha_RNA = alpha_cell[, "RNA"]
)


meta_df <- obj@meta.data %>%
  mutate(
    cell = rownames(.)
  )


plot_df <- alpha_RNA_df %>%
  left_join(
    meta_df[, c(
      "cell",
      "scml_cluster"
    )],
    by="cell"
  )



ggplot(
  plot_df,
  aes(
    x = scml_cluster,
    y = alpha_RNA,
    fill = scml_cluster
  )
)+
  geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.3,alpha = 1
  )+
  geom_boxplot(
    width = 0.12,
    outlier.size = 0.2,
    alpha = 0.9
  )+
  scale_fill_manual(
    values = cluster_colors
  )+
  theme_classic(
    base_size = 14
  )+
  labs(
    x = "SCML cluster",
    y = expression(alpha[RNA])
  )+
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    legend.position = "none"
  )

saveRDS(
  alpha_result,
  file = "/media/zenglab/data/jinpu/spaimpute/alpha_result0804.rds"
)
saveRDS(
  obj,
  file = "/media/zenglab/data/jinpu/spaimpute/spatialmethylome_result0804.rds"
)

# 3.Imputation ------------------------------------------------------------

W <- as(
  obj@graphs$scml_snn,
  "dgCMatrix"
)


spec <- ComputeGraphSpectrum2(
  W,
  n_eigs=300
)

saveRDS(
  spec,
  file = "/media/zenglab/data/jinpu/spaimpute/spectrum_result0804.rds"
)
#使用前进行：
L <- spec$L_sym
ord <- order(spec$lambda)
lambda <- spec$lambda[ord]
V <- spec$V_rw[,ord]
colnames(V) <- paste0(
  "v",
  seq_len(ncol(V))
)
lambda
summary(V[,1]) # should be all the same

#we also add function below for glance :

ComputeGraphSpectrum2 <- function(
    W,
    n_eigs = 300,
    transform_to_rw = TRUE,
    remove_zero_degree = TRUE
){
  
  ####################################
  ## degree
  ####################################
  
  D <- Matrix::rowSums(W)
  
  
  if(remove_zero_degree && any(D == 0)){
    
    stop(
      "Graph contains isolated nodes with degree zero."
    )
    
  }
  
  
  ####################################
  ## normalized Laplacian
  ####################################
  
  D_inv_sqrt <- Matrix::Diagonal(
    x = 1 / sqrt(D)
  )
  
  
  L_sym <- 
    Matrix::Diagonal(nrow(W)) -
    D_inv_sqrt %*%
    W %*%
    D_inv_sqrt
  
  
  
  ####################################
  ## eigendecomposition
  ####################################
  
  eig <- RSpectra::eigs_sym(
    L_sym,
    k = n_eigs,
    which = "SA"
  )
  
  
  V_sym <- eig$vectors
  
  lambda <- eig$values
  
  
  
  rownames(V_sym) <- rownames(W)
  
  
  
  ####################################
  ## transform eigenvectors
  ## L_rw space
  ####################################
  
  if(transform_to_rw){
    
    V_rw <-
      D_inv_sqrt %*%
      V_sym
    
    
    rownames(V_rw) <- rownames(W)
    
  }else{
    
    V_rw <- NULL
    
  }
  
  
  
  ####################################
  ## check DC component
  ####################################
  
  dc_check <- list(
    
    V_sym_first = V_sym[,1],
    
    V_rw_first = 
      if(transform_to_rw)
        V_rw[,1]
      else
        NULL,
    
    sqrt_degree =
      sqrt(D),
    
    constant =
      rep(1,nrow(W))
    
  )
  
  
  
  return(
    list(
      
      lambda = lambda,
      
      # normalized Laplacian eigenvectors
      V_sym = V_sym,
      
      # transformed eigenvectors
      # first vector approximately constant
      V_rw = V_rw,
      
      L_sym = L_sym,
      
      degree = D,
      
      dc_check = dc_check
      
    )
  )
}

ComputeNodeDirichletEnergy <- function(
    W,
    embedding,
    normalize = TRUE
){
  
  stopifnot(
    nrow(W)==nrow(embedding)
  )
  
  
  n <- nrow(W)
  
  
  deg <- Matrix::rowSums(W)
  
  
  Z <- as.matrix(embedding)
  
  
  energy <- numeric(n)
  
  
  for(i in seq_len(n)){
    
    nei <- which(W[i,] > 0)
    
    if(length(nei)==0){
      energy[i] <- 0
      next
    }
    
    
    diff <- sweep(
      Z[nei,,drop=FALSE],
      2,
      Z[i,],
      "-"
    )
    
    
    dist2 <- rowSums(diff^2)
    
    
    energy[i] <-
      0.5 *
      sum(
        W[i,nei] * dist2
      ) #energy equals 1/2 Wij * dist2(i,j)
    
  }
  
  
  if(normalize){
    
    energy <-
      energy/(deg+1e-10)
    
  } #node wise normalization
  
  
  return(energy)
}

ComputeCrossDirichletEnergy <- function(
    W_source,
    embedding_target,
    normalize = TRUE
){
  
  
  n <- nrow(W_source)
  
  
  deg <- Matrix::rowSums(W_source)
  
  
  Z <- as.matrix(
    embedding_target
  )
  
  
  energy <- numeric(n)
  
  
  for(i in seq_len(n)){
    
    
    nei <- which(
      W_source[i,] > 0
    )
    
    
    if(length(nei)==0){
      energy[i] <- 0
      next
    }
    
    
    diff <- sweep(
      Z[nei,,drop=FALSE],
      2,
      Z[i,],
      "-"
    )
    
    
    dist2 <- rowSums(diff^2)
    
    
    energy[i] <-
      0.5*
      sum(
        W_source[i,nei]*
          dist2
      )
    
  }
  
  
  if(normalize){
    
    energy <-
      energy/(deg+1e-10)
    
  }
  
  
  return(energy)
  
}

NormalizeEnergy <- function(E){
  
  E / (median(E) + 1e-10)
  
}

ComputeCellwiseAlpha <- function(
    graph_list,
    embedding_list,
    lambda_self = 1,
    beta = 1,
    normalize_alpha = TRUE
){
  
  
  layers <- names(graph_list)
  
  
  n <- nrow(
    graph_list[[1]]
  )
  
  
  M <- length(layers)
  
  
  score_mat <- matrix(
    0,
    nrow=n,
    ncol=M,
    dimnames=list(
      NULL,
      layers
    )
  )
  
  
  
  energy_self <- list()
  energy_self_norm <- list()
  
  
  energy_cross <- list()
  energy_cross_norm <- list()
  
  
  ################################
  # self energy
  ################################
  
  for(m in layers){
    
    energy_self[[m]] <-
      ComputeNodeDirichletEnergy(
        graph_list[[m]],
        embedding_list[[m]]
      )
    
    
    # modality-wise normalization
    energy_self_norm[[m]] <-
      NormalizeEnergy(
        energy_self[[m]]
      )
    
  }
  
  
  
  ################################
  # cross energy
  ################################
  
  for(m in layers){
    
    
    cross_total_norm <- rep(
      0,
      n
    )
    
    
    cross_total_raw <- rep(
      0,
      n
    )
    
    
    for(target in layers){
      
      if(target == m)
        next
      
      
      ################################
      # W_source from modality m
      # signal from target modality
      ################################
      
      cross_energy <-
        ComputeCrossDirichletEnergy(
          W_source =
            graph_list[[m]],
          embedding_target =
            embedding_list[[target]]
        )
      
      
      cross_total_raw <-
        cross_total_raw +
        cross_energy
      
      
      ################################
      # normalize according to target
      ################################
      
      cross_total_norm <-
        cross_total_norm +
        NormalizeEnergy(
          cross_energy
        )
      
      
    }
    
    
    energy_cross[[m]] <-
      cross_total_raw
    
    
    energy_cross_norm[[m]] <-
      cross_total_norm
    
    
    ################################
    # score
    ################################
    
    score_mat[,m] <-
      cross_total_norm -
      lambda_self *
      energy_self_norm[[m]]
    
    
  }
  
  
  
  ################################
  # score -> alpha
  ################################
  
  # lower score = better
  
  alpha <- exp(
    -beta *
      scale(score_mat)
  )
  
  
  alpha <-
    alpha /
    rowSums(alpha)
  
  
  
  ################################
  # avoid zero alpha
  ################################
  
  if(normalize_alpha){
    
    alpha[
      alpha < 1e-4
    ] <- 1e-4
    
    
    alpha <-
      alpha /
      rowSums(alpha)
    
  }
  
  graph_score_raw <- numeric(M)
  names(graph_score_raw) <- layers
  
  
  for(m in layers){
    
    # raw cross energy
    cross_raw <- energy_cross[[m]]
    
    # raw self energy
    self_raw <- energy_self[[m]]
    
    
    graph_score_raw[m] <-
      mean(cross_raw) -
      lambda_self * mean(self_raw)
    
  }
  
  graph_alpha_rawEnergy <- exp(
    -beta * scale(graph_score_raw)
  )
  
  
  graph_alpha_rawEnergy <-
    graph_alpha_rawEnergy /
    sum(graph_alpha_rawEnergy)
  
  
  
  return(
    list(
      
      alpha = alpha,
      
      graph_score_rawEnergy =
        graph_score_raw,
      
      graph_alpha_rawEnergy =
        graph_alpha_rawEnergy,
      
      self_energy =
        energy_self,
      
      self_energy_norm =
        energy_self_norm,
      
      cross_energy =
        energy_cross,
      
      cross_energy_norm =
        energy_cross_norm,
      
      score =
        score_mat
      
    )
  )
}

