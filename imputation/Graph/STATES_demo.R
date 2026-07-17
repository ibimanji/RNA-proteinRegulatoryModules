#this is for the demo of SCMLbased SPAIMPUTATION and function refine
#notice that this 'sampling' method may be risky since it is NOT suitable for the ZINB since we consider rare dropouts here.
#ALSO we can later see: scDesign2.

source("/home/jinpu/R_project/plot_in_R/R/spaimpute_skill_0708.R")
#!!!!!ALSO LOAD THE FUNCTION IN THE LAST PART. "APPE"
# 1.pseudomatrix-prepare --------------------------------------------------


ngenes <- 40
ncells <- 80

genes <- paste0("Gene",1:ngenes)
cells <- paste0("Cell",1:ncells)

celltype <- c(
  rep("OLG1",20),
  rep("OLG2",20),
  rep("MLG",20),
  rep("AST",20)
)


row <- numeric(ncells)
column <- numeric(ncells)



# ==========================
# OLG1
# random distribution
# ==========================

row[1:20] <-
  runif(20, 1, 40)

column[1:20] <-
  runif(20, 1, 40)



# ==========================
# OLG2
# spatial domain
# concentrated region
# ==========================

row[21:40] <-
  rnorm(
    20,
    mean = 32,
    sd = 3
  )

column[21:40] <-
  rnorm(
    20,
    mean = 32,
    sd = 3
  )


# 防止超出范围
row[21:40] <-
  pmin(
    pmax(row[21:40],20),
    40
  )

column[21:40] <-
  pmin(
    pmax(column[21:40],20),
    40
  )



# ==========================
# MLG
# random distribution
# ==========================

row[41:60] <-
  runif(20, 1, 40)

column[41:60] <-
  runif(20, 1, 40)



# ==========================
# AST
# random distribution
# ==========================

row[61:80] <-
  runif(20, 1, 40)

column[61:80] <-
  runif(20, 1, 40)



meta <- data.frame(
  row=row,
  column=column,
  replicate="B4",
  states_nn_alg1_label3=celltype,
  row.names=cells
)

OLG_marker <- c(
  "Gene1",
  "Gene2",
  "Gene3",
  "Gene4",
  "Gene5"
)
# OLG2 subtype marker
OLG2_marker <- c(
  "Gene6",
  "Gene7"
)
# MLG marker
MLG_marker <- c(
  "Gene8",
  "Gene9",
  "Gene10",
  "Gene11",
  "Gene12"
)
# AST marker
AST_marker <- c(
  "Gene13",
  "Gene14",
  "Gene15",
  "Gene16",
  "Gene17"
)

total_raw_mat <- matrix(
  0,
  ngenes,
  ncells,
  dimnames=list(
    genes,
    cells
  )
)


for(i in 1:ncells){
  
  type <- celltype[i]
  
  
  # background expression
  expr <- rgamma(
    ngenes,
    shape=1,
    rate=1
  )
  
  names(expr) <- genes
  
  if(type %in% c("OLG1","OLG2")){
    
    expr[OLG_marker] <-
      expr[OLG_marker] +
      rgamma(
        length(OLG_marker),
        shape=5,
        rate=1
      )
    
  }
  
  
  
  ###################################
  # OLG2 subtype program
  ###################################
  
  if(type=="OLG2"){
    
    expr[OLG2_marker] <-
      expr[OLG2_marker] +
      rgamma(
        length(OLG2_marker),
        shape=3,
        rate=1
      )
    
  }
  
  
  
  
  ###################################
  # MLG program
  ###################################
  
  if(type=="MLG"){
    
    expr[MLG_marker] <-
      expr[MLG_marker] +
      rgamma(
        length(MLG_marker),
        shape=8,
        rate=1
      )
    
  }
  
  
  
  ###################################
  # AST program
  ###################################
  
  if(type=="AST"){
    
    expr[AST_marker] <-
      expr[AST_marker] +
      rgamma(
        length(AST_marker),
        shape=8,
        rate=1
      )
    
  }
  
  
  
  lib_size <- rpois(
    1,
    5000
  )
  
  
  total_raw_mat[,i] <-
    rmultinom(
      1,
      size=lib_size,
      prob=expr/sum(expr)
    )
  
}

TE <- matrix(
  rbeta(
    ngenes*ncells,
    4,
    6
  ),
  nrow=ngenes
)

rb_raw_mat <- matrix(
  0,
  ngenes,
  ncells,
  dimnames=list(
    genes,
    cells
  )
)



for(i in 1:ncells){
  
  rb_raw_mat[,i] <-
    rbinom(
      ngenes,
      size=total_raw_mat[,i],
      prob=TE[,i]
    )
  
}


nt_raw_mat <-
  total_raw_mat-rb_raw_mat


states_demo <- CreateSeuratObject(
  counts=as(total_raw_mat,"dgCMatrix"),
  meta.data=meta
)

#Glance1:spatialMap
ggplot(
  states_demo@meta.data,
  aes(
    column,
    row,
    color=states_nn_alg1_label3
  )
)+
  geom_point(size=4)+
  coord_fixed()+
  scale_y_reverse()+
  theme_classic()

#Glance2: Expression HeatMap
states_demo <- NormalizeData(states_demo)
cell_order_type <- c(
  "OLG1",
  "OLG2",
  "MLG",
  "AST"
)
cell_order <- order(
  factor(
    states_demo$states_nn_alg1_label3,
    levels = cell_order_type
  )
)


gene_order <- c(
  # OLG shared
  "Gene1",
  "Gene2",
  "Gene3",
  "Gene4",
  "Gene5",
  
  # OLG2 specific
  "Gene6",
  "Gene7",
  
  # MLG marker
  "Gene8",
  "Gene9",
  "Gene10",
  "Gene11",
  "Gene12",
  
  # AST marker
  "Gene13",
  "Gene14",
  "Gene15",
  "Gene16",
  "Gene17",
  
  # remaining genes
  paste0(
    "Gene",
    18:40
  )
)


expr <- GetAssayData(
  states_demo,
  assay="RNA",
  layer="data"
)


expr <- expr[
  gene_order,
  cell_order
]

expr <- t(
  scale(
    t(
      as.matrix(expr)
    )
  )
)


cell_annotation <- states_demo$states_nn_alg1_label3[
  cell_order
]


ha <- HeatmapAnnotation(
  CellType = cell_annotation
)
celltype_colors <- c(
  OLG1 = "#4C78A8",
  OLG2 = "#72B7B2",
  MLG  = "#F58518",
  AST  = "#54A24B"
)


ha <- HeatmapAnnotation(
  CellType = celltype,
  col = list(
    CellType = celltype_colors
  ),
  annotation_name_gp = gpar(
    fontsize = 12
  )
)

expr_col <- colorRamp2(
  c(-2,0,2),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)

column_split <- factor(
  celltype,
  levels=c(
    "OLG1",
    "OLG2",
    "MLG",
    "AST"
  )
)

gene_split <- c(
  rep("OLG",5),
  rep("OLG2",2),
  rep("MLG",5),
  rep("AST",5),
  rep("BG",23)
)

Heatmap(
  expr,
  name="Expression\n(Z-score)",
  
  col=expr_col,
  
  top_annotation=ha,
  
  column_split=column_split,
  row_split=gene_split,
  
  cluster_columns=FALSE,
  cluster_rows=FALSE,
  
  show_column_names=FALSE,
  show_row_names=TRUE,
  
  row_names_gp=gpar(
    fontsize=9
  ),
  
  column_title="Simulated spatial transcriptomics",
  
  heatmap_legend_param=list(
    title_gp=gpar(fontsize=12),
    labels_gp=gpar(fontsize=10)
  ),
  
  border=TRUE,
  
  gap=unit(2,"mm")
)
states_demo[["rbRNA"]] <-
  CreateAssayObject(
    counts=as(rb_raw_mat,"dgCMatrix")
  )


states_demo[["ntRNA"]] <-
  CreateAssayObject(
    counts=as(nt_raw_mat,"dgCMatrix")
  )
DefaultAssay(states_demo) <- "ntRNA"


states_demo <- NormalizeData(
  states_demo,
  assay="ntRNA"
)


nt_expr <- GetAssayData(
  states_demo,
  assay="ntRNA",
  layer="data"
)



nt_expr <- nt_expr[
  gene_order,
  cell_order
]


nt_expr <- t(
  scale(
    t(
      as.matrix(nt_expr)
    )
  )
)



Heatmap(
  nt_expr,
  name="ntRNA\nZ-score",
  
  col=expr_col,
  
  top_annotation=ha,
  
  column_split=column_split,
  row_split=gene_split,
  
  cluster_columns=FALSE,
  cluster_rows=FALSE,
  
  show_column_names=FALSE,
  show_row_names=TRUE,
  
  column_title="ntRNA marker expression",
  
  border=TRUE,
  gap=unit(2,"mm")
)
states_demo


# 2.SCML glance -----------------------------------------------------------

library(igraph)
library(ggraph)
plot_graph_weight <- function(
    W,
    meta,
    title="Graph",
    ncell_show=80,
    use_spatial=FALSE
){
  set.seed(123)
  
  
  #----------------------------
  # ensure cell names
  #----------------------------
  
  if(is.null(colnames(W))){
    colnames(W) <- rownames(meta)
    rownames(W) <- rownames(meta)
  }
  
  
  cells_use <- sample(
    colnames(W),
    min(
      ncell_show,
      ncol(W)
    )
  )
  
  
  W_sub <- W[
    cells_use,
    cells_use
  ]
  
  
  meta_sub <- meta[
    cells_use,
    ,
    drop=FALSE
  ]
  
  
  
  #----------------------------
  # adjacency -> igraph
  #----------------------------
  
  g <- igraph::graph_from_adjacency_matrix(
    W_sub,
    mode="undirected",
    weighted=TRUE
  )
  
  
  #----------------------------
  # choose layout
  #----------------------------
  
  if(use_spatial){
    
    # real spatial coordinate
    
    layout_df <- data.frame(
      x=meta_sub$column,
      y=meta_sub$row
    )
    
    rownames(layout_df) <- cells_use
    
    layout_df <- layout_df[
      V(g)$name,
      ,
      drop=FALSE
    ]
    
    
    lay <- as.matrix(layout_df)
    
  }else{
    
    # graph topology layout
    
    lay <- igraph::layout_with_fr(
      g
    )
    
  }
  
  
  
  #----------------------------
  # plot
  #----------------------------
  
  ggraph(
    g,
    layout="manual",
    x=lay[,1],
    y=lay[,2]
  )+
    geom_edge_link(
      aes(
        width=weight,
        alpha=weight
      )
    )+
    
    geom_node_point(
      aes(
        color=
          meta_sub[
            V(g)$name,
            "states_nn_alg1_label3"
          ]
      ),
      size=4
    )+
    
    {
      if(use_spatial){
        scale_y_reverse()
      }
    }+
    
    coord_fixed()+
    
    theme_void()+
    
    ggtitle(title)
}




DefaultAssay(states_demo)<-"RNA"


states_demo <- NormalizeData(
  states_demo,
  normalization.method="RC",
  scale.factor=
    median(
      colSums(
        total_raw_mat
      )
    )
)



states_demo <- FindVariableFeatures(
  states_demo,
  nfeatures=30
)



states_demo <- ScaleData(
  states_demo
)



states_demo <- RunPCA(
  states_demo,
  npcs=10,
  reduction.name="rna.pca"
)



rb_counts <- GetAssayData(
  states_demo,
  assay="rbRNA",
  layer="counts"
)


rb_norm <- sweep(
  rb_counts,
  2,
  colSums(rb_counts),
  "/"
)*
  median(
    colSums(rb_counts)
  )



states_demo <- SetAssayData(
  states_demo,
  assay="rbRNA",
  layer="data",
  new.data=
    as(
      rb_norm,
      "dgCMatrix"
    )
)



states_demo <- ScaleData(
  states_demo,
  assay="rbRNA"
)

states_demo <- FindVariableFeatures(
  states_demo,
  assay="rbRNA",
  selection.method="vst",
  nfeatures=30
)

states_demo <- RunPCA(
  states_demo,
  assay="rbRNA",
  npcs=10,
  reduction.name="rbrna.pca"
)



nt_counts <- GetAssayData(
  states_demo,
  assay="ntRNA",
  layer="counts"
)



nt_norm <- sweep(
  nt_counts,
  2,
  colSums(nt_counts),
  "/"
)*
  median(
    colSums(nt_counts)
  )



states_demo <- SetAssayData(
  states_demo,
  assay="ntRNA",
  layer="data",
  new.data=
    as(
      nt_norm,
      "dgCMatrix"
    )
)


states_demo <- FindVariableFeatures(
  states_demo,
  assay="ntRNA",
  selection.method="vst",
  nfeatures=30
)

states_demo <- ScaleData(
  states_demo,
  assay="ntRNA"
)

states_demo <- FindVariableFeatures(
  states_demo,
  nfeatures=30
)

states_demo <- RunPCA(
  states_demo,
  assay="ntRNA",
  npcs=10,
  reduction.name="ntrna.pca"
)


Emb_rb <- Embeddings(
  states_demo,
  "rbrna.pca"
)[,1:10]


Emb_nt <- Embeddings(
  states_demo,
  "ntrna.pca"
)[,1:10]



W_rb <- build_knn_graph(
  emb=Emb_rb,
  k=10,
  sigma_nn=5,
  symmetrize="max"
)
rownames(W_rb) <- colnames(states_demo)
colnames(W_rb) <- colnames(states_demo)

p_rb = plot_graph_weight(
  W_rb,
  states_demo@meta.data,
  "rbRNA graph",
  use_spatial=FALSE
)

p_rb





W_nt <- build_knn_graph(
  emb=Emb_nt,
  k=10,
  sigma_nn=5,
  symmetrize="max"
)
rownames(W_rb) <- colnames(states_demo)
colnames(W_rb) <- colnames(states_demo)

p_nt <- plot_graph_weight(
  W_nt,
  states_demo@meta.data,
  "ntRNA graph",
  use_spatial=FALSE
)


p_nt





coords <- states_demo@meta.data[
  ,
  c(
    "column",
    "row"
  )
]


coords <- as.matrix(coords)

rownames(coords)<-
  colnames(states_demo)

coords <- states_demo@meta.data[,c("column","row")]

D <- as.matrix(dist(coords))

hist(
  D[D>0],
  breaks=50,
  main="Pairwise spatial distance",
  xlab="Distance"
)


W_spatial_orginal <- build_radius_graph(
  coords=coords,radius = 10,
  sigma_floor=1e-3
)#可以改成delaunay
rownames(W_spatial_orginal) <- colnames(states_demo)
colnames(W_spatial_orginal) <- colnames(states_demo)
p_spatial_ori <- plot_graph_weight(
  W_spatial_orginal,
  states_demo@meta.data,
  "original radius graph, \n r = 10",use_spatial=TRUE
)


p_spatial_ori

#2 methods for spatial enhance: (1) RNA layer enhance (2)considering spatial-neighborhood expressing
W_spatial_refine = SpatialEnhanceRNAgraph_sparse(W_rna = W_rb,
                                                 W_spatial_original = W_spatial_orginal,
                                                 symmetrize = 'max',lambda = 2)

p_spatial <- plot_graph_weight(
  W_spatial_refine,
  states_demo@meta.data,
  "Adaptive spatial graph_refined",use_spatial=FALSE
)
p_rb|p_spatial


neighbor_expr <- NeighborExpressionEmbedding(
  W_spatial_orginal,
  expr = t(rb_norm),
  alpha=0.1
)

neighbor_pca <- prcomp(
  neighbor_expr,
  scale.=TRUE
)$x[,1:10]

W_neighbor <- build_knn_graph(k=10,sigma_nn=5,neighbor_pca)
p_spatial <- plot_graph_weight(
  W_neighbor,
  states_demo@meta.data,
  "Neighborhood graph_refined",use_spatial=FALSE)
p_spatial


############################################################
#SCML
############################################################
scml_result <- run_scml_core( #存在scml_core_new,就是选取从"SM"改为"SA",Lmodified可能存在负数特征值
  graph_list=list(
    ntRNA=W_nt,
    rbRNA=W_rb,
    spatial=W_neighbor
  ),
  ndim=10,
  default_alpha=1,
  spatial_alpha=0.5
)

scml_result$eigenvalues

L_scml <- scml_result$L_modified
U_embedding = scml_result$U
rownames(U_embedding)
W_U <- build_knn_graph(
  emb=U_embedding,
  k=10,
  sigma_nn=5,
  symmetrize="max",
)
rownames(W_U) <- colnames(states_demo)
colnames(W_U) <- colnames(states_demo)

p_U <- plot_graph_weight(
  W_U,
  states_demo@meta.data,
  "U graph",
  use_spatial=FALSE
)
p_rb | p_U
rownames(W_U)

# 3. Uembedding UMAP scanning ----------------------------------------------


dim(U_embedding)
rownames(U_embedding) <- colnames(states_demo)

DefaultAssay(states_demo) <- "RNA"
states_demo <- FindNeighbors(
  states_demo,
  reduction="rna.pca",
  dims=1:10,
  graph.name="RNA_snn"
)
states_demo <- FindClusters(
  states_demo,
  graph.name="RNA_snn",
  resolution=0.5
)
states_demo <- RunUMAP(
  states_demo,
  reduction="rna.pca",
  dims=1:10,
  reduction.name="rna.umap"
)
p1 <- DimPlot(
  states_demo,
  reduction="rna.umap",
  group.by="states_nn_alg1_label3"
)
p1

states_demo[["scml"]] <-
  CreateDimReducObject(
    embeddings = U_embedding,
    key="SCML_",
    assay="RNA"
  )
states_demo <- FindNeighbors(
  states_demo,
  reduction="scml",
  dims=1:ncol(U_embedding),
  k.param=10,
  graph.name="SCML_snn"
)
states_demo <- FindClusters(
  states_demo,
  graph.name="SCML_snn",
  resolution=0.5
)
states_demo <- RunUMAP(
  states_demo,
  reduction="scml",
  dims=1:ncol(U_embedding),
  reduction.name="scml.umap"
)
p2 <- DimPlot(
  states_demo,
  reduction="scml.umap",
  group.by="states_nn_alg1_label3"
)
p1 + p2


# 4. IMPUTATION -----------------------------------------------------------

W <- BuildSCMLGraph(U_embedding, k=10, sigma_nn = 5)

spec <- ComputeGraphSpectrum(
  W,
  n_eigs=10
)

lambda <- spec$lambda
lambda
V <- spec$V
ord <- order(lambda)
lambda_sort <- lambda[ord]
V_sort <- V[, ord, drop = FALSE]
g_heat <- GraphFilterKernel(
  lambda_sort,
  type = "heat",beta = 10,
  plot = TRUE
)

rb_raw_mat = GetAssayData(states_demo, layer = 'counts', assay = 'rbRNA')
total_raw_mat = GetAssayData(states_demo, layer  = 'counts', assay = 'RNA')
nt_raw_mat = GetAssayData(states_demo, layer = 'counts', assay = 'ntRNA')

total_norm <- GetAssayData(states_demo, layer = 'data',assay = 'RNA')
rb_norm <- GetAssayData(states_demo, layer = 'data',assay = 'rbRNA')
nt_norm <- GetAssayData(states_demo, layer = 'data',assay = 'ntRNA')


rb_raw_f <- GraphWaveletDenoise(
  expr = rb_raw_mat,
  V = V_sort,
  g = g_heat
)
colnames(rb_raw_f) <- colnames(rb_raw_mat)
all.equal(rownames(rb_raw_f),rownames(rb_raw_mat))


nt_raw_f <- GraphWaveletDenoise(
  expr = nt_raw_mat,
  V = V_sort,
  g = g_heat
)
colnames(nt_raw_f) <- colnames(nt_raw_mat)
all.equal(rownames(nt_raw_f),rownames(nt_raw_mat))

total_raw_f <- GraphWaveletDenoise(
  expr = total_raw_mat,
  V = V_sort,
  g = g_heat
)
colnames(total_raw_f) <- colnames(total_raw_mat)
all.equal(rownames(total_raw_f),rownames(total_raw_mat))


rb_norm_f <- GraphWaveletDenoise(
  expr = rb_norm,
  V = V_sort,
  g = g_heat
)
colnames(rb_norm_f) <- colnames(rb_norm)
all.equal(
  rownames(rb_norm_f),
  rownames(rb_norm)
)


nt_norm_f <- GraphWaveletDenoise(
  expr = nt_norm,
  V = V_sort,
  g = g_heat
)
colnames(nt_norm_f) <- colnames(nt_norm)
all.equal(
  rownames(nt_norm_f),
  rownames(nt_norm)
)


total_norm_f <- GraphWaveletDenoise(
  expr = total_norm,
  V = V_sort,
  g = g_heat
)
colnames(total_norm_f) <- colnames(total_norm)
all.equal(
  rownames(total_norm_f),
  rownames(total_norm)
)


states_demo[["total_wavelet"]] <-
  CreateAssayObject(
    counts = total_raw_f
  )

states_demo[["total_wavelet"]] <-
  SetAssayData(
    states_demo[["total_wavelet"]],
    layer = "data",
    new.data = total_norm_f
  )

states_demo[["rb_wavelet"]] <-
  CreateAssayObject(
    counts = rb_raw_f
  )

states_demo[["rb_wavelet"]] <-
  SetAssayData(
    states_demo[["rb_wavelet"]],
    layer="data",
    new.data=rb_norm_f
  )


states_demo[["nt_wavelet"]] <-
  CreateAssayObject(
    counts = nt_raw_f
  )

states_demo[["nt_wavelet"]] <-
  SetAssayData(
    states_demo[["nt_wavelet"]],
    layer="data",
    new.data=nt_norm_f
  )

DefaultAssay(states_demo) <- "rb_wavelet"


FeaturePlot(
  states_demo,
  features=c("Gene6","Gene1","Gene11","Gene17","Gene29","Gene30"),
  reduction="scml.umap"
)

rb_raw <- GetAssayData(
  states_demo,
  assay="rbRNA",
  layer="data"
)

nt_raw <- GetAssayData(
  states_demo,
  assay="ntRNA",
  layer="data"
)
TE_cell_raw <-
  Matrix::colSums(rb_raw) /
  (
    Matrix::colSums(rb_raw) +
      Matrix::colSums(nt_raw) +
      1e-8
  )

rb_wave <- GetAssayData(
  states_demo,
  assay="rb_wavelet",
  layer="data"
)


nt_wave <- GetAssayData(
  states_demo,
  assay="nt_wavelet",
  layer="data"
)
TE_cell_wave <-
  Matrix::colSums(rb_wave) /
  (
    Matrix::colSums(rb_wave) +
      Matrix::colSums(nt_wave) +
      1e-8
  )
states_demo$TE_raw <- TE_cell_raw
states_demo$TE_wave <- TE_cell_wave

p3 <- FeaturePlot(
  states_demo,
  features = "TE_raw",
  reduction = "scml.umap",
  cols = c("lightgrey", "darkred")
) +
  ggtitle("Raw TE")
p4 <- FeaturePlot(
  states_demo,
  features = "TE_wave",
  reduction = "scml.umap",
  cols = c("lightgrey", "darkred")
) +
  ggtitle("Wavelet TE")
p3 + p4
