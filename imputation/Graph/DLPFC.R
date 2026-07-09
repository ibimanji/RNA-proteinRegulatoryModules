library(patchwork)
source("/home/jinpu/R_project/plot_in_R/R/spaimpute_skill_0708.R") #function loading
plot_dlpfc_map <- function(df, label_col, title, palette, show_legend = TRUE) {
  ggplot(df, aes(x = imagecol, y = -imagerow, color = .data[[label_col]])) +
    geom_point(size = 2.0, stroke = 0) +
    scale_color_manual(values = palette, na.value = "grey85", drop = FALSE) +
    coord_fixed() +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = 15),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.key.height = grid::unit(0.42, "cm"),
      legend.key.width = grid::unit(0.42, "cm"),
      plot.margin = margin(6, 6, 6, 6)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(title = title)
} #目前还没用
FindClustersTarget <- function(
    seu,
    target_clusters = 7,
    graph.name = "scml_snn",
    algorithm = 4,
    resolution.range = seq(0.05, 2, by = 0.01),
    random.seed = 12,
    verbose = FALSE
){
  
  best_diff <- Inf
  best_obj <- NULL
  best_resolution <- NA
  best_ncluster <- NA
  
  for(res in resolution.range){
    
    tmp <- suppressMessages(
      FindClusters(
        object = seu,
        graph.name = graph.name,
        algorithm = algorithm,
        resolution = res,
        random.seed = random.seed,
        verbose = verbose
      )
    )
    
    pred <- as.character(tmp$seurat_clusters)
    
    ncluster <- length(unique(pred))
    
    diff <- abs(ncluster - target_clusters)
    
    if(diff < best_diff){
      
      best_diff <- diff
      best_obj <- tmp
      best_resolution <- res
      best_ncluster <- ncluster
      
    }
    
    ## 完全命中直接返回
    if(ncluster == target_clusters){
      
      return(list(
        seu = tmp,
        resolution = res,
        ncluster = ncluster
      ))
      
    }
    
  }
  
  warning(
    paste0(
      "Cannot obtain exactly ",
      target_clusters,
      " clusters. Closest = ",
      best_ncluster
    )
  )
  
  list(
    seu = best_obj,
    resolution = best_resolution,
    ncluster = best_ncluster
  )
  
}# cluster into specific number of cluster
dlpfc_rdata_path <- "/media/zenglab/result/lingyuan/STEM/imputation/10x/data/spatialLIBD/Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata"
load(dlpfc_rdata_path) #load sce
dlpfc_sample <- "151673"
sce <- sce[, sce$sample_name == dlpfc_sample]
cat("Sample", dlpfc_sample, "\n")
cat("Dimension:", dim(sce), "\n")

counts <- assay(sce, "counts")
meta <- as.data.frame(colData(sce))
rownames(meta) <- colnames(sce)

dlpfc_seu <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  project = paste0("DLPFC_", dlpfc_sample)
)
dlpfc_seu <- subset(
  dlpfc_seu,
  subset = !is.na(layer_guess_reordered)
)#删除那些NA的type

DefaultAssay(dlpfc_seu) <- "RNA"
rna_counts <- get_assay_matrix(dlpfc_seu, assay = "RNA", layer = "counts")
rna_lib_sizes <- colSums(rna_counts)
cat("Running RNA normalization and PCA...\n")
dlpfc_seu <- NormalizeData(
  dlpfc_seu,
  normalization.method = "RC",
  scale.factor = median(rna_lib_sizes)
)
dlpfc_seu <- FindVariableFeatures(dlpfc_seu, selection.method = "vst", nfeatures = 2000)#这个Variable feature一般默认2000即可
dlpfc_seu <- ScaleData(dlpfc_seu)
dlpfc_seu <- RunPCA(
  dlpfc_seu,
  npcs = 50,
  reduction.name = "rna.pca",
  reduction.key = "rnaPC_"
)
dlpfc_ndim = 30 # n_dims选择的影响可见spaimpute-slide的第36页，这里的dim选择是用来构图的，跟runSCML的n_dims相比可以比他大也可以比他小都行
Emb_rna <- Embeddings(dlpfc_seu, reduction = "rna.pca")[, seq_len(dlpfc_ndim), drop = FALSE]

cat("Building RNA graph...\n")
W_rna <- build_knn_graph(
  emb = Emb_rna,
  k = 30, symmetrize = 'max', sigma_nn = 18
) # 这里使用的是Jaccard-based的KNN构图计算weight的方式了，基础的adaptive-sigma kernel可见build_knn_graph_tradition
#这里的选择影响的超参数检验可以见spaimpute-slide的第35页

cat("Building spatial graph...\n")
if (all(c("imagecol", "imagerow") %in% colnames(meta))) {
  coords <- as.matrix(meta[colnames(dlpfc_seu), c("imagecol", "imagerow")])
} else if (all(c("col", "row") %in% colnames(meta))) {
  coords <- as.matrix(meta[colnames(dlpfc_seu), c("col", "row")])
} else {
  coords <- spatialCoords(spe)
  coords <- coords[colnames(dlpfc_seu), , drop = FALSE]
}
rownames(coords) <- colnames(dlpfc_seu)
dim(coords)
coords # 获取coords
W_spatial <- build_delaunay_graph(coords) #对于10X visium也可以KNN，但是更高分辨率的ST就别用KNN了

W_spatial_refine <- ReweightSpatialGraph(
  W_spatial = W_spatial,
  Emb_rb = Emb_rna,
  Emb_nt = Emb_rna,#没得nt，rb区分，直接传入两次Emb_rna，理论上没问题（可见该函数的function，就是把Emb_rb跟Emb_nt col binds在一起然后算similarity）
  center = 0.34, #0.34是cosine similarity的一个分位数数值
  beta = 3,
  gate_quantile = 0.1, #10X visium不用删太多边，但对于高分辨率的，很多空间边无益于cell annotation and information diffusion，可以多删去一些
  symmetrize = "max"
) #因为这里还不是自适应，所以传入了center，该函数运行完会给予一个类似报告的内容，输出在console内，可以查看cosine-similarity的分布，然后联合sigmoid函数特质，确定center以及beta
#center，beta以及自适应后beta_scale, 见slides
############################################################
## Graph comparison: only for glance here 可以跳过
############################################################
deg0 <- Matrix::rowSums(W_rna)
deg1 <- Matrix::rowSums(W_spatial)
deg2 <- Matrix::rowSums(W_spatial_refine)

deg_mat <- cbind(
  RNA = deg0,
  Spatial = deg1,
  Spatial_refine = deg2
)

print(cor(deg_mat, method = "spearman"))

cat("Original edges :", length(W_spatial@x), "\n")
cat("Refined edges  :", length(W_spatial_refine@x), "\n")
cat(
  "Retention      :",
  round(
    length(W_spatial_refine@x) /
      length(W_spatial@x) * 100,
    2
  ),
  "%\n"
)

L0 <- Matrix::Diagonal(x = Matrix::rowSums(W_spatial)) - W_spatial
L1 <- Matrix::Diagonal(x = Matrix::rowSums(W_spatial_refine)) - W_spatial_refine

eig0 <- eigen(as.matrix(L0), only.values = TRUE)$values
eig1 <- eigen(as.matrix(L1), only.values = TRUE)$values

print(cor(eig0, eig1))

############################################################
## SCML
############################################################

stopifnot(nrow(W_rna) == ncol(dlpfc_seu))
stopifnot(nrow(W_spatial_refine) == ncol(dlpfc_seu))
cat("Running SCML with RNA + spatial graphs...\n")

scml_result <- run_scml_core(
  graph_list = list(
    RNA = W_rna,
    spatial = W_spatial_refine #这里选择W_spatial还是W_spatial_refine可以看slide的第33页，讨论在第34页
  ),
  ndim = 50,
  default_alpha = 1.00,
  spatial_alpha = 0.2 #这个alpha的超参数检验可以见slide的第32页
)

SCML_embedding <- scml_result$U
rownames(SCML_embedding) <- colnames(dlpfc_seu)

dlpfc_seu[["scml"]] <- CreateDimReducObject(
  embeddings = SCML_embedding,
  key = "SCML_",
  assay = "RNA"
)
dlpfc_seu <- FindNeighbors(
  dlpfc_seu,
  reduction = "scml",
  dims = 1:50,
  graph.name = "scml_snn",
  k.param = 20 #k的超参数检验可以看slide的第32页
)
dlpfc_seu <- FindClusters(
  dlpfc_seu,
  graph.name = "scml_snn",
  algorithm = 4,          # Leiden
  resolution = 0.54, #这里最好要调整一下，我也写了一个指定聚类到7类的function，可以考虑替换
  random.seed = 12,
  verbose = FALSE
)

############################################################
## ARI
############################################################
pred <- as.character(
  dlpfc_seu$seurat_clusters
)
unique(pred)
pred <- as.character(
  dlpfc_seu$seurat_clusters
)
truth <- dlpfc_seu$layer_guess_reordered #无NA之乱耳就很好了，当时有NA这里乱得很

library(mclust)
ari <- adjustedRandIndex(
  as.character(truth),
  pred
)
cat(
  sprintf(
    "SCML Leiden ARI = %.4f\n",
    ari
  )
)



############################################################
## PLOTTING for glance：颜色尚未对应
############################################################
layer_levels <- c(
  "Layer1",
  "Layer2",
  "Layer3",
  "Layer4",
  "Layer5",
  "Layer6",
  "WM",
  "NA"
)


layer_palette <- c(
  Layer1 = "#9ecae1",
  Layer2 = "#1f78b4",
  Layer3 = "#b2df8a",
  Layer4 = "#31a354",
  Layer5 = "#fb9a99",
  Layer6 = "#e31a1c",
  WM = "#fdbf6f",
  `NA` = "grey60"
  
)
truth_plot <- as.character(truth)
truth_plot[is.na(truth_plot)] <- "NA"#  这里其实没有NA了只是一个保险


cluster_levels <- sort(
  unique(pred)
)
cluster_palette <- c(
  
  "1" = "#9ecae1",
  "2" = "#1f78b4",
  "3" = "#b2df8a",
  "4" = "#31a354",
  "5" = "#fb9a99",
  "6" = "#e31a1c",
  "7" = "#fdbf6f",
  "8" = "grey60"
  
)


plot_df <- data.frame(
  barcode = colnames(dlpfc_seu),
  imagecol = as.numeric(
    dlpfc_seu$imagecol
  ),
  imagerow = as.numeric(
    dlpfc_seu$imagerow
  ),
  reference = factor(
    truth_plot,
    levels = layer_levels
  ),
  scml = factor(
    pred,
    levels = cluster_levels
  ),
  stringsAsFactors = FALSE
  
)



############################################################
## Reference plot
############################################################

p_ref <- plot_dlpfc_map(
  
  plot_df,
  
  label_col = "reference",
  
  title = paste0(
    "10x Visium human DLPFC\n",
    "sample no. ",
    unique(dlpfc_seu$sample_name),
    "\nreference annotations"
  ),
  
  palette = layer_palette,
  
  show_legend = TRUE
  
)



############################################################
## SCML cluster plot
############################################################

p_scml <- plot_dlpfc_map(
  
  plot_df,
  
  label_col = "scml",
  
  title = sprintf(
    "SCML + Leiden clusters\nARI = %.3f",
    ari
  ),
  
  palette = cluster_palette,
  
  show_legend = TRUE
  
)

p <- p_ref + p_scml +
  
  plot_layout(
    widths = c(1.2,1.2)
  )

print(p)
