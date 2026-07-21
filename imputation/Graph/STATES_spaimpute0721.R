#you should load functions in spaimpute_skill_0718.R also we change C1 into B4
source("/home/jinpu/R_project/plot_in_R/R/spaimpute_skill_0718.R")
load("/media/zenglab/data/jinpu/statesB4/states_B4.RData")

# 1.preprocessing ---------------------------------------------------------
DefaultAssay(states_seu) <- "RNA"
rna_counts <- GetAssayData(
  states_seu,
  assay = "RNA",
  layer = "counts"
)
rna_lib_sizes <- colSums(rna_counts)
states_seu <- NormalizeData(
  states_seu,
  normalization.method = "RC",
  scale.factor = median(rna_lib_sizes)
)
states_seu <- FindVariableFeatures(
  states_seu,
  selection.method = "vst",
  nfeatures = 1500
)
states_seu <- ScaleData(states_seu)
states_seu <- RunPCA(
  states_seu,
  npcs = 50,
  reduction.name = "rna.pca",
  reduction.key = "rnaPC_"
)
total_norm = GetAssayData(states_seu, assay = 'RNA', layer = 'data')

cat("Running rbRNA normalization...\n")
rb_counts <- GetAssayData(
  states_seu,
  assay = "rbRNA",
  layer = "counts"
)
scale_factor_ref_rb <- median(
  colSums(rb_counts)
)
rb_norm <- sweep(
  rb_counts,
  2,
  rna_lib_sizes,
  "/"
) * scale_factor_ref_rb
states_seu <- SetAssayData(
  states_seu,
  assay = "rbRNA",
  layer = "data",
  new.data = as(rb_norm, "dgCMatrix")
)
states_seu <- FindVariableFeatures(
  states_seu,
  assay = "rbRNA",
  selection.method = "vst",
  nfeatures = 1500
)
states_seu <- ScaleData(
  states_seu,
  assay = "rbRNA"
)
states_seu <- RunPCA(
  states_seu,
  assay = "rbRNA",
  npcs = 50,
  reduction.name = "rbrna.pca",
  reduction.key = "rbrnaPC_"
)


cat("Running ntRNA normalization...\n")
nt_counts <- GetAssayData(
  states_seu,
  assay = "ntRNA",
  layer = "counts"
)
scale_factor_ref_nt <- median(
  colSums(nt_counts)
)
nt_norm <- sweep(
  nt_counts,
  2,
  rna_lib_sizes,
  "/"
) * scale_factor_ref_nt
states_seu <- SetAssayData(
  states_seu,
  assay = "ntRNA",
  layer = "data",
  new.data = as(nt_norm, "dgCMatrix")
)
states_seu <- FindVariableFeatures(
  states_seu,
  assay = "ntRNA",
  selection.method = "vst",
  nfeatures = 1500
)
states_seu <- ScaleData(
  states_seu,
  assay = "ntRNA"
)
states_seu <- RunPCA(
  states_seu,
  assay = "ntRNA",
  npcs = 50,
  reduction.name = "ntrna.pca",
  reduction.key = "ntrnaPC_"
)

Emb_rb <- Embeddings(
  states_seu,
  reduction = "rbrna.pca"
)[,1:30]
Emb_nt <- Embeddings(
  states_seu,
  reduction = "ntrna.pca"
)[,1:30]

cat("Building rbRNA graph...\n")
states_seu <- FindNeighbors(
  object = states_seu,
  reduction = "rbrna.pca",
  dims = 1:30,
  k.param = 30,
  compute.SNN = TRUE,
  prune.SNN = 1/15,
  nn.method = "annoy",
  n.trees = 100,
  annoy.metric = "euclidean",
  l2.norm = TRUE,
  graph.name = c(
    "rbRNA_nn",
    "rbRNA_snn"
  ),
  verbose = TRUE
)
W_rb <- states_seu@graphs$rbRNA_snn
W_rb <- as(W_rb, "dgCMatrix")
#max(abs(W_rb - t(W_rb)))

cat("Building ntRNA graph...\n")
states_seu <- FindNeighbors(
  object = states_seu,
  reduction = "ntrna.pca",
  dims = 1:30,
  k.param = 30,
  compute.SNN = TRUE,
  prune.SNN = 1/15,
  nn.method = "annoy",
  n.trees = 100,
  annoy.metric = "euclidean",
  l2.norm = TRUE,
  graph.name = c(
    "ntRNA_nn",
    "ntRNA_snn"
  ),
  verbose = TRUE
)
W_nt <- states_seu@graphs$ntRNA_snn
W_nt <- as(W_nt, "dgCMatrix")
dim(W_nt)
max(abs(W_nt - t(W_nt)))

cat("Building spatial graph...\n")
coords <- states_seu@meta.data[
  c("column","row")
]
coords <- as.matrix(coords)
rownames(coords)
rownames(coords)<-
  colnames(states_seu)

gc()

W_spatial_orginal <- build_radius_graph(
  coords=coords,radius = 2000,
  sigma_floor=1e-3
)
rownames(W_spatial_orginal) <- colnames(states_seu)
colnames(W_spatial_orginal) <- colnames(states_seu)

total_mat = GetAssayData(states_seu, assay = 'RNA',layer = "counts")
neighbor_expr <- NeighborExpressionEmbedding(
  W_spatial_orginal,
  expr = t(total_mat),
  alpha=0.2
)
dim(neighbor_expr)
neighbor_expr_data = t(neighbor_expr)
colnames(neighbor_expr_data) = colnames(rb_norm)
all.equal(rownames(neighbor_expr_data),rownames(rb_norm))

neighbor_assay <- CreateAssayObject(
  counts = neighbor_expr_data
)
states_seu[["Neighbor"]] <- neighbor_assay
DefaultAssay(states_seu) <- "Neighbor"
states_seu = NormalizeData(states_seu,
                           assay = 'Neighbor',
                           normalization.method = "RC",
                           scale.factor = median(colSums(neighbor_expr_data)))
states_seu <- FindVariableFeatures(
  states_seu,
  assay="Neighbor",
  selection.method="vst",
  nfeatures= 1500
)
states_seu <- ScaleData(
  states_seu,
  assay="Neighbor"
)
states_seu <- RunPCA(
  states_seu,
  assay="Neighbor",
  features=VariableFeatures(states_seu),
  npcs=50,
  reduction.name="neighbor.pca",
  reduction.key="NEIGHBORPC_"
)

states_seu <- FindNeighbors(
  object = states_seu,
  reduction = "neighbor.pca",
  dims = 1:30,
  k.param = 30,
  compute.SNN = TRUE,
  prune.SNN = 1/15,
  nn.method = "annoy",
  n.trees = 100,
  annoy.metric = "euclidean",
  l2.norm = TRUE,
  graph.name = c(
    "nb_nn",
    "nb_snn"
  ),
  verbose = TRUE
)
W_neighbor <- states_seu@graphs$nb_snn
W_neighbor <- as(W_neighbor, "dgCMatrix")
dim(W_neighbor)
max(abs(W_neighbor - t(W_neighbor)))

stopifnot(
  nrow(W_neighbor) == ncol(states_seu)
)
stopifnot(
  nrow(W_rb) == ncol(states_seu)
)
stopifnot(
  nrow(W_nt) == ncol(states_seu)
)


# 2.SCML_core -------------------------------------------------------------
cat("Running SCML_core")
scml_result <- run_scml_core_pan(
  graph_list = list(
    ntRNA   = W_nt,
    rbRNA   = W_rb,
    spatial = W_neighbor
  ),alpha = c(ntRNA = 1,rbRNA = 1 , spatial = 0), 
  layer_weight = c(ntRNA = 1,rbRNA = 1, spatial = 0.5),
  ndim = 30
)
D <- ComputeGrassmannDistance(
  U_consensus = scml_result$U,
  U_layers = scml_result$U_layers
)
D
gc()
library(tidyverse)
scml_result$eigenvalues_layers

scml_result$eigenvalues
SCML_embedding <- scml_result$U
L_scml <- scml_result$L_modified 
rownames(SCML_embedding) <- colnames(states_seu)
DefaultAssay(states_seu) <- "RNA"
states_seu[["scml"]] <- CreateDimReducObject(
  embeddings = SCML_embedding,
  key = "SCML_",
  assay = "RNA"
)
states_seu <- FindNeighbors(
  states_seu,
  reduction = "scml",
  dims = 1:30,
  k.param = 30,
  graph.name = c(
    "scml_nn",
    "scml_snn"
  )
)

states_seu <- FindClusters(
  states_seu, graph.name = 'scml_snn',
  resolution = 0.4
)

states_seu <- RunUMAP(
  states_seu,
  reduction = "scml",reduction.name = "scml.umap",
  dims = 1:25, spread = 0.5, min.dist = 0.05
)

# 3.plot for Uembedding ---------------------------------------------------


p1 = DimPlot(
  states_seu,
  reduction = "scml.umap",
  label = TRUE
)
p2 = DimPlot(
  states_seu,
  reduction = "scml.umap",group.by = 'mousebrain_label2',
  label = FALSE
)
p1|p2

#
library(scales)
library(Polychrome)
meta <- states_seu@meta.data
meta$seurat_clusters <- factor(meta$seurat_clusters)
cluster_levels <- levels(meta$seurat_clusters)
ncluster <- length(cluster_levels)
cols <- DiscretePalette(
  n = ncluster,
  palette = "polychrome"
)

names(cols) <- levels(meta$seurat_clusters)
ggplot(
  meta,
  aes(
    x = -column,
    y = row,
    color = seurat_clusters
  )
) +
  geom_point(size = 0.8) +
  coord_fixed() +
  scale_color_manual(values = cols) +
  theme_classic() +
  labs(color = "Cluster")



cluster_levels <- levels(meta$seurat_clusters)
xlim_range <- range(-meta$column, na.rm = TRUE)
ylim_range <- range(meta$row, na.rm = TRUE)
plots <- list()
for (cl in cluster_levels) {
  p <- ggplot(
    meta,
    aes(
      x = -column,
      y = row
    )
  ) +
    geom_point(
      color = "white",
      size = 0.2
    ) +
    geom_point(
      data = meta[meta$seurat_clusters == cl, ],
      color = cols[cl],
      size = 0.4
    ) +
    coord_fixed(
      xlim = xlim_range,
      ylim = ylim_range
    ) +
    
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(
        hjust = 0.5,
        size = 11
      ),
      plot.margin = margin(2,2,2,2)
    ) +
    
    labs(
      title = paste0("Cluster ", cl)
    )
  
  plots[[cl]] <- p
}

combined_plot <- wrap_plots(
  plots,
  nrow = 7,
  ncol = 3
)
combined_plot


# 4. Imputation -----------------------------------------------------------

rb_raw_mat = GetAssayData(states_seu, layer = 'counts', assay = 'rbRNA')
total_raw_mat = GetAssayData(states_seu, layer  = 'counts', assay = 'RNA')
nt_raw_mat = GetAssayData(states_seu, layer = 'counts', assay = 'ntRNA')

total_norm <- GetAssayData(states_seu, layer = 'data',assay = 'RNA')
rb_norm <- GetAssayData(states_seu, layer = 'data',assay = 'rbRNA')
nt_norm <- GetAssayData(states_seu, layer = 'data',assay = 'ntRNA')

W <- build_knn_graph(SCML_embedding, k=30)

spec <- ComputeGraphSpectrum(
  W,
  n_eigs=30
)

lambda <- spec$lambda
V <- spec$V
ord <- order(lambda)

lambda <- lambda[ord]
V <- V[, ord]
lambda <- lambda[-1]
V <- V[, -1]

lambda_sort <- sort(lambda)

g_heat <- GraphFilterKernel(
  lambda,
  type = "heat",beta = 20,
  plot = TRUE
)

##rb filtered
GraphSpectralDenoise <- function(
    expr,
    V,
    g
){
  
  g <- g
  
  X <- t(as.matrix(expr))
  
  X_hat <- t(V) %*% X
  
  X_hat_filt <- g * X_hat
  
  X_filt <- V %*% X_hat_filt
  
  X_filt[X_filt < 0] <- 0
  
  t(X_filt)
}
rb_filtered <- GraphSpectralDenoise(
  expr = rb_norm,
  V = V,
  g = g_heat
)
colnames(rb_filtered) <- colnames(rb_norm)
all.equal(rownames(rb_filtered),rownames(rb_norm))

##nt filtered
nt_filtered <- GraphSpectralDenoise(
  expr = nt_norm,
  V = V,
  g = g_heat
)
colnames(nt_filtered) <- colnames(nt_norm)
all.equal(rownames(nt_filtered),rownames(nt_norm))

##totalRNA filtered
total_filtered <- GraphSpectralDenoise(
  expr = total_norm,
  V = V,
  g = g_heat
)
colnames(total_filtered) <- colnames(total_norm)

states_seu[["rb_filtered"]] <-
  CreateAssayObject(
    counts =  rb_filtered
  )
states_seu[["nt_filtered"]] <-
  CreateAssayObject(
    counts = nt_filtered
  )
states_seu[["total_filtered"]] <-
  CreateAssayObject(
    counts = total_filtered
  )

te_raw <- rb_raw_mat / (total_raw_mat + 1e-6)

te_norm <- rb_filtered / (rb_filtered + nt_filtered + 1e-6)

PlotSpatialGene2 <- function(
    seu,
    gene,
    assay,
    layer = "data",
    vmax,
    title = NULL,
    point_size = 0.10,
    n_colors = 256
) {
  
  # 1. 获取表达数据并应用上限截断 (保留原本不共享色标的逻辑)
  expr <- GetAssayData(
    seu,
    assay = assay,
    layer = layer
  )[gene, ]
  
  expr_plot <- pmin(expr, vmax)
  
  # 2. 构建绘图数据框
  df <- data.frame(
    spatial_x = -seu$column,
    spatial_y = seu$row,
    Value = expr_plot
  )
  
  # 获取坐标极值用于等比例计算及坐标轴限制
  x_min <- min(df$spatial_x, na.rm = TRUE)
  x_max <- max(df$spatial_x, na.rm = TRUE)
  y_min <- min(df$spatial_y, na.rm = TRUE)
  y_max <- max(df$spatial_y, na.rm = TRUE)
  
  ar <- (y_max - y_min) / (x_max - x_min)
  
  # 3. 细腻的红色渐变色标 (同步目标函数)
  reds <- colorRampPalette(
    c(
      "grey", #可以改成grey
      "#fee8c8",
      "#fdbb84",
      "#fc8d59",
      "#d7301f",
      "#7f0000"
    )
  )(n_colors)
  
  # 默认标题
  if (is.null(title)) {
    title <- sprintf("Spatial expression of %s", gene)
  }
  
  # 4. ggplot 绘图
  p <- ggplot(df, aes(x = spatial_x, y = spatial_y, color = Value)) +
    geom_point(size = point_size) +
    
    # 坐标轴缩放与留白 (同步目标函数: expand = 0.02)
    scale_x_continuous(limits = c(x_min, x_max), expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(y_min, y_max), expand = expansion(mult = 0.02)) +
    
    # 独立渐变色标配置
    scale_color_gradientn(
      colours = reds,
      limits = c(0, vmax),
      breaks = c(0, vmax),
      labels = c(
        sprintf("%.2f", 0),
        sprintf("%.2f", vmax)
      ),
      guide = guide_colorbar(
        title = "Signal intensity",
        title.position = "top",
        barwidth = 0.6,   # 移除原本的 unit() 包裹，使用目标函数的简洁数字
        barheight = 6
      )
    ) +
    
    # 主题设定 (theme_bw + 细节微调)
    theme_bw(base_size = 12) +
    
    theme(
      aspect.ratio = ar,
      panel.grid = element_blank(),
      
      # 同步目标函数：轴标题增加恰当的 Margin，避免文字紧贴坐标轴
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6)),
      
      # 保持原代码的清爽：去掉不必要的数字刻度，突出空间结构
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      
      plot.title = element_text(hjust = 0.5, face = "bold")
      # 注：已移除原代码中的 legend.title = element_blank()，以确保 Signal intensity 正常显示
    ) +
    
    # 坐标轴与图例名称 (同步目标函数)
    labs(
      title = title,
      x = "Column ", 
      y = "Row", 
      color = NULL
    )
  
  return(p)
}
GetVmax <- function(mat, gene, qmax = 0.99){
  
  quantile(
    mat[gene, ],
    probs = qmax,
    na.rm = TRUE
  )
}

genes_use <- c("Prkcd","Sparcl1")


plot_total <- list()

for(gene in genes_use){
  
  vmax_raw  <- GetVmax(total_raw_mat, gene)
  vmax_wave <- GetVmax(GetAssayData(states_seu, assay = 'total_filtered',layer = 'counts'), gene)
  
  plot_total[[length(plot_total)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "RNA",layer = 'counts',
      vmax = vmax_raw,
      title = paste0(gene, "\nTotal Raw")
    )
  
  plot_total[[length(plot_total)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "total_filtered",layer = 'counts',
      vmax = vmax_wave,
      title = paste0(gene, "\nTotal Filtered")
    )
}

p_total <- wrap_plots(plot_total, ncol = 2) +
  plot_annotation(title = "Total RNA (ntRNA + rbRNA)")
p_total

plot_rb <- list()

for(gene in genes_use){
  
  vmax_raw  <- GetVmax(rb_raw_mat, gene)
  vmax_wave <- GetVmax(GetAssayData(states_seu, assay = 'rb_filtered',layer = 'counts'), gene)
  
  plot_rb[[length(plot_rb)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "rbRNA",layer = 'counts',
      vmax = vmax_raw,
      title = paste0(gene, "\nRb Raw")
    )
  
  plot_rb[[length(plot_rb)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "rb_filtered",layer = 'counts',
      vmax = vmax_wave,
      title = paste0(gene, "\nRb Filtered")
    )
}

p_rb <- wrap_plots(plot_rb, ncol = 2) +
  plot_annotation(title = "rbRNA")
p_rb

