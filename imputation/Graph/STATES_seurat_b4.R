#you should load functions in spaimpute_skill_0718.R also we change C1 into B4
source("/home/jinpu/R_project/plot_in_R/R/spaimpute_skill_0718.R")
load("/media/zenglab/data/jinpu/statesB4/states_B4.RData")
library(tidyverse)
library(scales)
library(Polychrome)  #plotting color

# 1.preprocessing ---------------------------------------------------------
#for expression assays(RNA, ntRNA rbRNA), normalization and scale, PCA and snnGraph construction via seurat
cat("Running totalRNA normalization...\n")
DefaultAssay(states_seu) <- "RNA"

rna_counts <- GetAssayData(
  states_seu,
  assay = "RNA",
  layer = "counts"
)

rna_lib_sizes <- colSums(rna_counts)

states_seu <- NormalizeData(
  states_seu,
  normalization.method = "RC", # Relative counts. Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor 
  scale.factor = median(rna_lib_sizes)
)

rna_norm <- GetAssayData(
  states_seu,
  assay = "RNA",
  layer = "data"
)

rna_norm <- log1p(rna_norm) #Log(RC_normalized + 1), for relatively stable variance

states_seu <- SetAssayData(
  states_seu,
  assay = "RNA",
  layer = "data",
  new.data = rna_norm
)


states_seu <- FindVariableFeatures(
  states_seu,
  assay = "RNA",
  selection.method = "vst", #variance stablize transform
  nfeatures = 1500
) # 如果仅仅使用了RC而且没有log1p，这一步需要注意，这一步决定了哪些基因进PCA，现在加log1p这一步其实怎么做都还好

states_seu <- ScaleData(
  states_seu,
  assay = "RNA"
)

states_seu <- RunPCA(
  states_seu,
  assay = "RNA",
  npcs = 50,
  reduction.name = "rna.pca",
  reduction.key = "rnaPC_"
)


total_norm <- GetAssayData(
  states_seu,
  assay = "RNA",
  layer = "data"
)



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


# log1p
rb_norm <- log1p(rb_norm)


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


# log1p
nt_norm <- log1p(nt_norm)


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
)[,1:25]
Emb_nt <- Embeddings(
  states_seu,
  reduction = "ntrna.pca"
)[,1:25]

cat("Building rbRNA graph...\n")
states_seu <- FindNeighbors(
  object = states_seu,
  reduction = "rbrna.pca",
  dims = 1:25,
  k.param = 50,
  compute.SNN = TRUE,
  prune.SNN = 1/20, #抹杀小于1/20的权重
  nn.method = "annoy",
  n.trees = 100,
  annoy.metric = "euclidean", 
  l2.norm = TRUE,
  graph.name = c(
    "rbRNA_nn",
    "rbRNA_snn"
  ), #构图需要给俩名字，不然它默认给你一个knn图，weight都是1的那种knn
  verbose = TRUE
)
W_rb <- states_seu@graphs$rbRNA_snn
W_rb <- as(W_rb, "dgCMatrix")
max(abs(W_rb - t(W_rb))) #检查是否对称

cat("Building ntRNA graph...\n")
states_seu <- FindNeighbors(
  object = states_seu,
  reduction = "ntrna.pca",
  dims = 1:25,
  k.param = 50,
  compute.SNN = TRUE,
  prune.SNN = 1/20,
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
  coords=coords,radius = 2000, #半径，与coords保持一致
  sigma_floor=1e-3
)
rownames(W_spatial_orginal) <- colnames(states_seu)
colnames(W_spatial_orginal) <- colnames(states_seu)

total_mat = GetAssayData(states_seu, assay = 'RNA',layer = "counts") #拿的是totalRNA的count
neighbor_expr <- NeighborExpressionEmbedding(
  W_spatial_orginal,
  expr = t(total_mat),
  alpha=0.4 #这里的alpha参数为”自身表达的一个占比权重，alpha = 0代表着完全由邻居决定”，可以自行调整，甚至可以设置成0
)
dim(neighbor_expr)
neighbor_expr_data = t(neighbor_expr)
colnames(neighbor_expr_data) = colnames(total_mat)
all.equal(rownames(neighbor_expr_data),rownames(total_mat))

neighbor_assay <- CreateAssayObject(
  counts = neighbor_expr_data
)
states_seu[["Neighbor"]] <- neighbor_assay
DefaultAssay(states_seu) <- "Neighbor"
states_seu = NormalizeData(states_seu,
                           assay = 'Neighbor')
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
  dims = 1:25,
  k.param = 50,
  compute.SNN = TRUE,
  prune.SNN = 1/20,
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
  ),alpha = c(ntRNA = 1, rbRNA = 1 , spatial = 0.4),  #alpha spatial别太大不然cluster就分块，0-0.5都差不多
  layer_weight = c(ntRNA = 1, rbRNA = 1, spatial = 1),#默认都为1
  ndim = 50 #注意原始输出的1：50个Uembedding，理论上最能区分各个细胞的方向在U[50]上，对应的最小特征值也是eigs[50]
)
D <- ComputeGrassmannDistance(
  U_consensus = scml_result$U,
  U_layers = scml_result$U_layers
)#先不管这个function可以跳过，这个是用于计算多个模态间的k子维度（k也即是ndim的值）在grassmann流形上的一个投影距离。。。。。（摘自SCML,2014,IEEE）
D
gc()

eigs = scml_result$eigenvalues
SCML_embedding <- scml_result$U

ord <- order(eigs)

eigs<- eigs[ord]#reorder
eigs[1]
rownames(SCML_embedding) <- colnames(states_seu)
SCML_embedding <- SCML_embedding[, ord]

colnames(SCML_embedding) <- paste0(
  "SCML_",
  seq_len(ncol(SCML_embedding))
)
L_scml <- scml_result$L_modified 

DefaultAssay(states_seu) <- "RNA"

states_seu[["scml"]] <- CreateDimReducObject(
  embeddings = SCML_embedding,
  key = "SCML_",
  assay = "RNA"
)
states_seu <- FindNeighbors(
  states_seu,
  reduction = "scml",
  dims = 1:50,
  k.param = 50, annoy.metric = "cosine", #也可以选择欧氏距离，可以前后统一。我觉得scml后的Uembedding维度更高，方向是否相似可能（这种可能极其主观色彩）稍微好一点？但其实差不多理论上
  graph.name = c(
    "scml_nn",
    "scml_snn"
  )
)
states_seu <- FindClusters(
  states_seu, graph.name = 'scml_snn',
  resolution = 0.5
)

states_seu <- RunUMAP(
  states_seu,
  reduction = "scml",reduction.name = "scml.umap",
  dims = 1:50, spread = 2.0, min.dist = 0.1
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



# 4. Imputation -----------------------------------------------------------

rb_raw_mat = GetAssayData(states_seu, layer = 'counts', assay = 'rbRNA')
total_raw_mat = GetAssayData(states_seu, layer  = 'counts', assay = 'RNA')
nt_raw_mat = GetAssayData(states_seu, layer = 'counts', assay = 'ntRNA')

total_norm <- GetAssayData(states_seu, layer = 'data',assay = 'RNA')
rb_norm <- GetAssayData(states_seu, layer = 'data',assay = 'rbRNA')
nt_norm <- GetAssayData(states_seu, layer = 'data',assay = 'ntRNA')

#W <- build_knn_graph(SCML_embedding, k=50)

W = states_seu@graphs$scml_snn

spec <- ComputeGraphSpectrum(
  W,
  n_eigs=300 #这里的选择跟后面的滤波可以适配起来
)

L = spec$L
lambda <- spec$lambda
V <- spec$V

ord <- order(lambda) # reorder
lambda <- lambda[ord]
V <- V[, ord]
lambda[1] #应该为一个极其接近0的数字

rb_result <- RunGeneWiseGraphFilter(
  expr = rb_norm,
  V = V,
  lambda = lambda,
  L = L,
  beta_min = 2,
  beta_max = 2
)
rb_filtered <- rb_result$filtered
G_rb = rb_result$G
colnames(rb_filtered) <- colnames(rb_norm)
all.equal(rownames(rb_filtered),rownames(rb_norm))


nt_result <- RunGeneWiseGraphFilter(
  expr = nt_norm,
  V = V,
  lambda = lambda,
  L = L,
  beta_min = 2,
  beta_max = 2
)
nt_filtered <- nt_result$filtered
colnames(nt_filtered) <- colnames(nt_norm)
all.equal(rownames(nt_filtered),rownames(nt_norm))

total_result <- RunGeneWiseGraphFilter(
  expr = total_norm,
  V = V,
  lambda = lambda,
  L = L,
  beta_min = 2,
  beta_max = 2
)
total_filtered <- total_result$filtered
colnames(total_filtered) <- colnames(total_norm)

##rescale
rb_q99 <- QuantileRescaleGene(
  expr_raw = rb_norm,
  expr_denoise = rb_filtered,
  q = 0.99
)

rb_filtered_rescale <- rb_q99$expr

nt_q99 <- QuantileRescaleGene(
  expr_raw = nt_norm,
  expr_denoise = nt_filtered,
  q = 0.99
)

nt_filtered_rescale <- nt_q99$expr

total_q99 <- QuantileRescaleGene(
  expr_raw = total_norm,
  expr_denoise = total_filtered,
  q = 0.99
)

total_filtered_rescale <- total_q99$expr

states_seu[["rb_filtered"]] <-
  CreateAssayObject(
    counts =  rb_filtered_rescale
  )
states_seu[["nt_filtered"]] <-
  CreateAssayObject(
    counts = nt_filtered_rescale
  )
states_seu[["total_filtered"]] <-
  CreateAssayObject(
    counts = total_filtered_rescale
  )

te_raw <- rb_raw_mat / (total_raw_mat + 1e-6)
te_filtered <- rb_filtered / (rb_filtered + nt_filtered+ 1e-6) #注意，这里TE没有加权，对于TE 的robust的建模以及可视化仍然需要考量，但nt，rb，total是可以直接用的

states_seu[["te_raw"]] <-
  CreateAssayObject(
    counts = te_raw
  )
states_seu[["te_filtered"]] <-
  CreateAssayObject(
    counts = te_filtered
  )

PlotSpatialGene2 <- function(
    seu,
    gene,
    assay,
    layer = "data",
    vmax,
    title = NULL,
    point_size = 0.10,
    n_colors = 128
) {
  #绘图函数，color需要改
  expr <- GetAssayData(
    seu,
    assay = assay,
    layer = layer
  )[gene, ]
  
  expr_plot <- pmin(expr, vmax)
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
  reds <-  colorRampPalette(
  c(
    "#453781",
    "#31688e",
    "#21908d",
    "#6ece58",
    "#b8de29"
  ),
  space = "Lab"
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


genes_use <- c("Lamp5","Rgs4")
#"Enpp2","Rarres2","Ptgds",,"Chgb","C1ql2","Cux2","Cplx3","Ppp1r1b","Clu"

plot_total <- list()

for(gene in genes_use){
  
  vmax_raw  <- GetVmax(total_norm, gene)
  vmax_wave <- GetVmax(GetAssayData(states_seu, assay = 'total_filtered',layer = 'counts'), gene)
  
  plot_total[[length(plot_total)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "RNA",layer = 'data',
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
  plot_annotation(title = "Total RNA")
p_total

plot_rb <- list()

for(gene in genes_use){
  
  vmax_raw  <- GetVmax(rb_norm, gene)
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

plot_te <- list()

for(gene in genes_use){
  
  vmax_raw  <- GetVmax(te_raw, gene)
  vmax_wave <- GetVmax(te_filtered, gene)
  
  plot_te[[length(plot_te)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "te_raw",layer = 'counts',
      vmax = vmax_raw,
      title = paste0(gene, "\nTE Raw")
    )
  
  plot_te[[length(plot_te)+1]] <-
    PlotSpatialGene2(
      states_seu,
      gene = gene,
      assay = "te_filtered", layer = 'counts', 
      vmax = vmax_wave,
      title = paste0(gene, "\nTE Filtered")
    )
}

p_te <- wrap_plots(plot_te, ncol = 2) +
  plot_annotation(title = "Translation Efficiency")
p_te

#插补后dotplot
genes_use <- c(
  "Ly6c1",
  "Ly6e",
  "Ptgds",
  "Rgs5",
  "Rarres2",
  "Cllc6",
  "Folr1",
  "Enpp2",
  "Ttr",
  "Ctss",
  "C1qc",
  "C1qa",
  "Hexb",
  "Cacng4",
  "Olig2",
  "Sox10",
  "Pdgfra",
  "Ptprz1",
  "Trf",
  "Cnp",
  "Fth1",
  "Mobp",
  "Mbp",
  "Plp1",
  "Clu",
  "Mt2",
  "Atp1b2",
  "Gja1",
  "Aldoc",
  "Resp18",
  "Hap1",
  "Cadps2",
  "Gng8",
  "Tac2",
  "Tubb5",
  "Ntng1",
  "Prkcd",
  "Sparc",
  "Sparcl1",
  "Npy",
  "Sst",
  "Pvalb",
  "Gad2",
  "Gad1",
  "Hpca",
  "Nell2",
  "Mapk1",
  "Ppp3r1",
  "Nrgn"
)
DefaultAssay(states_seu) <- "total_filtered"
p <- DotPlot(
  states_seu,
  features = genes_use,
  group.by = "mousebrain_label2",
  assay = "total_filtered"
) +
  scale_color_gradientn(
    colours = c(
      "#3B4CC0",
      "white",
      "#B40426"
    )
  ) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )

p







