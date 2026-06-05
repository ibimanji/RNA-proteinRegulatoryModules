#seurat v5. based metacell construction demo[need refining]
getwd()
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
  library(data.table)
})
states_seu = readRDS("/media/zenglab/data/jinpu/0529statesmouse14month_ADWT/14mWT_refined.rds") #only a demo here since the annotation should be re-refined
DimPlot(states_seu, reduction = "umap", group.by = "celltype") 

#here we introduce supervised and unsupervised manner, semi-supervised manner which is bypassed here requires 
#(1)you are not sure for specific celltype annotation, this is ok in several situations
#(2)a purity test for supervised manner, like that you hide several labels for test.

############################
####1. supervised manner:###
############################
all_celltypes <- unique(states_seu$celltype)
global_metacell_labels <- rep(NA, ncol(states_seu))
names(global_metacell_labels) <- colnames(states_seu)
for (ct in all_celltypes) {
  cat(paste0("\n==================== celltype here is : ", ct, " ====================\n"))
  sub_seu <- subset(states_seu, subset = celltype == ct)
  n_cells <- ncol(sub_seu)
  cat(sprintf("   --> 当前亚群包含单细胞数量: %d\n", n_cells))
  if (n_cells <= 150) {
    cat(
      "   [提示] 细胞数不足150，全部细胞合并为单个 Metacell。\n"
    )
    formatted_ids <- paste0(
      ct,
      "_MC_001"
    ) #here is useful for CHO-PEP subtype which contains only 130s cells
    
    global_metacell_labels[
      colnames(sub_seu)
    ] <- formatted_ids
    cat(
      sprintf(
        "   [统计] 成功生成 %d 个 Metacell\n",
        1
      )
    )
    next
  }
  ########################################################
  ## STEP1 WNN-KNN
  ########################################################
  cat("   --> [STEP 1] 运行 k.nn=30 WNN 局部近邻构图...\n") #all the celltype shares one k may not be the best option here
  sub_seu <- FindMultiModalNeighbors(
    object = sub_seu,
    reduction.list = list(
      "rna.pca", #can be replaced by ntrna.pca, but in this seuratobj(refined states_seu), this is not given
      "rbrna.pca",
      "spatial.coords" #spatial.coords, ah, indeed can be integrated into reduction layer, but can it be dealt more ... finely?
    ),
    dims.list = list(
      1:30,
      1:30,
      1:2
    ),
    knn.graph.name = "local.wnn.knn",
    snn.graph.name = "local.wnn.snn",
    weighted.nn.name = "local.wnn",
    k.nn = 30
  )
  cat("   --> [STEP 2] 提取WNN-SNN图...\n")
  snn_mat <- sub_seu@graphs$local.wnn.snn
  g_rw <- igraph::graph_from_adjacency_matrix(
    snn_mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  ########################################################
  ## STEP3 Walktrap
  ########################################################
  cat("   --> [STEP 3] Walktrap clustering...\n")
  walktrap_tree <- igraph::cluster_walktrap(
    g_rw,
    weights = igraph::E(g_rw)$weight,
    steps = 4
  )
  ########################################################
  ## STEP4 gamma压缩
  ########################################################
  gamma_param <- 25 
  target_no <- floor(
    n_cells / gamma_param
  )#total cell number in this section, like if a defined gamma is given, metacell number is given
  if (target_no >= 2) {
    
    metacell_labels <- igraph::cut_at(
      walktrap_tree,
      no = target_no
    )
  } else {
    cat(
      "   [提示] 亚群细胞总数较少，采用Walktrap自适应切分...\n"
    )
    metacell_labels <- walktrap_tree$membership
  }
  formatted_ids <- paste0(
    ct,
    "_MC_",
    sprintf(
      "%03d",
      metacell_labels
    )
  )
  global_metacell_labels[
    colnames(sub_seu)
  ] <- formatted_ids
  
  mod_score <- igraph::modularity(
    g_rw,
    metacell_labels,
    weights = igraph::E(g_rw)$weight
  )#这个modularity: The modularity of a graph with respect to some division (or vertex types) measures how good the division is, or how separated are the different vertex types from each other. 
#like the argued max para in Leiden and Louvain and other cluster  
  cat(
    sprintf(
      "   [统计] 成功生成 %d 个 Metacell | Modularity = %.3f\n",
      length(unique(formatted_ids)),
      mod_score
    )
  )
}

###then is a function for metacell construction glance:
states_seu$metacell_id <- global_metacell_labels
rna_mat <- GetAssayData(states_seu, assay = "RNA", layer = "counts")
rb_mat <- GetAssayData(states_seu, assay = "rbRNA", layer = "counts")
mc_ids <- na.omit(unique(states_seu$metacell_id))
n_metacell <- length(mc_ids)
cat("Number of metacells:", n_metacell, "\n")
mc_size <- table(states_seu$metacell_id)
mc_size_df <- data.frame(
  metacell_id = names(mc_size),
  n_cells = as.integer(mc_size)
)#information for metacell glance
head(mc_size_df)
ggplot(
  mc_size_df,
  aes(n_cells)
) +
  geom_histogram(
    bins = 50
  ) +
  theme_bw() +
  labs(
    x = "Cells per metacell",
    y = "Count"
  )#distribution for mc_cellnumber ,1A

#metacell TE calculation
mc_factor <- factor(
  states_seu$metacell_id
)
mc_model <- sparse.model.matrix(
  ~0 + mc_factor
)
colnames(mc_model) <- levels(mc_factor)
mc_rna <- rna_mat %*% mc_model # [gene x cell] X [cell x metacell] 
mc_rb <- rb_mat %*% mc_model
dim(mc_rna)
dim(mc_rb)
all.equal(colnames(mc_rna),colnames(mc_rb))
all.equal(rownames(mc_rna),rownames(mc_rb))
mc_te <- mc_rb / (mc_rna+1e-06)

#metacell UMAP glanced:
library(viridis)
umap_coords <- as.data.frame(
  Embeddings(states_seu, reduction = "umap")
)
meta_df <- states_seu@meta.data
plot_data <- cbind(
  umap_coords,
  meta_df
)
plot_data_mc <- plot_data[
  !is.na(plot_data$metacell_id),
]#single cell level UMAP
mc_centroids <- plot_data_mc %>%
  group_by(
    metacell_id,
    celltype
  ) %>%
  summarise(
    UMAP_1 = mean(UMAP_1),
    UMAP_2 = mean(UMAP_2),
    Size   = n(),
    .groups = "drop"
  )#mc_UMAP坐标为改集合内所有

p_supercell <- ggplot() +
  geom_point(
    data = plot_data,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = celltype
    ),
    size = 0.4,
    alpha = 0.3,
    stroke = 0
  ) +
  geom_point(
    data = mc_centroids,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = celltype,
      size = Size
    ),
    alpha = 0.95,
    shape = 16
  ) +
  
  scale_color_viridis_d(
    option = "turbo",
    name = "Cell Type"
  ) +
  
  scale_size_continuous(
    range = c(0.8, 4.0),
    name = "Cells per Metacell"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    legend.position = "right"
  ) +
  
  labs(
    title = "Metacell Manifold Projection",
    subtitle = paste0(
      "Metacells = ",
      nrow(mc_centroids),
      " | Cells covered = ",
      round(
        100 *
          sum(!is.na(states_seu$metacell_id)) /
          ncol(states_seu),
        2
      ),
      "%"
    ),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )
print(p_supercell)




##############################
####2. unsupervised manner:###
##############################
#in this section we sealed a function, and added several glance parts for metacells
BuildGlobalMetacell <- function(
    states_seu,
    reduction.list = list(
      "rna.pca",
      "rbrna.pca",
      "spatial.coords"
    ),
    dims.list = list(
      1:30,
      1:30,
      1:2
    ),
    k.nn = 50,
    gamma_param = 35,
    walktrap_step = 4,
    umap_reduction = "umap",
    celltype_col = "celltype"
){
  cat("\n==================== 全局构建 Metacell ====================\n")
  n_cells <- ncol(states_seu)
  
  ############################################################
  ## STEP1 WNN
  ############################################################
  
  cat("   --> [STEP 1] 构建全局WNN图...\n")
  
  states_seu <- FindMultiModalNeighbors(
    object = states_seu,
    
    reduction.list = reduction.list,
    
    dims.list = dims.list,
    
    weighted.nn.name = "global.wnn",
    
    knn.graph.name = "global.wnn.knn",
    
    snn.graph.name = "global.wnn.snn",
    
    k.nn = k.nn
  )
  
  ############################################################
  ## STEP2 SNN
  ############################################################
  
  cat("   --> [STEP 2] 提取WNN-SNN图...\n")
  
  snn_mat <- states_seu@graphs$global.wnn.snn
  
  g_rw <- graph_from_adjacency_matrix(
    snn_mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  ############################################################
  ## STEP3 Walktrap
  ############################################################
  
  cat("   --> [STEP 3] Walktrap clustering...\n")
  
  walktrap_tree <- cluster_walktrap(
    g_rw,
    weights = E(g_rw)$weight,
    steps = walktrap_step
  )
  
  ############################################################
  ## STEP4 Gamma compression
  ############################################################
  
  target_no <- floor(
    n_cells / gamma_param
  )
  
  target_no <- max(
    target_no,
    2
  )
  
  cat(
    sprintf(
      "   --> Target metacells = %d\n",
      target_no
    )
  )
  
  metacell_labels <- cut_at(
    walktrap_tree,
    no = target_no
  )
  
  states_seu$metacell_id <- paste0(
    "MC_",
    sprintf(
      "%04d",
      metacell_labels
    )
  )
  
  ############################################################
  ## statistics
  ############################################################
  
  mod_score <- modularity(
    g_rw,
    metacell_labels,
    weights = E(g_rw)$weight
  )  #modularity score calculation
  
  mc_size <- table(
    states_seu$metacell_id
  ) 
  
  size_cv <- sd(mc_size) /
    mean(mc_size)
  
  ############################################################
  ## Aggregate matrix
  ############################################################
  
  mc_counts_list <- AggregateExpression(
    states_seu,
    group.by = "metacell_id",
    assays = "RNA",
    slot = "counts",
    return.seurat = FALSE
  )
  
  mc_counts <- mc_counts_list$RNA
  
  ############################################################
  ## sparsity
  ############################################################
  
  sc_counts <- GetAssayData(
    states_seu,
    assay = "RNA",
    layer = "counts"
  )
  
  sc_sparsity <-
    
    sum(sc_counts == 0) /
    
    (nrow(sc_counts) *
       ncol(sc_counts))
  
  mc_sparsity <-
    
    sum(mc_counts == 0) /
    
    (nrow(mc_counts) *
       ncol(mc_counts))
  
  ############################################################
  ## UMAP dataframe
  ############################################################
  
  umap_df <- as.data.frame(
    Embeddings(
      states_seu,
      umap_reduction
    )
  )
  
  meta_df <- states_seu@meta.data
  
  plot_df <- cbind(
    umap_df,
    meta_df
  )
  
  ############################################################
  ## purity
  ############################################################
  
  mc_vote <- plot_df %>%
    group_by(
      metacell_id,
      .data[[celltype_col]]
    ) %>%
    summarise(
      n = n(),
      .groups = "drop"
    )
  
  mc_info <- mc_vote %>%
    group_by(
      metacell_id
    ) %>%
    summarise(
      
      dominant_celltype =
        .data[[celltype_col]][
          which.max(n)
        ],
      
      purity =
        max(n) / sum(n),
      
      size =
        sum(n),
      
      .groups = "drop"
    )
  
  weighted_purity <-
    
    sum(
      mc_info$purity *
        mc_info$size
    ) /
    
    sum(
      mc_info$size
    )
  
  ############################################################
  ## purity histogram
  ############################################################
  
  p_purity <- ggplot(
    mc_info,
    aes(purity)
  ) +
    
    geom_histogram(
      bins = 30
    ) +
    
    geom_vline(
      xintercept = weighted_purity,
      linetype = 2
    ) +
    
    theme_classic() +
    
    labs(
      title = paste0(
        "Purity Distribution\n",
        "Weighted Purity = ",
        round(
          weighted_purity,
          3
        )
      )
    )
  
  ############################################################
  ## centroid
  ############################################################
  
  mc_centroid <- plot_df %>%
    group_by(
      metacell_id
    ) %>%
    summarise(
      
      UMAP_1 =
        mean(UMAP_1),
      
      UMAP_2 =
        mean(UMAP_2),
      
      .groups = "drop"
    )
  
  mc_centroid <- left_join(
    mc_centroid,
    mc_info,
    by = "metacell_id"
  )
  
  ############################################################
  ## metacell UMAP
  ############################################################
  
  p_umap <- ggplot() +
    
    geom_point(
      data = plot_df,
      aes(
        UMAP_1,
        UMAP_2
      ),
      color = "grey85",
      size = 0.15,
      alpha = 0.25
    ) +
    
    geom_point(
      data = mc_centroid,
      
      aes(
        UMAP_1,
        UMAP_2,
        
        color = dominant_celltype,
        
        size = size,
        
        alpha = purity
      )
    ) +
    
    scale_color_viridis_d(
      option = "turbo"
    ) +
    
    scale_alpha_continuous(
      range = c(
        0.3,
        1
      )
    ) +
    
    theme_classic() +
    
    labs(
      title = paste0(
        "Global Metacell UMAP\n",
        "k=",
        k.nn,
        " gamma=",
        gamma_param,
        " step=",
        walktrap_step,
        "\nWeighted Purity=",
        round(
          weighted_purity,
          3
        )
      )
    )
  
  ############################################################
  ## summary
  ############################################################
  
  summary_dt <- data.frame(
    
    k.nn = k.nn,
    
    gamma = gamma_param,
    
    walktrap_step = walktrap_step,
    
    n_cells = n_cells,
    
    n_metacells =
      length(
        unique(
          states_seu$metacell_id
        )
      ),
    
    mean_size =
      mean(mc_size),
    
    median_size =
      median(mc_size),
    
    size_cv =
      size_cv,
    
    modularity =
      mod_score,
    
    weighted_purity =
      weighted_purity,
    
    mean_purity =
      mean(mc_info$purity),
    
    median_purity =
      median(mc_info$purity),
    
    sd_purity =
      sd(mc_info$purity),
    
    sc_sparsity =
      sc_sparsity,
    
    mc_sparsity =
      mc_sparsity
  )
  
  cat(
    sprintf(
      "\n[统计] Metacells = %d\n",
      summary_dt$n_metacells
    )
  )
  
  cat(
    sprintf(
      "[统计] Modularity = %.3f\n",
      mod_score
    )
  )
  
  cat(
    sprintf(
      "[统计] Weighted Purity = %.3f\n",
      weighted_purity
    )
  )
  
  return(
    list(
      
      seu = states_seu,
      
      metacell_counts =
        mc_counts,
      
      summary =
        summary_dt,
      
      purity =
        mc_info,
      
      purity_plot =
        p_purity,
      
      metacell_plot =
        p_umap
      
    )
  )
}
res <- BuildGlobalMetacell(
  states_seu,
  
  k.nn = 50,
  
  gamma_param = 35,
  
  walktrap_step = 4
)

res$summary

print(
  res$metacell_plot
)
print(res$purity_plot)
