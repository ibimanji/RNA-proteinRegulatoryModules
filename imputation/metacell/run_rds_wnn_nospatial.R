#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(igraph)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

parse_args <- function(args) {
  cfg <- list(
    input = "states_with_plaque_info.rds",
    output_rds = "states_with_plaque_info_supercell2_like_metacells.rds",
    membership_csv = "states_with_plaque_info_supercell2_like_membership.csv",
    label = "states_nn_alg1_label3",
    split_by = "type",
    split_values = "all",
    total_assay = "totalRNA",
    rb_assay = "rbRNA",
    te_assay = "TE",
    total_reduction = "totalRNA_pca",
    rb_reduction = "rbRNA_pca",
    total_dims = "1:30",
    rb_dims = "1:30",
    gamma = 30,
    k_knn = 30,
    kith = NA,
    kernel = TRUE,
    nfeatures = 3000,
    min_cells_one_metacell = 150,
    seed = 12345
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      cfg$input <- key
      i <- i + 1
      next
    }

    name <- sub("^--", "", key)
    if (!name %in% names(cfg)) stop("Unknown argument: ", key)
    if (i == length(args)) stop("Missing value for argument: ", key)
    cfg[[name]] <- args[[i + 1]]
    i <- i + 2
  }

  cfg$gamma <- as.numeric(cfg$gamma)
  cfg$k_knn <- as.integer(cfg$k_knn)
  cfg$kith <- as.integer(cfg$kith)
  if (is.na(cfg$kith)) cfg$kith <- NULL
  cfg$kernel <- as.logical(cfg$kernel)
  cfg$nfeatures <- as.integer(cfg$nfeatures)
  cfg$min_cells_one_metacell <- as.integer(cfg$min_cells_one_metacell)
  cfg$seed <- as.integer(cfg$seed)
  cfg
}

parse_dims <- function(x) {
  x <- gsub("\\s+", "", x)
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    parts <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    return(seq(parts[[1]], parts[[2]]))
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_split_values <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return(NULL)
  x <- trimws(as.character(x))
  if (tolower(x) %in% c("all", "*")) return(NULL)
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  values <- values[nzchar(values)]
  if (length(values) == 0) return(NULL)
  unique(values)
}

get_assay_matrix <- function(object, assay, slot_or_layer) {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    return(SeuratObject::GetAssayData(object, assay = assay, layer = slot_or_layer))
  }
  SeuratObject::GetAssayData(object, assay = assay, slot = slot_or_layer)
}

as_seurat_matrix <- function(mat) {
  if (inherits(mat, "dgCMatrix")) return(mat)
  if (inherits(mat, "sparseMatrix")) return(methods::as(mat, "dgCMatrix"))
  if (inherits(mat, "Matrix")) return(methods::as(mat, "dgCMatrix"))
  mat
}

create_assay_with_data <- function(counts, data) {
  counts <- as_seurat_matrix(counts)
  data <- as_seurat_matrix(data)

  if (exists("CreateAssay5Object", envir = asNamespace("SeuratObject"), inherits = FALSE)) {
    return(SeuratObject::CreateAssay5Object(counts = counts, data = data))
  }

  assay <- SeuratObject::CreateAssayObject(counts = counts)
  assay <- SeuratObject::SetAssayData(assay, slot = "data", new.data = data)
  assay
}

replace_nonfinite_with_zero <- function(mat, label = "matrix") {
  if (inherits(mat, "sparseMatrix")) {
    bad <- which(!is.finite(mat@x))
    if (length(bad) > 0) {
      message("Replacing ", length(bad), " non-finite values with 0 in ", label, ".")
      mat@x[bad] <- 0
      mat <- Matrix::drop0(mat)
    }
    return(mat)
  }

  bad <- which(!is.finite(mat), arr.ind = TRUE)
  if (nrow(bad) > 0) {
    message("Replacing ", nrow(bad), " non-finite values with 0 in ", label, ".")
    mat[bad] <- 0
  }
  mat
}

get_data_matrix <- function(object, assay, slot_or_layer) {
  get_assay_matrix(object, assay, slot_or_layer)
}

set_assay_data_matrix <- function(assay, slot_or_layer, data) {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    return(SeuratObject::SetAssayData(assay, layer = slot_or_layer, new.data = data))
  }
  SeuratObject::SetAssayData(assay, slot = slot_or_layer, new.data = data)
}

align_matrices <- function(..., names_for_error = NULL) {
  mats <- list(...)
  common_features <- Reduce(intersect, lapply(mats, rownames))
  common_cells <- Reduce(intersect, lapply(mats, colnames))

  if (length(common_features) == 0) {
    stop("No shared features across matrices", names_for_error %||% ".")
  }
  if (length(common_cells) == 0) {
    stop("No shared cells across matrices", names_for_error %||% ".")
  }

  lapply(mats, function(mat) mat[common_features, common_cells, drop = FALSE])
}

build_analysis_seurat_from_rds <- function(cfg) {
  message("Reading Seurat RDS: ", cfg$input)
  object <- readRDS(cfg$input)
  if (!inherits(object, "Seurat")) stop("Input RDS is not a Seurat object.")

  required_assays <- c(cfg$total_assay, cfg$rb_assay, cfg$te_assay)
  missing_assays <- setdiff(required_assays, names(object@assays))
  if (length(missing_assays) > 0) {
    stop("Missing assays in input object: ", paste(missing_assays, collapse = ", "))
  }
  if (!cfg$label %in% colnames(object@meta.data)) {
    stop("Label column not found in meta.data: ", cfg$label)
  }
  if (!is.null(cfg$split_by) && !cfg$split_by %in% colnames(object@meta.data)) {
    stop("split_by column not found in meta.data: ", cfg$split_by)
  }

  message("Extracting totalRNA counts/data, rbRNA counts, and TE counts...")
  total_counts <- get_assay_matrix(object, cfg$total_assay, "counts")
  total_data <- get_assay_matrix(object, cfg$total_assay, "data")
  rb_counts <- get_assay_matrix(object, cfg$rb_assay, "counts")
  te_counts <- get_assay_matrix(object, cfg$te_assay, "counts")

  aligned <- align_matrices(
    total_counts,
    total_data,
    rb_counts,
    te_counts,
    names_for_error = " totalRNA counts/data, rbRNA counts, and TE counts."
  )
  total_counts <- aligned[[1]]
  total_data <- aligned[[2]]
  rb_counts <- aligned[[3]]
  te_counts <- aligned[[4]]

  te_counts <- replace_nonfinite_with_zero(te_counts, label = paste0(cfg$te_assay, " counts"))
  total_data <- replace_nonfinite_with_zero(total_data, label = paste0(cfg$total_assay, " data"))

  message("Computing normalized rbRNA as TE counts * totalRNA data...")
  rb_data <- te_counts * total_data
  rb_data <- replace_nonfinite_with_zero(rb_data, label = paste0(cfg$rb_assay, " data"))
  total_counts <- as_seurat_matrix(total_counts)
  total_data <- as_seurat_matrix(total_data)
  rb_counts <- as_seurat_matrix(rb_counts)
  rb_data <- as_seurat_matrix(rb_data)

  meta <- object@meta.data[colnames(total_counts), , drop = FALSE]
  seurat <- CreateSeuratObject(
    counts = total_counts,
    assay = cfg$total_assay,
    meta.data = meta
  )
  seurat[[cfg$total_assay]] <- create_assay_with_data(total_counts, total_data)
  seurat[[cfg$rb_assay]] <- create_assay_with_data(rb_counts, rb_data)
  DefaultAssay(seurat) <- cfg$total_assay
  seurat
}

subset_seurat_cells <- function(seurat, cells) {
  cells <- intersect(cells, colnames(seurat))
  if (length(cells) == 0) stop("No cells left after subsetting.")
  subset(seurat, cells = cells)
}

top_variable_features <- function(mat, nfeatures) {
  nfeatures <- min(nfeatures, nrow(mat))
  if (nfeatures >= nrow(mat)) return(rownames(mat))

  if (inherits(mat, "Matrix")) {
    mu <- Matrix::rowMeans(mat)
    mu2 <- Matrix::rowMeans(mat^2)
    vars <- as.numeric(mu2 - mu^2)
  } else {
    vars <- apply(mat, 1, stats::var)
  }
  vars[!is.finite(vars)] <- 0
  keep <- which(vars > 0)
  if (length(keep) == 0) {
    stop("No positive-variance features found for assay preprocessing.")
  }
  keep <- keep[order(vars[keep], decreasing = TRUE)]
  rownames(mat)[keep[seq_len(min(nfeatures, length(keep)))]]
}

sanitize_reduction <- function(seurat, reduction) {
  emb <- Seurat::Embeddings(seurat, reduction = reduction)
  bad <- which(!is.finite(emb), arr.ind = TRUE)
  if (nrow(bad) > 0) {
    message("Replacing ", nrow(bad), " non-finite values with 0 in reduction ", reduction, ".")
    emb[bad] <- 0
    seurat@reductions[[reduction]]@cell.embeddings <- emb
  }
  seurat
}

preprocess_assay <- function(seurat, assay, reduction_name, dims, nfeatures, seed) {
  DefaultAssay(seurat) <- assay
  data_mat <- get_data_matrix(seurat, assay = assay, slot_or_layer = "data")
  VariableFeatures(seurat) <- top_variable_features(data_mat, nfeatures)
  seurat <- ScaleData(seurat, assay = assay, features = VariableFeatures(seurat), verbose = FALSE)
  seurat <- RunPCA(
    seurat,
    assay = assay,
    features = VariableFeatures(seurat),
    reduction.name = reduction_name,
    npcs = max(dims),
    seed.use = seed,
    verbose = FALSE
  )
  seurat <- sanitize_reduction(seurat, reduction_name)
  seurat
}

make_complete_graph <- function(cells, weight = 1) {
  graph <- igraph::make_full_graph(length(cells), directed = FALSE, loops = FALSE)
  igraph::V(graph)$name <- cells
  if (igraph::ecount(graph) > 0) igraph::E(graph)$weight <- weight
  graph
}

print_graph_diagnostics <- function(graph, label = "graph") {
  weights <- igraph::E(graph)$weight
  if (is.null(weights)) weights <- rep(1, igraph::ecount(graph))
  strength <- igraph::strength(graph, weights = weights)
  zero_strength <- sum(!is.finite(strength) | strength <= 0)

  message(
    "[graph diagnostics] ", label,
    " | vertices=", igraph::vcount(graph),
    " edges=", igraph::ecount(graph),
    " zero_strength=", zero_strength
  )
  if (length(strength) > 0) {
    q <- stats::quantile(strength, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    message(
      "[graph diagnostics] ", label,
      " | strength min/q1/median/q3/max = ",
      paste(signif(as.numeric(q), 4), collapse = " / ")
    )
  }
}

compute_unimodal_knn <- function(
    seurat,
    assay,
    reduction,
    dims,
    k_knn = 30,
    kith = NULL,
    kernel = TRUE,
    label = NULL,
    subset_label = NULL,
    verbose = FALSE) {
  if (!is.null(label) && !is.null(subset_label)) {
    labels <- seurat[[label]][, 1]
    keep <- rownames(seurat@meta.data)[!is.na(labels) & labels == subset_label]
    seurat <- seurat[, keep]
  }

  if (ncol(seurat) <= 1) {
    graph <- igraph::make_empty_graph(n = ncol(seurat), directed = FALSE)
    igraph::V(graph)$name <- colnames(seurat)
    return(graph)
  }

  k_knn <- min(k_knn, ncol(seurat) - 1L)
  DefaultAssay(seurat) <- assay
  seurat <- Seurat::FindNeighbors(
    seurat,
    reduction = reduction,
    dims = dims,
    k.param = k_knn,
    return.neighbor = TRUE,
    verbose = verbose
  )

  graph_name <- paste0(assay, ".nn")
  nn <- seurat@neighbors[[graph_name]]

  if (kernel) {
    kith <- kith %||% max(2L, k_knn %/% 2L)
    kith <- min(kith, ncol(nn@nn.dist))
    j <- as.numeric(t(nn@nn.idx))
    i <- ((seq_along(j) - 1L) %/% k_knn) + 1L
    x <- as.numeric(t(nn@nn.dist))

    adj <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(ncol(seurat), ncol(seurat)))
    rownames(adj) <- colnames(seurat)
    colnames(adj) <- colnames(seurat)
    graph <- igraph::graph_from_adjacency_matrix(adj, diag = FALSE, mode = "directed", weighted = TRUE)
    graph <- igraph::as.undirected(graph, edge.attr.comb = "mean")

    adj <- igraph::as_adjacency_matrix(graph, attr = "weight")
    sigmas <- nn@nn.dist[, kith]
    sigmas[sigmas <= 0] <- stats::median(sigmas[sigmas > 0])
    sigmas[is.na(sigmas)] <- 1
    adj <- expm1(-((adj %*% Matrix::Diagonal(x = 1 / sigmas))^2)) + igraph::as_adjacency_matrix(graph, sparse = TRUE)
    adj <- adj + Matrix::t(adj)
    graph <- igraph::graph_from_adjacency_matrix(adj, diag = FALSE, mode = "undirected", weighted = TRUE)
  } else {
    graph_name <- paste0(assay, "_nn")
    adj <- seurat@graphs[[graph_name]] + Matrix::t(seurat@graphs[[graph_name]])
    graph <- igraph::graph_from_adjacency_matrix(adj, diag = FALSE, mode = "undirected", weighted = TRUE)
    igraph::E(graph)$weight <- 1
  }

  igraph::V(graph)$name <- colnames(seurat)
  graph
}

compute_multimodal_knn <- function(
    seurat,
    assays,
    reductions,
    dims,
    k_knn = 30,
    kernel = TRUE,
    label = NULL,
    subset_label = NULL,
    verbose = FALSE) {
  if (!is.null(label) && !is.null(subset_label)) {
    labels <- seurat[[label]][, 1]
    keep <- rownames(seurat@meta.data)[!is.na(labels) & labels == subset_label]
    seurat <- seurat[, keep]
  }

  if (ncol(seurat) <= 1) {
    graph <- igraph::make_empty_graph(n = ncol(seurat), directed = FALSE)
    igraph::V(graph)$name <- colnames(seurat)
    return(graph)
  }
  if (ncol(seurat) <= 3) {
    return(make_complete_graph(colnames(seurat), weight = 1))
  }

  k_knn <- min(k_knn, ncol(seurat) - 1L)
  knn_range <- min(200L, ncol(seurat))
  knn_graph_name <- "localwnnknn"
  snn_graph_name <- "localwnnsnn"
  weighted_nn_name <- "localwnn"
  searching <- TRUE
  last_error <- NULL
  while (searching) {
    searching <- FALSE
    seurat <- tryCatch(
      Seurat::FindMultiModalNeighbors(
        seurat,
        reduction.list = reductions,
        dims.list = dims,
        k.nn = k_knn,
        knn.range = knn_range,
        knn.graph.name = knn_graph_name,
        snn.graph.name = snn_graph_name,
        weighted.nn.name = weighted_nn_name,
        verbose = verbose
      ),
      error = function(e) {
        last_error <<- e
        searching <<- TRUE
        seurat
      }
    )
    if (searching) {
      knn_range <- knn_range - 1L
      if (k_knn >= knn_range) k_knn <- knn_range - 1L
      if (k_knn < 2L || knn_range < 3L) {
        message(
          "FindMultiModalNeighbors failed for a small/unstable subset; ",
          "using a complete graph fallback. Last error: ",
          conditionMessage(last_error)
        )
        return(make_complete_graph(colnames(seurat), weight = 1))
      }
    }
  }

  if (snn_graph_name %in% names(seurat@graphs)) {
    adj <- seurat@graphs[[snn_graph_name]]
    if (length(adj@x) == 0) {
      message("WNN-SNN graph has no edges; using a complete graph fallback.")
      return(make_complete_graph(colnames(seurat), weight = 1))
    }
    adj@x[!is.finite(adj@x) | adj@x <= 0] <- 1e-16
    graph <- igraph::graph_from_adjacency_matrix(
      adj,
      diag = FALSE,
      mode = "undirected",
      weighted = TRUE
    )
    igraph::V(graph)$name <- colnames(seurat)
    return(graph)
  }

  if (kernel) {
    nn <- seurat@neighbors[[weighted_nn_name]]
    if (is.null(nn) && "weighted.nn" %in% names(seurat@neighbors)) {
      nn <- seurat@neighbors[["weighted.nn"]]
    }
    if (is.null(nn)) {
      message("Weighted nearest-neighbor object not found; using a complete graph fallback.")
      return(make_complete_graph(colnames(seurat), weight = 1))
    }
    j <- as.numeric(t(nn@nn.idx))
    i <- ((seq_along(j) - 1L) %/% k_knn) + 1L

    x <- as.numeric(t(1 - 2 * (nn@nn.dist^2)))
    keep <- is.finite(i) & is.finite(j) & is.finite(x) & !is.na(i) & !is.na(j)
    i <- i[keep]
    j <- j[keep]
    x <- x[keep]
    x[!is.finite(x) | x <= 0] <- 1e-16
    if (length(i) == 0 || length(j) == 0 || length(x) == 0) {
      message("Weighted nearest-neighbor graph has no valid edges; using a complete graph fallback.")
      return(make_complete_graph(colnames(seurat), weight = 1))
    }

    adj <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(ncol(seurat), ncol(seurat)))
    rownames(adj) <- colnames(seurat)
    colnames(adj) <- colnames(seurat)
    adj <- adj + Matrix::t(adj)
    graph <- igraph::graph_from_adjacency_matrix(adj, diag = FALSE, mode = "undirected", weighted = TRUE)
  } else {
    adj <- seurat@graphs[["wknn"]] + Matrix::t(seurat@graphs[["wknn"]])
    graph <- igraph::graph_from_adjacency_matrix(adj, diag = FALSE, mode = "undirected", weighted = TRUE)
    igraph::E(graph)$weight <- 1
  }

  igraph::V(graph)$name <- colnames(seurat)
  graph
}

union_graphs_by_weight <- function(graphs, all_cells) {
  all_cells <- as.character(all_cells)
  all_edges <- list()
  all_weights <- list()

  for (idx in seq_along(graphs)) {
    graph <- graphs[[idx]]
    if (igraph::vcount(graph) == 0 || igraph::ecount(graph) == 0) next

    edge_list <- igraph::as_edgelist(graph, names = TRUE)
    weights <- igraph::E(graph)$weight
    if (is.null(weights)) weights <- rep(1, nrow(edge_list))
    weights[!is.finite(weights) | weights <= 0] <- 1e-16

    all_edges[[length(all_edges) + 1L]] <- edge_list
    all_weights[[length(all_weights) + 1L]] <- weights
  }

  if (length(all_edges) == 0) {
    message("No edges found while merging label graphs; using complete graph fallback for all cells.")
    return(make_complete_graph(all_cells, weight = 1))
  }

  edge_df <- as.data.frame(do.call(rbind, all_edges), stringsAsFactors = FALSE)
  colnames(edge_df) <- c("from", "to")
  edge_df$weight <- unlist(all_weights, use.names = FALSE)
  edge_df <- edge_df[edge_df$from %in% all_cells & edge_df$to %in% all_cells & edge_df$from != edge_df$to, , drop = FALSE]

  edge_df$a <- pmin(edge_df$from, edge_df$to)
  edge_df$b <- pmax(edge_df$from, edge_df$to)
  edge_df <- stats::aggregate(weight ~ a + b, data = edge_df, FUN = sum)

  cell_index <- stats::setNames(seq_along(all_cells), all_cells)
  i <- unname(cell_index[edge_df$a])
  j <- unname(cell_index[edge_df$b])
  keep <- !is.na(i) & !is.na(j) & is.finite(edge_df$weight) & edge_df$weight > 0
  i <- i[keep]
  j <- j[keep]
  x <- edge_df$weight[keep]

  if (length(i) == 0) {
    message("Merged label graph has no valid edges; using complete graph fallback for all cells.")
    return(make_complete_graph(all_cells, weight = 1))
  }

  adj <- Matrix::sparseMatrix(
    i = c(i, j),
    j = c(j, i),
    x = c(x, x),
    dims = c(length(all_cells), length(all_cells)),
    dimnames = list(all_cells, all_cells)
  )
  graph <- igraph::graph_from_adjacency_matrix(
    adj,
    diag = FALSE,
    mode = "undirected",
    weighted = TRUE
  )
  igraph::V(graph)$name <- all_cells
  graph
}

build_label_aware_graph <- function(seurat, assays, reductions, dims, k_knn, kernel, label, verbose, min_cells_one_metacell = 150) {
  graphs <- build_label_aware_graphs(
    seurat = seurat,
    assays = assays,
    reductions = reductions,
    dims = dims,
    k_knn = k_knn,
    kernel = kernel,
    label = label,
    verbose = verbose,
    min_cells_one_metacell = min_cells_one_metacell
  )

  graph <- union_graphs_by_weight(graphs, colnames(seurat))
  print_graph_diagnostics(graph, label = paste0(label, "=union"))
  graph
}

build_label_aware_graphs <- function(seurat, assays, reductions, dims, k_knn, kernel, label, verbose, min_cells_one_metacell = 150) {
  labels <- seurat[[label]][, 1]
  known_labels <- na.omit(unique(labels))

  graphs <- lapply(known_labels, function(one_label) {
    n_label_cells <- sum(!is.na(labels) & labels == one_label)
    message("[label-aware] ", label, "=", one_label, " | cells=", n_label_cells)
    if (n_label_cells <= min_cells_one_metacell) {
      message(
        "[label-aware] ", label, "=", one_label,
        " has <= ", min_cells_one_metacell,
        " cells; using complete graph fallback for this label."
      )
      label_cells <- rownames(seurat@meta.data)[!is.na(labels) & labels == one_label]
      graph <- make_complete_graph(label_cells, weight = 1)
    } else {
      graph <- compute_multimodal_knn(
        seurat = seurat,
        assays = assays,
        reductions = reductions,
        dims = dims,
        k_knn = k_knn,
        kernel = kernel,
        label = label,
        subset_label = one_label,
        verbose = verbose
      )
    }
    print_graph_diagnostics(graph, label = paste0(label, "=", one_label))
    graph
  })
  names(graphs) <- paste0(label, "=", make.names(as.character(known_labels), unique = TRUE))

  if (any(is.na(labels))) {
    unknowns <- rownames(seurat@meta.data)[is.na(labels)]
    message("[label-aware] ", label, "=NA | cells=", length(unknowns))
    if (length(unknowns) <= min_cells_one_metacell) {
      message(
        "[label-aware] ", label, "=NA",
        " has <= ", min_cells_one_metacell,
        " cells; using complete graph fallback for this label."
      )
      unknown_graph <- make_complete_graph(unknowns, weight = 1)
    } else {
      unknown_seurat <- subset_seurat_cells(seurat, unknowns)
      unknown_graph <- compute_multimodal_knn(
        seurat = unknown_seurat,
        assays = assays,
        reductions = reductions,
        dims = dims,
        k_knn = k_knn,
        kernel = kernel,
        verbose = verbose
      )
    }
    print_graph_diagnostics(unknown_graph, label = paste0(label, "=unknown"))
    graphs[[length(graphs) + 1L]] <- unknown_graph
    names(graphs)[length(graphs)] <- paste0(label, "=NA")
  }

  graphs
}

supercell_assign <- function(values, membership) {
  split_values <- split(values, membership)
  vapply(split_values, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1]
  }, character(1))
}

supercell_purity <- function(values, membership) {
  split_values <- split(values, membership)
  vapply(split_values, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    max(table(x)) / length(x)
  }, numeric(1))
}

cluster_one_graph_walktrap <- function(graph, n_metacells) {
  cells <- igraph::V(graph)$name
  membership <- rep(NA_integer_, length(cells))
  names(membership) <- cells
  zero_strength <- character(0)
  walktrap <- NULL
  next_cluster <- 1L

  if (length(cells) == 0) {
    return(list(
      membership = membership,
      walktrap = walktrap,
      zero_strength = zero_strength
    ))
  }

  weights <- igraph::E(graph)$weight
  if (is.null(weights)) weights <- rep(1, igraph::ecount(graph))
  strength <- igraph::strength(graph, weights = weights)
  zero_strength <- names(strength)[!is.finite(strength) | strength <= 0]
  core_cells <- setdiff(cells, zero_strength)

  if (length(core_cells) >= 2) {
    core_graph <- igraph::induced_subgraph(graph, vids = core_cells)
    core_weights <- igraph::E(core_graph)$weight
    if (is.null(core_weights)) core_weights <- rep(1, igraph::ecount(core_graph))
    core_strength <- igraph::strength(core_graph, weights = core_weights)
    still_bad <- names(core_strength)[!is.finite(core_strength) | core_strength <= 0]
    if (length(still_bad) > 0) {
      core_cells <- setdiff(core_cells, still_bad)
      zero_strength <- c(zero_strength, still_bad)
      core_graph <- igraph::induced_subgraph(graph, vids = core_cells)
    }

    if (length(core_cells) >= 2) {
      n_core_metacells <- min(length(core_cells), max(1L, n_metacells - length(zero_strength)))
      if (n_core_metacells <= 1L) {
        core_membership <- stats::setNames(rep(1L, length(core_cells)), core_cells)
      } else {
        walktrap <- igraph::cluster_walktrap(core_graph, weights = igraph::E(core_graph)$weight)
        core_membership <- igraph::cut_at(walktrap, no = n_core_metacells)
        names(core_membership) <- igraph::V(core_graph)$name
        core_membership_names <- names(core_membership)
        core_membership <- as.integer(factor(core_membership, levels = sort(unique(core_membership))))
        names(core_membership) <- core_membership_names
      }
      membership[names(core_membership)] <- core_membership
      next_cluster <- max(core_membership) + 1L
    }
  }

  if (length(core_cells) == 1) {
    membership[core_cells] <- next_cluster
    next_cluster <- next_cluster + 1L
  }

  if (length(zero_strength) > 0) {
    for (cell in zero_strength) {
      if (is.na(membership[[cell]])) {
        membership[[cell]] <- next_cluster
        next_cluster <- next_cluster + 1L
      }
    }
  }

  if (any(is.na(membership))) {
    for (cell in names(membership)[is.na(membership)]) {
      membership[[cell]] <- next_cluster
      next_cluster <- next_cluster + 1L
    }
  }

  membership <- as.integer(factor(membership, levels = sort(unique(membership))))
  names(membership) <- cells

  list(
    membership = membership,
    walktrap = walktrap,
    zero_strength = zero_strength
  )
}

cluster_graphs_walktrap <- function(graphs, gamma) {
  membership <- integer(0)
  walktraps <- vector("list", length(graphs))
  names(walktraps) <- names(graphs)
  zero_strength_cells <- vector("list", length(graphs))
  names(zero_strength_cells) <- names(graphs)
  offset <- 0L
  graph_names <- names(graphs)

  for (idx in seq_along(graphs)) {
    graph <- graphs[[idx]]
    graph_name <- if (!is.null(graph_names) && length(graph_names) >= idx && nzchar(graph_names[[idx]])) {
      graph_names[[idx]]
    } else {
      paste0("graph_", idx)
    }
    n_cells <- igraph::vcount(graph)
    n_metacells <- max(1L, floor(n_cells / gamma))
    message(
      "Running igraph::cluster_walktrap within ", graph_name,
      " | cells=", n_cells,
      " target_metacells=", n_metacells
    )

    clustered <- cluster_one_graph_walktrap(graph, n_metacells)
    graph_membership <- clustered$membership
    graph_membership <- graph_membership + offset
    membership <- c(membership, graph_membership)
    if (length(membership) > 0) offset <- max(membership)
    walktraps[[idx]] <- clustered$walktrap
    zero_strength_cells[[idx]] <- clustered$zero_strength

    if (length(clustered$zero_strength) > 0) {
      message(
        "Found ", length(clustered$zero_strength),
        " zero-strength cells in ", graph_name,
        "; assigning them to singleton metacells."
      )
    }
  }

  membership <- membership[!duplicated(names(membership))]
  membership_names <- names(membership)
  membership <- as.integer(factor(membership, levels = sort(unique(membership))))
  names(membership) <- membership_names

  list(
    membership = membership,
    walktraps = walktraps,
    zero_strength_cells = zero_strength_cells
  )
}

metacell_expression <- function(seurat, assays, group_by, layer = "counts", method = "aggregate") {
  groups <- seurat[[group_by]][, 1]
  groups <- factor(groups, levels = sort(unique(groups)))
  category <- Matrix::sparse.model.matrix(~ 0 + groups)
  rownames(category) <- colnames(seurat)
  colnames(category) <- paste0("Metacell_", seq_len(ncol(category)))

  col_sums <- Matrix::colSums(category)
  if (method == "average") {
    category <- category %*% Matrix::Diagonal(x = 1 / col_sums)
  }

  result <- list()
  for (assay in assays) {
    mat <- get_data_matrix(seurat, assay = assay, slot_or_layer = layer)
    out <- mat %*% category
    colnames(out) <- paste0("Metacell_", seq_len(ncol(out)))
    result[[assay]] <- as_seurat_matrix(out)
  }
  result
}

scimplify_for_seurat_like <- function(
    seurat,
    gamma,
    assays,
    reductions,
    dims,
    k_knn = 30,
    kith = NULL,
    kernel = TRUE,
    label = NULL,
    min_cells_one_metacell = 150,
    verbose = FALSE) {
  if (!is.null(label)) {
    message("Using label-aware graph construction with column: ", label)
    label_graphs <- build_label_aware_graphs(
      seurat = seurat,
      assays = assays,
      reductions = reductions,
      dims = dims,
      k_knn = k_knn,
      kernel = kernel,
      label = label,
      verbose = verbose,
      min_cells_one_metacell = min_cells_one_metacell
    )
    graph <- label_graphs
  } else {
    graph <- compute_multimodal_knn(seurat, assays, reductions, dims, k_knn, kernel, verbose = verbose)
  }

  if (!is.null(label)) {
    clustered <- cluster_graphs_walktrap(label_graphs, gamma)
    membership <- clustered$membership
    membership <- membership[colnames(seurat)]
    missing_cells <- colnames(seurat)[is.na(membership)]
    if (length(missing_cells) > 0) {
      next_cluster <- max(membership, na.rm = TRUE) + 1L
      for (cell in missing_cells) {
        membership[[cell]] <- next_cluster
        next_cluster <- next_cluster + 1L
      }
    }
    walktraps <- clustered$walktraps
    zero_strength <- clustered$zero_strength_cells
  } else {
    message("Running igraph::cluster_walktrap...")
    n_metacells <- max(2L, floor(ncol(seurat) / gamma))
    clustered <- cluster_one_graph_walktrap(graph, n_metacells)
    membership <- clustered$membership[colnames(seurat)]
    walktraps <- list(all = clustered$walktrap)
    zero_strength <- clustered$zero_strength

    if (length(zero_strength) > 0) {
      message(
        "Found ", length(zero_strength),
        " zero-strength cells; assigning them to singleton metacells."
      )
    }
  }

  membership <- as.integer(factor(membership, levels = sort(unique(membership))))
  names(membership) <- colnames(seurat)

  group_col <- paste0("metacell_g", gamma)
  seurat[[group_col]] <- membership

  counts_list <- metacell_expression(seurat, assays = assays, group_by = group_col, layer = "counts", method = "aggregate")
  data_list <- metacell_expression(seurat, assays = assays, group_by = group_col, layer = "data", method = "average")

  first_assay <- assays[[1]]
  seurat_mc <- CreateSeuratObject(
    counts = as_seurat_matrix(counts_list[[first_assay]]),
    assay = first_assay
  )
  seurat_mc[[first_assay]] <- create_assay_with_data(
    counts = counts_list[[first_assay]],
    data = data_list[[first_assay]]
  )
  for (assay in assays[-1]) {
    seurat_mc[[assay]] <- create_assay_with_data(
      counts = counts_list[[assay]],
      data = data_list[[assay]]
    )
  }

  categorical <- vapply(seurat@meta.data, function(x) is.character(x) || is.factor(x), logical(1))
  for (field in names(categorical)[categorical]) {
    seurat_mc[[field]] <- supercell_assign(seurat[[field]][, 1], paste0("Metacell_", membership))
    seurat_mc[[paste0(field, "_purity")]] <- supercell_purity(seurat[[field]][, 1], paste0("Metacell_", membership))
  }

  seurat_mc$size <- as.numeric(table(membership))
  seurat_mc@misc$membership <- membership
  seurat_mc@misc$metacells_hierarchy <- walktraps
  seurat_mc@misc$walktrap_clusters <- membership
  seurat_mc@misc$zero_strength_cells <- zero_strength
  seurat_mc@misc$gamma <- gamma
  seurat_mc@misc$graph <- graph
  seurat_mc
}

merge_metacell_objects <- function(metacell_objects, sample_names, assays) {
  if (length(metacell_objects) == 1) return(metacell_objects[[1]])

  renamed <- vector("list", length(metacell_objects))
  for (i in seq_along(metacell_objects)) {
    obj <- metacell_objects[[i]]
    prefix <- make.names(sample_names[[i]])
    obj <- RenameCells(obj, new.names = paste(prefix, colnames(obj), sep = "_"))
    obj$sample_type <- sample_names[[i]]
    renamed[[i]] <- obj
  }

  merged <- Reduce(function(x, y) merge(x, y), renamed)
  merged@misc$split_by <- "type"
  merged@misc$per_sample_gamma <- lapply(metacell_objects, function(x) x@misc$gamma)
  merged@misc$per_sample_n_metacells <- stats::setNames(
    vapply(metacell_objects, ncol, numeric(1)),
    sample_names
  )
  merged
}

run_one_sample <- function(seurat, cfg, sample_name = NULL) {
  total_dims <- parse_dims(cfg$total_dims)
  rb_dims <- parse_dims(cfg$rb_dims)

  if (!is.null(sample_name)) {
    message("========== Sample: ", sample_name, " ==========")
  }
  message("Cells in run: ", ncol(seurat))
  if (!is.null(cfg$label) && cfg$label %in% colnames(seurat@meta.data)) {
    label_counts <- sort(table(seurat[[cfg$label]][, 1]), decreasing = TRUE)
    message("[label counts] ", cfg$label, ":")
    print(label_counts)
  }

  if (!is.null(cfg$label) && cfg$label %in% colnames(seurat@meta.data) &&
      ncol(seurat) <= cfg$min_cells_one_metacell) {
    message(
      "Cells <= min_cells_one_metacell (", cfg$min_cells_one_metacell,
      "); assigning all cells in this run to one metacell."
    )
    membership <- stats::setNames(rep(1L, ncol(seurat)), colnames(seurat))
    group_col <- paste0("metacell_g", cfg$gamma)
    seurat[[group_col]] <- membership
    counts_list <- metacell_expression(seurat, assays = c(cfg$total_assay, cfg$rb_assay), group_by = group_col, layer = "counts", method = "aggregate")
    data_list <- metacell_expression(seurat, assays = c(cfg$total_assay, cfg$rb_assay), group_by = group_col, layer = "data", method = "average")
    seurat_mc <- CreateSeuratObject(counts = as_seurat_matrix(counts_list[[cfg$total_assay]]), assay = cfg$total_assay)
    seurat_mc[[cfg$total_assay]] <- create_assay_with_data(counts_list[[cfg$total_assay]], data_list[[cfg$total_assay]])
    seurat_mc[[cfg$rb_assay]] <- create_assay_with_data(counts_list[[cfg$rb_assay]], data_list[[cfg$rb_assay]])
    seurat_mc$size <- ncol(seurat)
    seurat_mc@misc$membership <- membership
    seurat_mc@misc$gamma <- cfg$gamma
    if (!is.null(sample_name)) {
      seurat_mc$sample_type <- sample_name
      seurat_mc@misc$sample_type <- sample_name
    }
    return(seurat_mc)
  }

  message("Preprocessing assay: ", cfg$total_assay)
  seurat <- preprocess_assay(seurat, cfg$total_assay, cfg$total_reduction, total_dims, cfg$nfeatures, cfg$seed)

  message("Preprocessing assay: ", cfg$rb_assay)
  seurat <- preprocess_assay(seurat, cfg$rb_assay, cfg$rb_reduction, rb_dims, cfg$nfeatures, cfg$seed)

  seurat_mc <- scimplify_for_seurat_like(
    seurat = seurat,
    gamma = cfg$gamma,
    assays = c(cfg$total_assay, cfg$rb_assay),
    reductions = list(cfg$total_reduction, cfg$rb_reduction),
    dims = list(total_dims, rb_dims),
    k_knn = cfg$k_knn,
    kith = cfg$kith,
    kernel = cfg$kernel,
    label = cfg$label,
    min_cells_one_metacell = cfg$min_cells_one_metacell,
    verbose = TRUE
  )
  if (!is.null(sample_name)) {
    seurat_mc$sample_type <- sample_name
    seurat_mc@misc$sample_type <- sample_name
  }
  seurat_mc
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  set.seed(cfg$seed)

  seurat <- build_analysis_seurat_from_rds(cfg)
  if (!is.null(cfg$split_by)) {
    split_values <- seurat[[cfg$split_by]][, 1]
    available_sample_names <- sort(unique(as.character(split_values[!is.na(split_values)])))
    if (length(available_sample_names) == 0) {
      stop("split_by column has no non-NA values: ", cfg$split_by)
    }
    requested_sample_names <- parse_split_values(cfg$split_values)
    if (is.null(requested_sample_names)) {
      sample_names <- available_sample_names
      message("Requested split values: all")
    } else {
      missing_sample_names <- setdiff(requested_sample_names, available_sample_names)
      if (length(missing_sample_names) > 0) {
        stop(
          "Requested split_values not found in meta.data$", cfg$split_by, ": ",
          paste(missing_sample_names, collapse = ", "),
          ". Available values: ",
          paste(available_sample_names, collapse = ", ")
        )
      }
      sample_names <- requested_sample_names
      message("Requested split values: ", paste(sample_names, collapse = ", "))
    }

    message("Running independently by meta.data$", cfg$split_by, ": ", paste(sample_names, collapse = ", "))
    metacell_objects <- lapply(sample_names, function(sample_name) {
      cells <- rownames(seurat@meta.data)[as.character(split_values) == sample_name]
      sample_seurat <- subset_seurat_cells(seurat, cells)
      run_one_sample(sample_seurat, cfg, sample_name = sample_name)
    })
    seurat_mc <- merge_metacell_objects(
      metacell_objects = metacell_objects,
      sample_names = sample_names,
      assays = c(cfg$total_assay, cfg$rb_assay)
    )
    seurat_mc@misc$per_sample_membership <- stats::setNames(
      lapply(metacell_objects, function(x) x@misc$membership),
      sample_names
    )
  } else {
    seurat_mc <- run_one_sample(seurat, cfg)
  }

  saveRDS(seurat_mc, cfg$output_rds)
  if (!is.null(cfg$split_by)) {
    membership_df <- do.call(rbind, lapply(names(seurat_mc@misc$per_sample_membership), function(sample_name) {
      membership <- seurat_mc@misc$per_sample_membership[[sample_name]]
      data.frame(
        sample_type = sample_name,
        cell = names(membership),
        metacell = paste(make.names(sample_name), paste0("Metacell_", as.integer(membership)), sep = "_"),
        membership = as.integer(membership),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    membership <- seurat_mc@misc$membership
    membership_df <- data.frame(
      cell = names(membership),
      metacell = paste0("Metacell_", as.integer(membership)),
      membership = as.integer(membership),
      stringsAsFactors = FALSE
    )
  }
  utils::write.csv(membership_df, cfg$membership_csv, row.names = FALSE)

  message("Done.")
  message("Saved Seurat metacell object: ", cfg$output_rds)
  message("Saved membership CSV: ", cfg$membership_csv)
  message("Input cells: ", ncol(seurat))
  message("Metacells: ", ncol(seurat_mc))
}

main()
