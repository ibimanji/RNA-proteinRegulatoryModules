PlotSpatialGene2 <- function(
    seu,
    gene,
    assay,
    vmax,
    title = NULL,
    point_size = 1
){
  
  expr <- FetchData(
    seu,
    vars = gene,
    assay = assay
  )[,1]
  
  expr_plot <- pmin(expr, vmax)
  
  df <- data.frame(
    row  = seu$row,
    col  = seu$column,
    expr = expr_plot
  )
  
  ggplot(
    df,
    aes(
      x = col,
      y = -row,
      color = expr
    )
  ) +
    geom_point(size = point_size) +
    scale_color_gradientn(
      colours = c(
        "white",
        "#FEE8C8",
        "#FDBB84",
        "#E34A33",
        "#B30000"
      ),
      limits = c(0, vmax),
      oob = scales::squish
    ) +
    coord_fixed() +
    theme_classic(base_size = 14) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = title,
      color = gene
    )
}

############################################################
## Helper: gene-specific vmax (NO shared vmax anymore)
############################################################

GetVmax <- function(mat, gene, qmax = 0.99){
  
  quantile(
    mat[gene, ],
    probs = qmax,
    na.rm = TRUE
  )
}
