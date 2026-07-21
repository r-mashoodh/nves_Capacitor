plot_grn <- function(
    influence_file,
    diffnet_file,
    outfile = NULL,
    edge_info = "weight",
    n_tfs = 20,
    edge_min = 0.1,
    cmap = "viridis",
    layout = "graphopt",
    title = NULL,
    full_output = FALSE,
    outdegree_limits = NULL,
    save_plot = TRUE,
    width = 10,
    height = 10,
    dpi = 300
) {
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(scales)
  
  # ---- 1. Read data ----
  influence <- read_tsv(influence_file, show_col_types = FALSE)
  diff_network <- read_tsv(diffnet_file, show_col_types = FALSE)
  
  # ---- 2. Validate edge_info ----
  if (!full_output && edge_info != "weight") {
    message("full_output is FALSE; falling back to edge_info = 'weight'")
    edge_info <- "weight"
  }
  if (!edge_info %in% colnames(diff_network)) {
    stop(paste0("Column '", edge_info, "' not found in diff_network file."))
  }
  
  # ---- 3. Select top differential TFs ----
  top_tfs <- influence %>%
    filter(G_score_scaled > 0) %>%
    arrange(influence_score) %>%
    slice_tail(n = n_tfs) %>%
    pull(factor)
  
  if (length(top_tfs) == 0) {
    warning("No differential TFs found.")
    return(NULL)
  }
  
  # ---- 4. Filter to TF-TF edges above cutoff ----
  tf_edges <- diff_network %>%
    filter(source %in% top_tfs, target %in% top_tfs) %>%
    filter(.data[[edge_info]] > edge_min)
  
  if (nrow(tf_edges) == 0) {
    warning("No edges remain after filtering.")
    return(NULL)
  }
  
  # ---- 5. Build directed graph ----
  tf_grn <- graph_from_data_frame(tf_edges, directed = TRUE)
  tf_grn <- delete_vertices(tf_grn, which(degree(tf_grn, mode = "all") == 0))
  
  if (vcount(tf_grn) == 0) {
    warning("No connected TFs remain.")
    return(NULL)
  }
  
  # ---- 6. Node attributes: weighted outdegree ----
  # ANANSE uses weighted outdegree (sum of edge_info for outgoing edges)
  weighted_out <- strength(tf_grn, mode = "out", weights = edge_attr(tf_grn, edge_info))
  V(tf_grn)$weighted_outdegree <- weighted_out
  
  # Node size: 600 + outdegree * 100 (ANANSE), rescaled for ggraph
  raw_size <- 600 + weighted_out * 100
  V(tf_grn)$node_size <- rescale(raw_size, to = c(4, 16))
  
  # ---- 7. Edge attributes ----
  E(tf_grn)$raw_weight <- edge_attr(tf_grn, edge_info)
  
  # ---- 8. Plot ----
  p <- ggraph(tf_grn, layout = layout) +
    geom_edge_arc(
      aes(width = raw_weight),
      strength = 0.1,
      arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
      end_cap = circle(5, "mm"),
      color = "grey60"
    ) +
    # Nodes colored by weighted outdegree
    geom_node_point(
      aes(size = node_size, fill = weighted_outdegree),
      shape = 21, color = "black", stroke = 0.5
    ) +
    # Labels
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      size = 3.5,
      fontface = "bold",
      family = "sans",
      bg.color = "white",
      bg.r = 0.15
    ) +
    # Node fill: continuous colorbar (ColorBrewer palettes)
    scale_fill_distiller(
      palette = cmap,
      direction = 1,
      name = "Outdegree\n(regulation\nof other TFs)",
      limits = outdegree_limits
    ) +
    scale_size_identity() +
    scale_edge_width_continuous(
      range = c(0.5, 1.5),
      name = edge_info
    ) +
    theme_graph(base_family = "sans", background = "white") +
    theme(
      legend.position = "right",
      text = element_text(family = "sans", size = 11),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11)
    ) +
    labs(
      title = title
      #subtitle = paste0("Top ", n_tfs, " TFs | edge_min = ", edge_min)
    )
  
  # ---- 9. Save ----
  if (save_plot && !is.null(outfile)) {
    if (grepl("\\.pdf$", outfile, ignore.case = TRUE)) {
      ggsave(outfile, p, width = width, height = height, dpi = dpi, device = cairo_pdf)
    } else {
      ggsave(outfile, p, width = width, height = height, dpi = dpi)
    }
    message("Saved to: ", outfile)
  }
  
  return(p)
}
