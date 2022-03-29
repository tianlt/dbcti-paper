
#monocle_plot##################
monocle_plot<-function (cds, x = 1, y = 2, color_by = "State", show_tree = TRUE, 
          show_backbone = TRUE, backbone_color = "black", markers = NULL, 
          use_color_gradient = FALSE, markers_linear = FALSE, show_cell_names = FALSE, 
          show_state_number = FALSE, cell_size = 1.5, cell_link_size = 0.75, 
          cell_name_size = 2, state_number_size = 2.9, show_branch_points = TRUE, 
          theta = 0, ...) 
{
  requireNamespace("igraph")
  gene_short_name <- NA
  sample_name <- NA
  sample_state <- pData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA
  lib_info_with_pseudo <- pData(cds)
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- reducedDimS(cds)
  }
  else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
    reduced_dim_coords <- reducedDimK(cds)
  }
  else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
    select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>% 
    mutate(sample_name = rownames(.), sample_state = rownames(.))
  dp_mst <- minSpanningTree(cds)
  if (is.null(dp_mst)) {
    stop("You must first call orderCells() before using this function")
  }
  edge_df <- dp_mst %>% igraph::as_data_frame() %>% select_(source = "from", 
                                                            target = "to") %>% left_join(ica_space_df %>% select_(source = "sample_name", 
                                                                                                                  source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                         by = "source") %>% left_join(ica_space_df %>% select_(target = "sample_name", 
                                                                                                                                               target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                                                      by = "target")
  data_df <- t(monocle::reducedDimS(cds)) %>% as.data.frame() %>% 
    select_(data_dim_1 = x, data_dim_2 = y) %>% rownames_to_column("sample_name") %>% 
    mutate(sample_state) %>% left_join(lib_info_with_pseudo %>% 
                                         rownames_to_column("sample_name"), by = "sample_name")
  return_rotation_mat <- function(theta) {
    theta <- theta/180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
           nrow = 2)
  }
  rot_mat <- return_rotation_mat(theta)
  cn1 <- c("data_dim_1", "data_dim_2")
  cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
  cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
  data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
  edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
  edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
  markers_exprs <- NULL
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% 
                              markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData), 
                                                          ])))
      colnames(markers_exprs)[1:2] <- c("feature_id", 
                                        "cell_id")
      markers_exprs <- merge(markers_exprs, markers_fData, 
                             by.x = "feature_id", by.y = "row.names")
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
      0) {
    data_df <- merge(data_df, markers_exprs, by.x = "sample_name", 
                     by.y = "cell_id")
    if (use_color_gradient) {
      if (markers_linear) {
        g <- ggplot(data = data_df, aes(x = data_dim_1, 
                                        y = data_dim_2)) + geom_point(aes(color = value), 
                                                                      size = 15, na.rm = TRUE) + scale_color_viridis(name = paste0("value"), 
                                                                                                                               ...) + facet_wrap(~feature_label)
      }
      else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, 
                                        y = data_dim_2)) + geom_point(aes(color = log10(value + 
                                                                                          0.1)), size = 15, na.rm = TRUE) + 
          scale_color_viridis(name = paste0("log10(value + 0.1)"), 
                              ...) + facet_wrap(~feature_label)
      }
    }
    else {
      if (markers_linear) {
        g <- ggplot(data = data_df, aes(x = data_dim_1, 
                                        y = data_dim_2, size = (value * 0.1))) + facet_wrap(~feature_label)
      }
      else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, 
                                        y = data_dim_2, size = log10(value + 0.1))) + 
          facet_wrap(~feature_label)
      }
    }
  }
  else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
  }
  if (show_tree) {
    g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                     yend = "target_prin_graph_dim_2"), size = 2, 
                          linetype = "solid", na.rm = TRUE, data = edge_df)
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
      0) {
    if (use_color_gradient) {
    }
    else {
      g <- g + geom_point(aes_string(color = color_by), 
                          na.rm = TRUE)
    }
  }
  else {
    if (use_color_gradient) {
    }
    else {
      g <- g + geom_point(aes_string(color = color_by), 
                          size = 4, na.rm = TRUE)#######################point size
    }
  }
  if (show_branch_points && cds@dim_reduce_type == "DDRTree") {
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- ica_space_df %>% slice(match(mst_branch_nodes, 
                                                    sample_name)) %>% mutate(branch_point_idx = seq_len(n()))
    g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
                                   y = "prin_graph_dim_2"), size = 10, na.rm = TRUE,  ########branching point size
                        branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                                y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                                     size = 4, color = "white", na.rm = TRUE, branch_point_df)
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  if (show_state_number) {
    g <- g + geom_text(aes(label = sample_state), size = state_number_size)
  }
  g <- g  + 
   theme_classic()+xlab('x')+ylab('y') + theme_classic()+xlab('x')+ylab('y') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                                                              axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),
                                                                                                                                              axis.title.y=element_blank(), legend.title=element_text(size=18, face="bold"),legend.text=element_text(size=16))
  g
}

#sco_plot################
sco_plot <- function (space, progression_group = NULL, path = NULL, contour = FALSE, 
          progression_group_palette = NULL, point_size = 2, point_alpha = 1, 
          path_size = 0.5, path_alpha = 1, contour_alpha = 0.2) 
{
  check_numeric_matrix(space, "space", finite = TRUE)
  check_numeric_matrix(path, "path", finite = TRUE, is_nullable = TRUE)
  check_logical_vector(contour, "contour", length = 1)
  check_numeric_vector(progression_group, "progression_group", 
                       is_nullable = TRUE, finite = TRUE, length = nrow(space), 
                       factor = TRUE)
  min <- min(space[, 1:2])
  max <- max(space[, 1:2])
  diff <- (max - min)/2
  space_df <- data.frame(space[, 1:2], check.rows = FALSE, 
                         check.names = FALSE, stringsAsFactors = FALSE)
  colnames(space_df) <- c("Comp1", "Comp2")
  if (!is.null(progression_group)) 
    space_df$progression_group <- progression_group
  lim <- if (contour) {
    c(min - 0.1 * diff, max + 0.1 * diff)
  }
  else {
    c(min, max)
  }
  g <- ggplot() + theme_classic() + labs(x = "Component 1", 
                                         y = "Component 2", colour = "Group", fill = "Group") + 
    xlim(min - diff, max + diff) + ylim(min - diff, max + 
                                          diff) + coord_equal(xlim = lim, ylim = lim)
  if (contour) {
    if (!is.null(progression_group) && is.numeric(progression_group)) {
      stop("If contour is TRUE, the progression group must be a factor or a character.")
    }
    aes_contour <- aes_string("Comp1", "Comp2", z = "density")
    if (!is.null(progression_group)) 
      aes_contour$fill <- quote(progression_group)
    groupings <- if (is.null(progression_group)) {
      list(group = seq_len(nrow(space_df)))
    }
    else {
      unique_groups <- unique(progression_group)
      gr <- lapply(unique_groups, function(col) which(col == 
                                                        progression_group))
      names(gr) <- unique_groups
      gr
    }
    density_df <- map_df(names(groupings), function(group_name) {
      group_ix <- groupings[[group_name]]
      kde_out <- MASS::kde2d(space_df[group_ix, 1], space_df[group_ix, 
                                                             2], lims = c(min - diff, max + diff, min - diff, 
                                                                          max + diff))
      rownames(kde_out$z) <- names(kde_out$x) <- paste0("row", 
                                                        seq_along(kde_out$x))
      colnames(kde_out$z) <- names(kde_out$y) <- paste0("col", 
                                                        seq_along(kde_out$y))
      names(dimnames(kde_out$z)) <- c("x", "y")
      kde_out$z %>% reshape2::melt(kde_out$z, varnames = c("x", 
                                                           "y"), value.name = "density") %>% as_tibble() %>% 
        transmute(progression_group = group_name, Comp1 = kde_out$x[.data$x], 
                  Comp2 = kde_out$y[.data$y], density = .data$density)
    })
    if (!is.null(progression_group) && is.factor(progression_group)) 
      density_df$progression_group <- factor(density_df$progression_group, 
                                             levels = levels(progression_group))
    g <- g + stat_contour(geom = "polygon", aes_contour, 
                          density_df, breaks = 1, alpha = contour_alpha)
  }
  aes_point <- aes_string("Comp1", "Comp2")
  if (!is.null(progression_group)) 
    aes_point$colour <- quote(progression_group)
  g <- g + geom_point(aes_point, space_df, size = 5, 
                      alpha = point_alpha)
  if (!is.null(path)) 
    g <- g + geom_path(aes_string("Comp1", "Comp2"), data.frame(path), 
                       size = 2, alpha = path_alpha) #################path thickness
  palette <- if (!is.null(progression_group_palette)) {
    progression_group_palette
  }
  else if (is.character(progression_group) || is.factor(progression_group)) {
    .default_discrete_palette(progression_group)
  }
  else if (is.numeric(progression_group)) {
    .default_continuous_palette()
  }
  if (is.character(progression_group) || is.factor(progression_group)) {
    if (!is.null(progression_group_palette)) {
      if (is.null(names(progression_group_palette)) || 
          !setequal(names(progression_group_palette), 
                    progression_group)) {
        stop("progression_group_palette must be a named vector of colours\n", 
             "where the names correspond to the unique values in progression_group")
      }
    }
    g <- g + scale_color_manual(values = palette)
    if (contour) {
      g <- g + scale_fill_manual(values = palette)
    }
  }
  else if (is.numeric(progression_group)) {
    g <- g + scale_color_gradientn(colours = palette)
  }
  g <- g+ 
    theme_classic()+xlab('x')+ylab('y') + theme_classic()+xlab('x')+ylab('y') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                      axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),
                                                                                      axis.title.y=element_blank(), legend.title=element_text(size=18, face="bold"),legend.text=element_text(size=16))
  g
}

check_numeric_matrix <- function(x, param_name, is_nullable = FALSE, finite = FALSE, sparse = FALSE) {
  if (is_nullable && is.null(x)) {
    return(invisible())
  }
  
  sparse <- sparse && dynutils::is_sparse(x)
  
  if (!sparse) {
    check <- is.matrix(x) || is.data.frame(x)
    
    j <- 1
    while (j <= ncol(x) && check) {
      check <- check && is.numeric(x[,j]) && (!finite || all(is.finite(x[,j])))
      j <- j + 1
    }
  } else {
    check <- is.numeric(x@x) && (!finite || all(is.finite(x@x)))
    
  }
  
  if (!check) {
    error <- paste0(
      sQuote(param_name),
      " must be ",
      ifelse(is_nullable, "NULL, ", ""),
      "a numeric matrix, ",
      ifelse(sparse, "a sparse numeric matrix, ", ""),
      " or a data frame containing only ",
      ifelse(finite, "finite ", ""),
      "numeric values."
    )
    stop(error)
  }
}

check_logical_vector <- function(x, param_name, is_nullable = TRUE, length = NULL) {
  if (is_nullable && is.null(x)) {
    return(invisible())
  }
  
  check <- is.logical(x)
  
  check <- check && (is.null(length) || length(x) == length)
  
  if (!check) {
    error <- paste0(
      sQuote(param_name),
      " must be a logical vector consisting of ",
      ifelse(!is.null(length), paste0(length, " "), ""),
      "logical(s)"
    )
    
    stop(error)
  }
}

check_numeric_vector <- function(x, param_name, is_nullable = TRUE, finite = FALSE, whole = FALSE, range = NULL, length = NULL, factor = FALSE) {
  if (is_nullable && is.null(x)) {
    return(invisible())
  }
  
  if (factor && is.factor(x)) {
    x <- as.numeric(x)
  }
  
  check <- is.numeric(x)
  
  check <- check && (!finite || all(is.finite(x)))
  check <- check && (!whole || all(round(x) == x))
  check <- check && (is.null(range) || all(range[[1]] <= x & x <= range[[2]]))
  check <- check && (is.null(length) || length(x) == length)
  
  if (!check) {
    error <- paste0(
      sQuote(param_name),
      " must be a numeric vector consisting of ",
      ifelse(!is.null(length), paste0(length, " "), ""),
      ifelse(finite, "finite ", ""),
      ifelse(whole, "whole ", ""),
      "number(s)",
      ifelse(!is.null(range), paste0(" within the range of [", range[[1]], ", ", range[[2]], "]"), "")
    )
    
    stop(error)
  }
}

.default_discrete_palette <- function(progression_group) {
  num_progressions <- length(levels(progression_group))
  progression_cols <-
    if (num_progressions <= 9) {
      RColorBrewer::brewer.pal(num_progressions, "Set1")
    } else {
      .gg_color_hue(num_progressions)
    }
  stats::setNames(progression_cols, levels(progression_group))
}


.default_continuous_palette <- function() {
  rev(RColorBrewer::brewer.pal(9, "RdYlBu"))
}

#tscan_plot##############
tscan_plot <- function (mclustobj, x = 1, y = 2, MSTorder = NULL, show_tree = T, 
          show_cell_names = T, cell_name_size = 3, markerexpr = NULL) 
{
  color_by = "State"
  lib_info_with_pseudo <- data.frame(State = mclustobj$clusterid, 
                                     sample_name = names(mclustobj$clusterid))
  lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
  S_matrix <- mclustobj$pcareduceres
  pca_space_df <- data.frame(S_matrix[, c(x, y)])
  colnames(pca_space_df) <- c("pca_dim_1", "pca_dim_2")
  pca_space_df$sample_name <- row.names(pca_space_df)
  edge_df <- merge(pca_space_df, lib_info_with_pseudo, by.x = "sample_name", 
                   by.y = "sample_name")
  edge_df$markerexpr <- markerexpr[edge_df$sample_name]
  if (!is.null(markerexpr)) {
    g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2, 
                                    size = markerexpr))
    g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE)
  }
  else {
    g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2))
    g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE, 
                        size = 5) #####################point size
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  if (show_tree) {
    clucenter <- mclustobj$clucenter[, c(x, y)]
    clulines <- NULL
    if (is.null(MSTorder)) {
      allsp <- shortest.paths(mclustobj$MSTtree)
      longestsp <- which(allsp == max(allsp), arr.ind = T)
      MSTorder <- get.shortest.paths(mclustobj$MSTtree, 
                                     longestsp[1, 1], longestsp[1, 2])$vpath[[1]]
    }
    for (i in 1:(length(MSTorder) - 1)) {
      clulines <- rbind(clulines, c(clucenter[MSTorder[i], 
                                              ], clucenter[MSTorder[i + 1], ]))
    }
    clulines <- data.frame(x = clulines[, 1], xend = clulines[, 
                                                              3], y = clulines[, 2], yend = clulines[, 4])
    g <- g + geom_segment(aes_string(x = "x", xend = "xend", 
                                     y = "y", yend = "yend", size = NULL), data = clulines, 
                          size = 3) #################################path thickness
    clucenter <- data.frame(x = clucenter[, 1], y = clucenter[, 
                                                              2], id = 1:nrow(clucenter))
    g <- g + geom_text(aes_string(label = "id", x = "x", 
                                  y = "y", size = NULL), data = clucenter, size = 10)
  }
  g <- g + guides(colour = guide_legend(override.aes = list(size = 5))) + 
   
     theme_classic()+xlab('x')+ylab('y') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                 axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),
                                                                                                 axis.title.y=element_blank(), legend.title=element_text(size=18, face="bold"),legend.text=element_text(size=16))
  g
}

#plot_index_tscan
plot_index_tscan <- function(object, index){
  alledges <- as.data.frame(igraph::get.edgelist(object$MSTtree), 
                            stringsAsFactors = F)
  alledges[, 1] <- as.numeric(alledges[, 1])
  alledges[, 2] <- as.numeric(alledges[, 2])
  ggdata = cbind(as.data.frame(object[["pcareduceres"]]), as.character(index))
  g <- ggplot2::ggplot() + ggplot2::geom_point(data = ggdata, aes(x = PC1, y = PC2, color = index), size = 5) + ###############point size
    
    theme_classic()+xlab('x')+ylab('y') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),
                                                axis.title.y=element_blank(), legend.title=element_text(size=18, face="bold"),legend.text=element_text(size=16))
  for (i in 1:nrow(alledges)) {
    linedata = data.frame(x = object$clucenter[c(alledges[i, 1],alledges[i, 2]), 1], y = object$clucenter[c(alledges[i, 1],alledges[i, 2]), 2])
    g <- g + ggplot2::geom_line(data = linedata, aes(x=x,y=y), size = 3)  ###############path thickness
  }
  plot(g)
}
