#' Plot Brain Connectivity Matrix
#'
#' This function plots a normalized connectivity matrix with network annotations.
#'
#' @param mat A square connectivity matrix.
#' @return A ggplot2 object displaying the heatmap.
#' @importFrom ggplot2 ggplot geom_raster scale_fill_gradientn theme element_blank element_text annotate
#' @importFrom reshape2 melt
#' @importFrom dplyr arrange pull
#' @importFrom tibble tibble
#' @export
plot_conn <- function(mat){
  # Find CSV file inside installed package
  csv_path <- system.file("extdata", "power264_sorted_node_information.csv", package = "locusCCA.CVRtesting")
  
  # Read CSV dynamically
  label <- read.csv(csv_path)
  original_label <- label[order(label$Original_ROI), ]
  original_label <- tibble(Node = original_label$smithRSN, label = original_label$SmithName)
  
  # Generate plot
 plot_network_heatmap1(mat / max(abs(mat)), 
                        network_information = original_label, 
                        network_label_vars = "label", 
                        network_line_thickness = 0.5)
}

# conn - P x P connectivity matrix
# network information dataframe/tibble with each row corresponding to a node
# network label vars - a single column name or vector of columns names containing the network labels
#    will be sorted outer to inner (arrange order)

# margin_multiplier increases the size of the margins where the network labels are place. This
#    can require some tuning if long network names are used. Default is 1


# include_network_labels is vector or scalar applies to all)

plot_network_heatmap1 <- function(conn,
                                  network_information = NULL,
                                  network_label_vars  = NULL,
                                  network_label_var_order = NULL,
                                  network_label_size  = 3,
                                  network_line_thickness = 1,
                                  network_line_color     = "black",
                                  include_network_labels = 1,
                                  margin_multiplier = 1.0,
                                  color_limits = c(-1, 1),
                                  diagonal_on = TRUE){
  
  # cleanup the arguments a little to make them more convenient for a 
  # loop over modules
  if (length(network_line_color) == 1){
    network_line_color <- rep(network_line_color, length(network_label_vars))
  }
  if (length(include_network_labels) == 1){
    include_network_labels <- rep(include_network_labels, length(network_label_vars))
  }
  
  # Make sure that the sorted order only contains variables that are found in the network
  # information table
  # TODO -> list structure for multiple variables?
  
  # First step is to sort by network information (if has been provided)
  if (!is.null(network_label_vars)){
    new_network_table         <- network_information
    new_network_table$tempInd <- 1:nrow(network_information)
    
    # Create factor versions of the label variables. This is to give us a lot of 
    # flexibility in how we sort the data
    # Determine the number of requested network annotations (for example, network, subnetwork)
    n_network_annotations <- length(network_label_vars)
    print(paste0("number of annotations: ", n_network_annotations))
    
    # Sort the labels in the order requested by the user\
    # TWO CASES
    # (1) Single vector for n_network_annotations == 1
    # (2) List of vectors -> needed if several layers of annotation (general case)
    if (is.list(network_label_var_order)){
      
      for (i_annotation in 1:n_network_annotations){
        
        # Check if an annotation order has been provided for this variable
        if (!is.null(network_label_var_order[[network_label_vars[i_annotation]]])){
          
          variable_ordering     <- network_label_var_order[[network_label_vars[i_annotation]]]
          network_labels        <- new_network_table[,network_label_vars[i_annotation]]
          unique_network_labels <- unique( pull(network_labels) )
          not_specified         <- unique_network_labels[(unique_network_labels %in% variable_ordering) == FALSE]
          unique_network_labels <- c(variable_ordering, not_specified)
          new_network_table[,network_label_vars[i_annotation]] = factor( pull(new_network_table[,network_label_vars[i_annotation]]), levels=unique_network_labels)
          
        }
        
      }  
      
      # Handle special case where there is only one annotation variable and user input a vector
      # instead of a list
    } else {
      
      if (n_network_annotations == 1){
        if (!is.null(network_label_var_order)){
          print("sorting - single vector")
          # This is at the level of each node
          print(network_label_vars[1])
          network_labels        <- new_network_table[,network_label_vars[1]]
          unique_network_labels <- unique( pull(network_labels) )
          not_specified         <- unique_network_labels[(unique_network_labels %in% network_label_var_order) == FALSE]
          unique_network_labels <- c(network_label_var_order, not_specified)
          new_network_table[,network_label_vars[1]] = factor( pull(new_network_table[,network_label_vars[1]]), levels=unique_network_labels)
        } 
      }
      
    }
    
    # Sort by the network label vars - TODO SPECIFY ORDERING OF NETWORKS
    sorted_network_table <- new_network_table %>% arrange( across(network_label_vars) )
    sorted_node_order <- sorted_network_table$tempInd
    
    # Convert the factors back to character variables
    sorted_network_table[,network_label_vars] = paste0( pull(sorted_network_table[,network_label_vars]))
  } else {
    sorted_node_order <- 1:nrow(conn)
  }
  
  # Get the node ordder list based on the sort above that will be used to
  # re-order conn argument
  conn <- conn[sorted_node_order, sorted_node_order]
  
  # Below is for ggplot version
  color_palette = colorRampPalette(c("blue","white", "red"))(1e3)
  
  D <- dim(conn)[1]
  
  if (diagonal_on == FALSE){
    diag(conn) <- 0.0
  }
  
  heatmap_data <- melt(conn)
  
  # Generate the plot
  FN_plot <- ggplot(heatmap_data, aes(x = Var1, y = Var2)) +
    geom_raster(aes(fill = value), interpolate=TRUE) +
    scale_fill_gradientn(name = "",
                         colours = color_palette,
                         limits = color_limits,
                         oob = scales::squish) +
    theme(legend.key.height = unit(dev.size()[2] / 10, "inches")) + 
    theme(axis.title.y=element_blank(),
          axis.text.y =element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x =element_blank(),
          axis.ticks.x = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    coord_cartesian(clip = "off") +
    theme(plot.margin=unit(c(1,1,1 * margin_multiplier,1 * margin_multiplier),"cm"))
  
  # Add network annotation if provided 
  # TODO add check that provided
  if (!is.null(network_information) && !is.null(network_label_vars) ){
    
    # Determine the number of requested network annotations (for example, network, subnetwork)
    n_network_annotations <- length(network_label_vars)
    
    for (i_annotation in 1:n_network_annotations){
      
      # This is at the level of each node
      network_labels        <- sorted_network_table[,network_label_vars[i_annotation]]
      
      network_start_indices <- as.numeric(paste0(rownames(sorted_network_table)[!duplicated(network_labels)]))
      network_line_locs     <- (network_start_indices - 1/2)[-1] # first line not needed (far LHS)
      network_labels[is.na(network_labels)] <- "Missing"
      line_data             <- data.frame(
        x = c(0,network_line_locs,D),
        xend = c(0,network_line_locs,D),
        y = rep(0, length(network_line_locs)+2),
        yend = rep(D,  length(network_line_locs)+2)
      )
      FN_plot = FN_plot + 
        geom_segment(data = line_data, aes(x=x, y=y, xend=xend, yend=yend), size = network_line_thickness, color = network_line_color[i_annotation]) +
        geom_segment(data = line_data, aes(x=y, y=x, xend=yend, yend=xend), size = network_line_thickness, color = network_line_color[i_annotation])
      
      # Add network labels
      # TODO option for which sides annotations are on
      if (include_network_labels[i_annotation] == 1){
        network_midpoints = (network_start_indices + c(network_start_indices[-1], D)) / 2
        FN_plot = FN_plot +
          annotate("text", x = network_midpoints+5, y = 0,
                   label = paste0(unique(pull(network_labels)),'  '),
                   size = network_label_size[i_annotation],
                   angle = 45, hjust=1,
                   color = network_line_color[i_annotation]) +
          annotate("text", y = network_midpoints, x = 0,  label = paste0(unique(pull(network_labels)),'  '),
                   hjust = 1,
                   size = network_label_size[i_annotation],
                   color = network_line_color[i_annotation]) 
      } # end of check that network labels were requested on the plot for this level
    }
    
  }
  
  
  return(FN_plot)
}



