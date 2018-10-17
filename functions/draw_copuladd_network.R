draw_copuladd_network <- function(x_all, idx_x_sub, 
                                  num_nodes = NULL, node_names = NULL, 
                                  set_seed = NULL, nodecol = NULL,
                                  edgetext = TRUE, graphlayout = 'circle', 
                                  file_pdf = NULL, 
                                  file_tiff = NULL, 
                                  overlay = FALSE,
                                  rho2_thresh = 0.2, node_order = NULL)
{ 
  ## Draw networks and save figures 
  # INPUTS:
  #  idx_x_sub     - resu_copuladd.sig <- resu_copuladd[idx_x_sub,,drop=FALSE]; index for subgraph
  #  x_all         - resu_copuladd; it is a whole graph
  #  num_nodes     - NULL or integer
  #  node_names    - NULL or a vector of characters of the same length as num_nodes
  #  set_seed      - NULL or integer
  #  nodecol       - NULL or a vector of colors (of the same length as num_nodes, or of length 1)
  #  edgetext      - TRUE if want to write DD on each edge
  #  graphlayout   - 'circle'
  #  file_pdf = NULL or file name string
  #  file_png = NULL or file name string
  #  overlay = TRUE if the larger x_all should be drawn under x = resu_copuladd.sig
  #
  # < CDD for FMRI >
  #   
  #   Copyright (C) 2018 Namgil Lee & Jong-Min Kim
  

  ## Assume x_all is not null, x_all is a matrix or data.frame
  if (is.null(x_all)){
    warning('Input x_all is NULL')
    return(-1)
  }
  
  ## Define num_nodes and node_names for plotting. 'x_all' is also used for this.
  # Skip checking the uniqueness
  nodes = union(x_all[,1], x_all[,2]) #collect in {node_from, node_to}
  if (is.null(num_nodes)) {
    if (is.null(node_names)) {
      num_nodes = length(nodes)
      node_names = as.character(nodes)
    } else {
      num_nodes = length(node_names)
      #to verify
      if (num_nodes < length(nodes)) {
        print('Input node_names has smaller number of variables than in Input x_all. Ignore node_names')
        num_nodes = length(nodes)
        node_names = as.character(nodes)
      }
    }
  } else {
    ##Check Input num_nodes
    if (num_nodes < length(nodes)) {
      print('Input num_nodes has smaller elements than in Input x_all. Ignore num_nodes and node_names')
      num_nodes = length(nodes)
      
      #I guess the input node_names may mismatch with num_nodes
      node_names = as.character(nodes)
    } else {
    
        if (is.null(node_names)) {
          ##note that as.character(nodes) can have smaller length than the desired 'num_nodes'
          node_names <- as.character(1:num_nodes) 
          
        } else {
          ######## Both num_nodes and node_names are given #########
          if (length(node_names) < num_nodes) {
            warning('Input node_names has shorter length than desired num_nodes')
          }
        }
      
    }
  }
  
  
  if (is.null(node_order)) {
    node_order <- 1:num_nodes
  }
  
  ## Set a random seed value
  if (!is.null(set_seed)) {
    set.seed(set_seed)
  }
  
  #*****************************************************************************#
  #             Make igraph graph object                                        #
  #*****************************************************************************#
  # See, e.g., igraph reference manual for 'make_graph'
  # 1) Define a graph with fixed number of vertices
  # 2) Set node labels
  if (!isTRUE(overlay)) {
    #Do not overlay
    #Only draw a subgraph
    #So, we overwrite the whole graph, x_all, but use x_all[idx_x_sub,]
    #We overwrite idx_x_sub as well.
    x_all <- x_all[idx_x_sub,,drop=FALSE]
    idx_x_sub <- rep(TRUE, nrow(x_all))
  } else {
    #Overlay networks
    #Draw whole graph and its subgraph
    #Use most of the whole graph, x_all. 
    #Do not overwrite but slightly modify idx_x_sub
    idx_x_correlated <-  (x_all$dd_from2to >= rho2_thresh)  ###!!!<--------------- THRESHOLD for OVERLAYED NETWORK ---
    x_all <- x_all[idx_x_correlated,,drop=FALSE]
    idx_x_sub <- idx_x_sub[idx_x_correlated]
  }
  
  
  myedge <- t(x_all[,1:2])   #Draw whole graph
  num_edges = ncol(myedge)
  mygraph <- make_graph(edges = myedge, n = num_nodes, directed = TRUE)
  V(mygraph)$name <- node_names


  ## Set node colors
  # SKIP checking length of col (length(col)==num_nodes || length(col)==1)
  if (!is.null(nodecol)) {
    V(mygraph)$color = adjustcolor(nodecol, 0.3)
  }


  ## Set options(width, color, label) for node edges
  E(mygraph)$width = rep(1, num_edges) ##default is 1
  if (length(idx_x_sub)>0)
    E(mygraph)$width[idx_x_sub] = 1 + x_all$dd_diff[idx_x_sub]*40
  E(mygraph)$color = rep('darkgrey', num_edges) ##default is darkgrey
  if (length(idx_x_sub)>0)
    E(mygraph)$color[idx_x_sub] = 'red'
  E(mygraph)$lty = rep(2, num_edges) ## default is dotted
  if (length(idx_x_sub)>0)
    E(mygraph)$lty[idx_x_sub] = 1 

  if (edgetext) { 
    E(mygraph)$label = rep("", num_edges)  ##default is ""
    if (length(idx_x_sub)>0)
      E(mygraph)$label[idx_x_sub] = round(x_all$dd_diff[idx_x_sub], 3)
    E(mygraph)$label.cex = 1.5
  }

  ## Node layout
  if (identical(tolower(graphlayout),'circle')) {
    coords <- layout_in_circle(mygraph, order = node_order)
  } else {
    #random 
    coords <- layout_with_fr(mygraph)
  }


  ######### Draw a network ##########  
  pdf(file_pdf, width=10, height=10)
  plot(mygraph, layout = coords, vertex.size=25, vertex.label.cex=2.5, edge.arrow.size=1)
  dev.off()
  
  #png(file_png, width=12*2.5,height=12*2.5,units='cm', res=300)
  tiff(file_tiff, width=10*2.5,height=10*2.5,units='cm', res=300)
  plot(mygraph, layout = coords, vertex.size=25, vertex.label.cex=2.5, edge.arrow.size=1)
  dev.off()

  return(0)
  
}

