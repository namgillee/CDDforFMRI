print_copuladd_network_centrality <- function(x_all, idx_x_sub, 
                                  num_nodes = NULL, node_names = NULL, 
                                  type = 'degree', overlay = TRUE)
{ 
  ## Compute centrality measures and print on screen. 
  ## We use bi-directional network.
  ##
  ## Reference: draw_copuladd_network_bidirectional.R
  
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
  

  #*****************************************************************************#
  #             Make igraph graph object                                        #
  #*****************************************************************************#
  # See, e.g., igraph reference manual for 'make_graph'
  # 1) Define a graph with fixed number of vertices
  # 2) Set node labels
  if (!isTRUE(overlay)) {
    #Do not overlay
    #Only draw a subgraph
    
    #So, we overwrite the whole graph, x_all, by  x_all[idx_x_sub,]
    #Then, idx_x_sub   is not valid, so set it TRUE
    x_all <- x_all[idx_x_sub,,drop=FALSE]
    idx_x_sub <- rep(TRUE, nrow(x_all))
  } else {
    #Overlay networks
    #Draw whole graph and its subgraph
    
    #Use (most of) the whole graph, x_all. 
    #Do not overwrite (but slightly modify) idx_x_sub
    idx_x_thresholded <-  (x_all$dd_from2to >= 0.2)  # Edge THRESHOLD for OVERLAYED NETWORK 
    x_all <- x_all[idx_x_thresholded,,drop=FALSE]
    idx_x_sub <- idx_x_sub[idx_x_thresholded]
  }
  
  
  ## Set Edges: bidirectional
  # The edges are a stack of one-side-edges and other-side-edges
  # We will concatenate the cdd matrix x_all with the other-side, 
  # and then build edges
  # Step-1)  We already know that x_all[,3] >= 0.2
  # Step-2)  Lets concatenate x_all with  x_all[,4] >= 0.2
  idx_otherdir = (x_all[,4] >= 0.2)
  x_otherdir = x_all[idx_otherdir,c(2,1,4,3,5:7)]
  x_all  =  rbind(as.matrix(x_all),    as.matrix(x_otherdir)) 
  x_all  =  as.data.frame(x_all, row.names = 1:nrow(x_all))
  idx_x_sub  =  c(idx_x_sub,  rep(FALSE,sum(idx_otherdir))) #edges to color by red
  idx_edge_curved = c(idx_otherdir, rep(TRUE,sum(idx_otherdir))) #edges to make curved
  
  
  myedge <- t(x_all[,1:2])   #Draw whole graph
  num_edges = ncol(myedge)
  mygraph <- make_graph(edges = myedge, n = num_nodes, directed = TRUE)
  V(mygraph)$name <- node_names

  
  # 1. Degree distribution
  print('*Degree*')
  print(degree(mygraph))
  print('*Betweenness')
  print(betweenness(mygraph))
  

  # ###### REDUNDANT OPTIONS WHICH IS FOR DRAWING ######
  # ## Set options(width, color, label) for node edges
  # E(mygraph)$width = rep(1, num_edges) ##default is 1
  # if (length(idx_x_sub)>0)
  #   E(mygraph)$width[idx_x_sub] = 1 + x_all$dd_diff[idx_x_sub]*40
  # E(mygraph)$color = rep('darkgrey', num_edges) ##default is darkgrey
  # if (length(idx_x_sub)>0)
  #   E(mygraph)$color[idx_x_sub] = 'red'
  # E(mygraph)$lty = rep(2, num_edges) ## default is dotted
  # if (length(idx_x_sub)>0)
  #   E(mygraph)$lty[idx_x_sub] = 1 
  
  # if (edgetext) { 
  #   E(mygraph)$label = rep("", num_edges)  ##default is ""
  #   if (length(idx_x_sub)>0)
  #     E(mygraph)$label[idx_x_sub] = round(x_all$dd_diff[idx_x_sub], 3)
  #   E(mygraph)$label.cex = 2.0
  # }

  # ## Node layout
  # if (identical(tolower(graphlayout),'circle')) {
  #   coords <- layout_in_circle(mygraph, order = 1:num_nodes)
  # } else {
  #   #random 
  #   coords <- layout_with_fr(mygraph)
  # }


  # ######### Draw a network ##########  
  # pdf(file_pdf, width=12,height=12)
  # plot(mygraph, layout = coords, vertex.size=30, vertex.label.cex=2, edge.arrow.size=1, 
  #      edge.curved=(idx_edge_curved)*0.2)
  # dev.off()
  # 
  # png(file_png, width=12*2.5,height=12*2.5,units='cm', res=300)
  # plot(mygraph, layout = coords, vertex.size=30, vertex.label.cex=2, edge.arrow.size=1, 
  #      edge.curved=(idx_edge_curved)*0.2)
  # dev.off()

  return(0)
  
}

