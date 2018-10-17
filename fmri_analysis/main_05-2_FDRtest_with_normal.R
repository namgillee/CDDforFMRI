## main_05-2 Perform the multiple hypothesis test by the false discovery rate (FDR)
## procedure of Strimmer (2008) on the Copula DD values to determine 
## significant connectivity. And draw the connectivity network of brain regions.
#
# Input: 
#    The copula DD values ('resu_copula_dd') saved in each file 
#    in the folder "resu_main_03_03/" is a table consisting of
#    rows(==edges) and columns(==features) where
#        > names(resu_copuladd)
#        [1] "node_from"  "node_to"    "dd_from2to" "dd_to2from" "dd_diff"    "dd_diff_lb" "dd_diff_ub"
# Output:
#    For each edge, we will append additional features:
#        FDR, fdr, prob.fdr
#    So the final table consists of the columns 
#        > names(resu_locfdr)
#        [1] "node_from"  "node_to"    "dd_from2to" "dd_to2from" "dd_diff"    "dd_diff_lb" "dd_diff_ub"
#            "FDR"   "fdr"   "prob.fdr"
#
# 2018.09.27. Namgil Lee & Jong-Min Kim


rm(list=ls())

## Load library and functions
library(GeneNet)
source('../functions/myFDRTOOL.R')
source('../functions/nullmodel.R')
source('../functions/ecdf.pval.R')
library(igraph)
source('../functions/draw_copuladd_network.R')
source('../functions/draw_copuladd_network_bidirectional.R')


## Set directories for reading and writing
resu_ts_dir = './resu_main_03_01/'
cdddata_dir <- './resu_main_03_03/'
out_dir <- './resu_main_05-2-normal_fdr_network/'


## Draw lfdr histogram and fitting curves
plot.locfdr = TRUE


## Iterate over all subjects   ###SUBJ 128 IS MISSING###
subj_start = 1
subj_end = 155
sink(file = file(paste0(out_dir, "f_output", subj_start, ".txt")))

for (ids in setdiff(subj_start:subj_end,128) ) {
  strsubj <- sprintf('%03d',ids)
  print(paste('** Subject', strsubj, '**'))
  
  
  ## Load CDD result 
  load(file = paste0(cdddata_dir, 'resu_main_', strsubj, '_copuladd.Rdata'))
  #resu_copuladd
  #> names(resu_copuladd)
  #[1] "node_from"  "node_to"    "dd_from2to" "dd_to2from" "dd_diff"    "dd_diff_lb" "dd_diff_ub"
  
  
  ## Build edges with CDD correlation, sqrt(rho2)
  ## Additionally, take negative values: -sqrt(rho2)
  ## to use FDRTOOL
  edge_list_u2v_pos <- cbind(sqrt(resu_copuladd$dd_from2to), resu_copuladd$node_from, resu_copuladd$node_to)
  edge_list_v2u_pos <- cbind(sqrt(resu_copuladd$dd_to2from), resu_copuladd$node_to, resu_copuladd$node_from)
  edge_list_u2v_neg <- cbind(-sqrt(resu_copuladd$dd_from2to), resu_copuladd$node_from, resu_copuladd$node_to)
  edge_list_v2u_neg <- cbind(-sqrt(resu_copuladd$dd_to2from), resu_copuladd$node_to, resu_copuladd$node_from)
  
  edge_list <- rbind(edge_list_u2v_pos, edge_list_u2v_neg, edge_list_v2u_pos, edge_list_v2u_neg)
  colnames(edge_list) <- c('rho', 'node_from', 'node_to')
  
  # Sort 
  sort_idx <- order(abs(edge_list[, 1]), decreasing = TRUE)
  edge_list <- edge_list[sort_idx, ]
  colnames(edge_list) <- c('rho', 'node_from', 'node_to')  

  
  ## Multiple testing with locfdr
  #   require('fdrtool')
  if (any(edge_list[,'rho'] > 1) || any(edge_list[,'rho'] < -1))
    stop("Data out of range: input correlations must be in [-1; 1]")
  #   out <-fdrtool(z.transform(pc),"correlation", plot=plot.locfdr)   
  if (plot.locfdr) tiff(paste0(out_dir,'hist_lfdr_', strsubj, '.tiff'))
  out <- myFDRTOOL(z.transform(edge_list[,'rho']),
                  "normal",
                  plot=plot.locfdr, 
                  verbose = FALSE)
  if (plot.locfdr) dev.off()
  # > names(out)
  # pval  qval  lfdr  statistic  param
  # NOTE: kappa ====  scale.param <- $param[1, 5]    : parameter of null model
  #       lambda ==== eta0 <- $param[1, 3]           : A weight to null model


  
  ## Save results   
  pval_rev           <- rep(0,length(out$pval))
  pval_rev[sort_idx] <- out$pval
  qval_rev           <- rep(0,length(out$qval))
  qval_rev[sort_idx] <- out$qval
  lfdr_rev           <- rep(0,length(out$lfdr))
  lfdr_rev[sort_idx] <- out$lfdr
  
  n = nrow(resu_copuladd)
  resu_locfdr <- cbind(resu_copuladd, pval = pval_rev[1:n], qval = qval_rev[1:n], lfdr = lfdr_rev[1:n], lfdr_to2from = lfdr_rev[(n+1):(2*n)])

  save(resu_locfdr, file=paste0(out_dir, "resu_main_locfdr_", strsubj, ".Rdata"))
  # > names(resu_locfdr)
  #[1] "node_from"  "node_to"    "dd_from2to" "dd_to2from" "dd_diff"    "dd_diff_lb" "dd_diff_ub" "pval"      
  #[9] "qval"       "lfdr" 
  
  
  
  
  ################# Plot Network #################
  ## Load time series data (for supplementary information) ##
  load(paste0(resu_ts_dir, 'resu_main_', strsubj, '_time_series.Rdata'))
  idxcol_nonnull <- (tsdat_roi_name != "") ##REMOVE MISSING DATA; SEE, "main_03_01...R"
  network_node_name <- c(paste("R.",LETTERS[1:7],sep=""), paste("L.",LETTERS[7:1],sep=""))
  ## tsdat_roi, tsdat_roi_name
  
  #**********************************************************#
  #             Select significant DD                        #
  #**********************************************************#
  resu_copuladd.sig <- resu_locfdr[resu_locfdr$lfdr<0.2,,drop=FALSE]
  
  print(paste('** Number of significant edges:', nrow(resu_copuladd.sig)))
  
  
  #**********************************************************#
  #             Draw significant DD and save                 #
  # NOTE:     At the moment, we defer to draw DD network    #
  #**********************************************************#
  mynodecolors <- c("green", "yellow", "cyan", "orange", "white")[c(2,2,2,3,4,1,4 , 4,1,4,3,2,2,2)]
  
  drawout <- draw_copuladd_network(x_all = resu_copuladd.sig, idx_x_sub = resu_copuladd.sig$dd_diff_lb>0,
                                   num_nodes = sum(idxcol_nonnull), node_names = network_node_name[idxcol_nonnull],
                                   set_seed = NULL, nodecol = mynodecolors[idxcol_nonnull],
                                   edgetext = TRUE, graphlayout = 'circle',
                                   file_pdf = paste0(out_dir, 'resu_net_', strsubj, '_copuladd.pdf'),
                                   file_tiff = paste0(out_dir, 'resu_net_', strsubj, '_copuladd.tiff'),
                                   overlay = TRUE, 
                                   rho2_thresh = 0, node_order = c(4:1,sum(idxcol_nonnull):5))
  
}
sink()
