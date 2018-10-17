copuladd_compute <- function(data, subset, weights = NULL, na.action, 
                             output.type = c('edgelist'), 
                             bootstrap.number = 1000, 
                             bootstrap.alpha = 0.05,
                             bootstrap.ratio = 0.7)
{
  # Computes copula directional dependence between each pair of variables
  # 
  # Input
  # 
  #   data          - [nxp] matrix with real valued entries
  #   output.type   - Default is 'edgelist'
  #   bootstrap.size - If >0, then subsamples of size [n x bootstrap.size]
  #                   are resampled.
  #
  # Output
  #   The dd values between variables are returned in the specified format
  #   
  #   1)'edgelist'  - a matrix with more than two columns: 
  #                   "node_from", "node_to", "dd_diff", 
  #                   "dd_diff_lb", "dd_diff_ub",
  #                   "vtou_rho2", "utov_rho2"
  #
  # < CDD for FMRI >
  #   
  #   Copyright (C) 2018 Namgil Lee & Jong-Min Kim
  
  
  require(copula)
  require(gcmr)
  #require(boot)
  #if (bootstrap.number>0) {
  #  require(gmodels)
  #}
  
  #-------------------------------- Inputs ------------------------------------#
  #      'data' is a [nxp] matrix
  #----------------------------------------------------------------------------#
  if (is.null(n <- nrow(data))) 
    stop("'data' must be a matrix")
  if (n == 0L) 
    stop("0 cases")
  p <- ncol(data)
  if (p == 0L) {
    return(data.frame(node_from = numeric(), node_to = numeric(), 
                      dd_from2to = numeric(), dd_to2from = numeric(),
                      dd_diff = numeric(), dd_diff_lb = numeric(), 
                      dd_diff_ub = numeric()))
  }
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  
  #-------------------------------- Output ------------------------------------#
  cl <- match.call()
  
  ## Node names (Not used)
  # node_names <- colnames(data)
  # if(is.null(node_names))
  #   node_names = as.character(1:p)
  
  #=============================================================================#
  #                    Transform original data to U(0,1)                        #
  #=============================================================================#
  
  Emp.index <- matrix(rep(0,n*p),n,p)
  for(i in 1:p)
    Emp.index[,i] <- empiric_df(data[,i],data[,i])
  
  #=============================================================================#
  #                    Compute Directional Dependence                           #
  #=============================================================================#
  num_edge <- p*(p-1)/2
  ret <- as.data.frame(matrix(0,num_edge,7))
  names(ret) <- c('node_from', 'node_to', 'dd_from2to', 'dd_to2from',
                  'dd_diff', 'dd_diff_lb', 'dd_diff_ub')
  boot_diff_all = vector('list',0) #list of all boot differences(rho2); length=pairs of (idu,idv)
  
  cur_row <- 0
  if (p>=2) {
    for (idu in 1:(p-1)) {
      for (idv in (idu+1):p) {
        cur_row <- cur_row + 1
        
        rslt <- copuladd_compute_2d(Emp.index[,c(idu,idv)], out_diff = FALSE)
        utov_rho2 <- rslt$utov_rho2
        vtou_rho2 <- rslt$vtou_rho2
        
        #---------- Record dd --------#
        if ( utov_rho2 > vtou_rho2 ) { #Direction is: u -> v 
          ret$node_from[cur_row] = idu
          ret$node_to[cur_row] = idv
          ret$dd_from2to[cur_row] = utov_rho2
          ret$dd_to2from[cur_row] = vtou_rho2
          ret$dd_diff[cur_row] = utov_rho2 - vtou_rho2
            
        } else { #Direction is: v -> u
          ret$node_from[cur_row] = idv
          ret$node_to[cur_row] = idu
          ret$dd_from2to[cur_row] = vtou_rho2
          ret$dd_to2from[cur_row] = utov_rho2
          ret$dd_diff[cur_row] = vtou_rho2 - utov_rho2
        }
        #------------------------------#
        
        #========================================================================#
        #                  Bootstrap Confidence Interval                         #
        #========================================================================#
        if (bootstrap.number > 0) {

          ##########OLD_BOOTSTRAP_ALGORITHM_FOR_CONFIDENCE_INTERVAL##########
          # boot.out <- boot(Emp.index[,c(idu,idv)], copuladd_compute_2d, R=bootstrap.number)
          # myci = boot.ci(boot.out, conf = 1-bootstrap.alpha, type = "basic")
          # 
          # #------- Record CI ----------#
          # if ( utov_rho2 > vtou_rho2 ) { #Direction is: u -> v 
          #   ret$dd_diff_lb[cur_row] = myci$basic[1,4]
          #   ret$dd_diff_ub[cur_row] = myci$basic[1,5]
          # 
          # } else {
          #   ret$dd_diff_lb[cur_row] = -myci$basic[1,5]
          #   ret$dd_diff_ub[cur_row] = -myci$basic[1,4]
          # }
          # #----------------------------#
          ############################NEW ALGORITHM##########################
          nsubset = floor(bootstrap.ratio * n)
          if (nsubset >= n) {
            print('In:copuladd_compute():bootstrap.ratio is too large (>=1)')
          }
          boot_diff <- rep(0,bootstrap.number)
          for (k in 1:bootstrap.number) {
            ## Select subset from data
            ## Compute differences in all bootstraps
            boot_diff[k] <- copuladd_compute_2d(Emp.index[sample(n,nsubset),c(idu,idv)], out_diff = TRUE)
          }
          #------- Record CI ----------#
          #myci = ci(boot_diff, confidence = 1-bootstrap.alpha)
          nnomiss = sum(!is.na(boot_diff))
          est_mn <- mean(boot_diff, na.rm = TRUE)  #mean
          est_se <- sd(boot_diff, na.rm = TRUE)/sqrt(nnomiss) #standard error
          ci.low  <- est_mn + qt(bootstrap.alpha/2, nnomiss-1) * est_se
          ci.high  <- est_mn - qt(bootstrap.alpha/2, nnomiss-1) * est_se
          if ( utov_rho2 > vtou_rho2 ) { #Direction is: u -> v
            ret$dd_diff_lb[cur_row] = ci.low
            ret$dd_diff_ub[cur_row] = ci.high

          } else {
            ret$dd_diff_lb[cur_row] = -ci.high
            ret$dd_diff_ub[cur_row] = -ci.low
          }
          #------- Record all boot differences------#
          boot_diff_all[[cur_row]] <- boot_diff
          #################################
          
          
        }#if(bootstrap.number)


      }#for(idv)
    }#for(idu)
  }#if(p>=2)
    
  #----------------------- Return value ----------------------# 
  class(ret) <- c(class(ret), 'copuladd')
  attr(ret,'call') <- cl
  attr(ret,'boot_diff') <- boot_diff_all
  ret
}