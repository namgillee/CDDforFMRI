copuladd_compute_pd <- function(dat, idx = 1:nrow(dat), out_diff = TRUE)
{
  # < CDD for FMRI >
  #   
  #   Copyright (C) 2018 Namgil Lee & Jong-Min Kim
  
  #dat is a [nxp] matrix
  #  with variables ordered as ""u", "v", "w.1", ..., "w.p-2"
  
  dat <- as.data.frame(dat[idx,,drop=FALSE])
  p <- ncol(dat)
  
  colnames(dat)[1:2] <- c("u","v")
  
  r12 <- gcmr( u~., data = dat, marginal = beta.marg(link = "logit"), cormat = arma.cormat(0, 0) )
  z12 <- exp( r12$estimate[1] + colSums(t(dat[, -1, drop=FALSE]) * r12$estimate[-c(1,p+1)]) )
  Er12 <- z12/(1 + z12)
  vtou_rho2 <- var(Er12)/var(dat$u)

  r21 <- gcmr( v~., data = dat, marginal = beta.marg(link = "logit"), cormat = arma.cormat(0, 0) )
  z21 <- exp( r21$estimate[1] + colSums(t(dat[, -2, drop=FALSE]) * r21$estimate[-c(1,p+1)]) )
  Er21 <- z21/(1 + z21)
  utov_rho2 <- var(Er21)/var(dat$v)

  if (out_diff) {
    #---------------------------------------#
    # Result output is the differnece of..  #
    #---------------------------------------#
    rslt <- utov_rho2 - vtou_rho2
  } else {
    rslt <- data.frame(utov_rho2, vtou_rho2)
  }
  rslt
}