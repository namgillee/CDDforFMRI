copuladd_compute_2d <- function(dat, idx = 1:nrow(dat), out_diff = TRUE)
{
  # < CDD for FMRI >
  #   
  #   Copyright (C) 2018 Namgil Lee & Jong-Min Kim
  
  #dat is a [nx2] matrix
  
  dat <- as.data.frame(dat[idx,1:2,drop=FALSE])
  
  colnames(dat) <- c("u","v")
  
  r12<-gcmr( u~v, data = dat, marginal = beta.marg(link = "logit"), cormat = arma.cormat(0, 0) )
  Er12<-exp(r12$estimate[1]+dat$v*r12$estimate[2])/(1+exp(r12$estimate[1]+dat$v*r12$estimate[2]))
  vtou_rho2<-var(Er12)/var(dat$u)

  r21<-gcmr( v~u, data = dat, marginal = beta.marg(link = "logit"), cormat = arma.cormat(0, 0) )
  Er21<-exp(r21$estimate[1]+dat$u*r21$estimate[2])/(1+exp(r21$estimate[1]+dat$u*r21$estimate[2]))
  utov_rho2<-var(Er21)/var(dat$v)

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