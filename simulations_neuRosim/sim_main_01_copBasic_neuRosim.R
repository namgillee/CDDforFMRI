## Simulation of resting-state fMRI (RS-fMRI) data based on R packages copBasic and neuRosim
## : Generate sample fMRI time series and fitted curves
## by Namgil Lee & Jong-Min Kim
## 2018.09.26.


rm(list=ls())

require(copBasic)
require(neuRosim)
library(copula)
library(gcmr)
library(boot)
#library(igraph)
source('../functions/copuladd_compute.R')
source('../functions/copuladd_compute_2d.R')
source('../functions/empiric_df.R')

set.seed(1111)


## Generate underlying neural activities in RS-fMRI which are correlated & directed
Betas  = c(0.1, 0.3, 0.5)

totaltime <- 60*15   # time in seconds
acc <- 0.1           # higher resolution
effectsize=1         # amplitude parameter to multiply after convolution (see, demo_neuRosim.R)
fwhm = 4
TR = 2

n = totaltime / acc  # number of time points to generate


DFMRIS = vector('list',length(Betas))

for (idb in 1:length(Betas)) {
  bet = Betas[idb]
  print(paste('beta =',bet))
  
  
  para <- list(cop1=PLACKETTcop,  cop2=PLACKETTcop, # Note the singularity
               para1=5000, para2=5, alpha=1-bet, beta=bet)
  
  
  # bivariate data on unit interval [0,1]
  for (k in 1:100) {
    print(k)
    D <- simCOP(n=n, cop=composite2COP, para=para, cex=0.5)
    
    #check if any value is infinite
    xtmp = c(D[,1],D[,2])
    if ( !any(is.infinite( xtmp ))  && min(xtmp)>0 && max(xtmp)<1  )
      break
  }
  ######################## save figure #####################################
  tiff(paste0('resu_sim_main_01_scatter_D_beta0',bet*10,'.tiff'))
  plot(D, cex.lab=1.3, cex=0.5, main=paste("Pearson's r =",round(cor(D)[2],2)))
  dev.off()
  #pdf(paste0('resu_sim_main_01_scatter_D_beta0',bet*10,'.pdf'))
  #plot(D, cex.lab=1.3, cex=0.5, main=paste("Pearson's r =",round(cor(D)[2],2)))
  #dev.off()
  #print(cor(D))
  ##########################################################################
  
  
  # bivaraite data transformed to standard normal distribution
  Dnorm = matrix(0, n, 2)
  for (i in 1:n) {
    for (j in 1:2) {
      Dnorm[i,j] <- qnorm(D[i,j])
    }
  }
  colnames(Dnorm) <- c('X','Y')
  print('Correlation before convolution:')
  print(cor(Dnorm))
  ######################## save figure ####################################
  #tiff(paste0('resu_sim_main_01_scatter_Dnorm_beta0',bet*10,'.tiff'))
  #plot(Dnorm, cex.lab=1.3, cex=0.5, main=paste("Pearson's r =",round(cor(Dnorm)[2],2)), xlab='X', ylab='Y')
  #dev.off()
  #pdf(paste0('resu_sim_main_01_scatter_Dnorm_beta0',bet*10,'.pdf'))
  #plot(Dnorm, cex.lab=1.3, cex=0.5, main=paste("Pearson's r =",round(cor(Dnorm)[2],2)), xlab='X', ylab='Y')
  #dev.off()
  ##########################################################################
  
  
  ## Simulate RS-fMRI BOLD signal by convolution with Gaussian HRF
  Dconv = matrix(0, n, 2)
  for (j in 1:2) {
    s.conv <- convolve(gammaHRF(seq(acc, totaltime, acc), fwhm, verbose = FALSE), rev(Dnorm[,j]))
    Dconv[,j] <- effectsize * s.conv
  }
  colnames(Dconv) <- c('X','Y')
  Dconv <- data.frame(Dconv)
  print('Correlation after convolution:')
  print(cor(Dconv))
  
  
  plot.ts(cbind(Dnorm,Dconv))
  
  
  ## Subsampling by fMRI scanning
  Dfmri <- Dconv[seq(TR/acc, n, TR/acc),]
  colnames(Dfmri) <- c("X","Y")
  print('Correlation of fMRI BOLD:')
  print(cor(Dfmri))
  
  
  DFMRIS[[idb]] <- Dfmri
  
}  


#*********************************************************************************************#
##### Plot simulated fMRI data  ###############################################################
tiff(paste0('resu_sim_main_01_fitted_v_to_u.tiff'), width=960, height=960)
par(mfrow=c(3,length(Betas)))
for (idb in 1:length(Betas)) {
  Dfmri = DFMRIS[[idb]]
  bet = Betas[idb]
  
  Emp.index = matrix(0,nrow(Dfmri),2)
  colnames(Emp.index) <- c('u','v')
  Emp.index <- data.frame(Emp.index)
  
  ## Compute PDF, CDF, and iCDF
  ## Compute CDF transformation of fMRI data
  pdf=vector('list',2); 
  cdf=vector('list',2);
  icdf=vector('list',2)
  for (i in 1:2) { # use i-th column of Dfmri
    #kernel density estimation/smoothing
    pdf[[i]] <- density(Dfmri[,i]) 
    # Interpolate the density
    f <- approxfun(pdf[[i]]$x, pdf[[i]]$y, yleft=0, yright=0) 
    # Get the cdf by numeric integration
    cdf[[i]] <- function(x){
      integrate(f, -Inf, x)$value
    }
    # Use a root finding function to invert the cdf
    icdf[[i]] <- function(q){
      xval = seq(-5,5,0.01)
      yval = NULL
      for (kk in 1:length(xval))
        yval[kk] <- cdf[[i]](xval[kk])
      r = max(which(yval <= q))
      if (!is.null(r))
        return(xval[r])
      else
        return(-10)
        
      # uniroot(function(x){cdf[[i]](x) - q}, range(x))$root
    }
    
    # Transformation
    for (j in 1:nrow(Dfmri)) 
      Emp.index[j,i] <- cdf[[i]](Dfmri[j,i])
  }
  
  
  ## Estimate CDD with scatterplot/fitted curve (2D case)
  r12<-gcmr( u~v, data = Emp.index, marginal =beta.marg(link = "logit"),
             cormat = arma.cormat(0, 0) )
  Er12<-exp(r12$estimate[1]+Emp.index$v*r12$estimate[2])/(1+exp(r12$estimate[1]+Emp.index$v*r12$estimate[2]))
  r21<-gcmr(v ~u,  data = Emp.index, marginal =beta.marg(link = "logit"),cormat = arma.cormat(0, 0) )
  Er21<-exp(r21$estimate[1]+Emp.index$u*r21$estimate[2])/(1+exp(r21$estimate[1]+Emp.index$u*r21$estimate[2]))
  
  #plot v->u
  #tiff('resu_sim_main_01_fitted_v_to_u_beta0',bet*10,'.tiff')
  par(mar=c(5.1,6.1,4.1,2.1))
  plot(Dfmri[,c(2,1)], main=paste('beta =',bet), cex.lab = 3.2, cex.main=2.4, cex=1.2, ylim=c(-5.0, 4.5), xlim=c(-4,4))
  x = seq(0.01, 0.99, by=0.01)
  m = r12$estimate[1]+x*r12$estimate[2]
  y = exp(m)/(1+exp(m))
  
  x2 = NULL
  for (kk in 1:length(x))
    x2[kk] <- icdf[[1]](x[kk])
  y2 = NULL
  for (kk in 1:length(y))
    y2[kk] <- icdf[[2]](y[kk])
  lines(x2  ,  y2,  col='red',lwd=3)
  #dev.off()
  
  
} 
for (idb in 1:length(Betas)) {
  Dfmri = DFMRIS[[idb]]
  par(mar=c(1.1,6.1,4.1,2.1)) #par(mar=c(5.1,4.1,4.1,2.1))
  plot(Dfmri[,1], type='l', main="", xlab="", ylab="X", cex.lab=2.6, ylim=c(-5.0, 4.5), xaxt = "n") #suppress x axis
  
}
for (idb in 1:length(Betas)) {
  Dfmri = DFMRIS[[idb]]
  par(mar=c(5.1,6.1,0.1,2.1)) #par(mar=c(5.1,4.1,4.1,2.1))
  plot(Dfmri[,2], type='l', main="", xlab="Time", ylab="Y", cex.lab=2.6, ylim=c(-4,4)) #suppress x axis
  #axis(side = 1, at = seq(0,500,100), labels = seq(0,1000,200))
}
dev.off()
#*********************************************************************************************#

#*********************************************************************************************#
## Plot simulated fMRI data 2: only fitting curves
## use Betas[] and DFMRIS[[]]
#pdf(paste0('resu_sim_main_01_fitted_v_to_u_threeblk.pdf'), width=96/8, height=30/8)
tiff(paste0('resu_sim_main_01_fitted_v_to_u_threeblk.tiff'), width=960, height=300)
par(mfrow=c(1,length(Betas)))
for (idb in 1:length(Betas)) {
  Dfmri = DFMRIS[[idb]]
  bet = Betas[idb]
  
  Emp.index = matrix(0,nrow(Dfmri),2)
  colnames(Emp.index) <- c('u','v')
  Emp.index <- data.frame(Emp.index)
  
  ## Compute PDF, CDF, and iCDF
  ## Compute CDF transformation of fMRI data
  pdf=vector('list',2); 
  cdf=vector('list',2);
  icdf=vector('list',2)
  for (i in 1:2) { # use i-th column of Dfmri
    #kernel density estimation/smoothing
    pdf[[i]] <- density(Dfmri[,i]) 
    # Interpolate the density
    f <- approxfun(pdf[[i]]$x, pdf[[i]]$y, yleft=0, yright=0) 
    # Get the cdf by numeric integration
    cdf[[i]] <- function(x){
      integrate(f, -Inf, x)$value
    }
    # Use a root finding function to invert the cdf
    icdf[[i]] <- function(q){
      xval = seq(-5,5,0.01)
      yval = NULL
      for (kk in 1:length(xval))
        yval[kk] <- cdf[[i]](xval[kk])
      r = max(which(yval <= q))
      if (!is.null(r))
        return(xval[r])
      else
        return(-10)
      
      # uniroot(function(x){cdf[[i]](x) - q}, range(x))$root
    }
    
    # Transformation
    for (j in 1:nrow(Dfmri)) 
      Emp.index[j,i] <- cdf[[i]](Dfmri[j,i])
  }
  
  
  ## Estimate CDD with scatterplot/fitted curve (2D case)
  r12<-gcmr( u~v, data = Emp.index, marginal =beta.marg(link = "logit"),
             cormat = arma.cormat(0, 0) )
  Er12<-exp(r12$estimate[1]+Emp.index$v*r12$estimate[2])/(1+exp(r12$estimate[1]+Emp.index$v*r12$estimate[2]))
  r21<-gcmr(v ~u,  data = Emp.index, marginal =beta.marg(link = "logit"),cormat = arma.cormat(0, 0) )
  Er21<-exp(r21$estimate[1]+Emp.index$u*r21$estimate[2])/(1+exp(r21$estimate[1]+Emp.index$u*r21$estimate[2]))

  #plot v->u
  #tiff('resu_sim_main_01_fitted_v_to_u_beta0',bet*10,'.tiff')
  par(mar=c(5.1,6.1,4.1,2.1))
  plot(Dfmri[,c(2,1)], main=paste('beta =',bet), cex.lab = 3.2, cex.main=2.4, cex=1.2, ylim=c(-5.0, 4.5), xlim=c(-4,4))
  x = seq(0.01, 0.99, by=0.01)
  m = r12$estimate[1]+x*r12$estimate[2]
  y = exp(m)/(1+exp(m))
  
  x2 = NULL
  for (kk in 1:length(x))
    x2[kk] <- icdf[[1]](x[kk])
  y2 = NULL
  for (kk in 1:length(y))
    y2[kk] <- icdf[[2]](y[kk])
  lines(x2  ,  y2,  col='red',lwd=3)
  #dev.off()
} 
dev.off()
#*********************************************************************************************#


#************************ Procedure of Transformations: fixed beta ***************************#
bet = 0.1
print(paste('beta =',bet))
para <- list(cop1=PLACKETTcop,  cop2=PLACKETTcop, # Note the singularity
             para1=5000, para2=5, alpha=1-bet, beta=bet)
# bivariate data on unit interval [0,1]
for (k in 1:100) {
  print(k)
  D <- simCOP(n=n, cop=composite2COP, para=para, cex=0.5)
  
  #check if any value is infinite
  xtmp = c(D[,1],D[,2])
  if ( !any(is.infinite( xtmp ))  && min(xtmp)>0 && max(xtmp)<1  )
    break
}

# bivaraite data transformed to standard normal distribution
Dnorm = matrix(0, n, 2)
for (i in 1:n) {
  for (j in 1:2) {
    Dnorm[i,j] <- qnorm(D[i,j])
  }
}
colnames(Dnorm) <- c('X','Y')

## Simulate RS-fMRI BOLD signal by convolution with Gaussian HRF
Dconv = matrix(0, n, 2)
for (j in 1:2) {
  s.conv <- convolve(gammaHRF(seq(acc, totaltime, acc), fwhm, verbose = FALSE), rev(Dnorm[,j]))
  Dconv[,j] <- effectsize * s.conv
}
colnames(Dconv) <- c('X','Y')
Dconv <- data.frame(Dconv)

## Subsampling by fMRI scanning
Dfmri <- Dconv[seq(TR/acc, n, TR/acc),]
colnames(Dfmri) <- c("X","Y")

### Plot D, Dnorm, Dconv, Dfmri ###
pdf(paste0('resu_sim_main_01_transform_beta0',bet*10,'.pdf'), width=64/7, height=96/7)
#tiff(paste0('resu_sim_main_01_transform_beta0',bet*10,'.tiff'), width=640, height=960)
par(mfrow=c(3,2))

par(mar=c(5.1,5.1,4.1,2.1))
plot(D, cex.lab=1.3, cex=0.5, main=paste("(a) Original data on unit interval"), cex.main=2, cex.lab=2)

par(mar=c(5.1,5.1,4.1,2.1))
plot(Dnorm, cex.lab=1.3, cex=0.5, main=paste("(b) Transformed data to Normal Dist."), xlab='X', ylab='Y', cex.main=2, cex.lab=2)

par(mar=c(1.1,5.1,4.1,2.1)) #par(mar=c(5.1,4.1,4.1,2.1))
plot(Dconv[,1], type='l', main="(c) Simulated BOLD via Gaussian HRF", xlab="", ylab="X", cex.lab=1.3, ylim=c(-5.0, 4.5), xaxt = "n", cex.main=2, cex.lab=2) 
par(mar=c(1.1,5.1,4.1,2.1)) #par(mar=c(5.1,4.1,4.1,2.1))
plot(Dfmri[,1], type='l', main="(d) Observed FMRI at TR=2(sec)", xlab="", ylab="X", cex.lab=1.3, ylim=c(-5.0, 4.5), xaxt = "n", cex.main=2, cex.lab=2)

par(mar=c(5.1,5.1,0.1,2.1)) #par(mar=c(5.1,4.1,4.1,2.1))
plot(Dconv[,2], type='l', main="", xlab="Time", ylab="Y", cex.lab=1.3, ylim=c(-4,4), cex.lab=2) 
par(mar=c(5.1,5.1,0.1,2.1)) #par(mar=c(5.1,4.1,4.1,2.1))
plot(Dfmri[,2], type='l', main="", xlab="Time", ylab="Y", cex.lab=1.3, ylim=c(-4,4), cex.lab=2) 

dev.off()
par(mfrow=c(1,1))
#*****************************************************************************************#



save(DFMRIS, Betas, 
      file = 'resu_sim_main_01_sampleplots.RData')


# 
# 
# ## Draw sample time series 
# pdf('resu_sim_main_01_timeseries_Dfmri_beta0',bet*10,'.pdf')
# dev.off()
# tiff('resu_sim_main_01_timeseries_Dfmri_beta0',bet*10,'.tiff')
# plot.ts(Dfmri, main="", xlab="Time", cex.lab=1.3)
# dev.off()


## Estimate Copula Directional Dependence (CDD) with bootstrap *******************#
#resu_copuladd <- copuladd_compute(data = Dfmri, bootstrap.number = 100)
#print(resu_copuladd)
# > print(resu_copuladd)
# node_from node_to dd_from2to dd_to2from    dd_diff  dd_diff_lb dd_diff_ub
# 1         2       1   0.370889  0.3586669 0.01222203 -0.03687544 0.07335026
#*********************************************************************************#








