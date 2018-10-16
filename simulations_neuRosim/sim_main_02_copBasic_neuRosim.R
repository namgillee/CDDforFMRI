## Simulation of resting-state fMRI (RS-fMRI) data based on R packages copBasic and neuRosim
## Sensitivity analysis for various simulation parameters: 
##     copula param2(==beta)
##     data length in seconds
##
## The results can be presented in 
##     table of MSEs and other statistics
##     plot of MSEs in CDD vs Correlation
##
##     time cources (in sim main 01)
##     Scatter plot and fitted curves (in sim main 01)
## 
## by Namgil Lee & Jong-Min Kim
## 2018.09.24.


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

set.seed(2222)

## Copula parameters 
Betas = c(0.1,0.2,0.3,0.4,0.5)   # asymmetry/directionality weaker; correlation stronger; 
Times = 60 * c(5,10,15,20,25,30)    # RS-fMRI time in seconds

## fMRI simulated data parameters 
acc <- 0.1           # higher resolution
effectsize=1         # amplitude parameter to multiply after convolution (see, demo_neuRosim.R)
fwhm = 4
TR = 2   #seconds


resuCorDfmri <- array(0,c(length(Betas),length(Times)))
resuCDDutov <- resuCDDvtou <- resuCDDdiff <- array(0,c(length(Betas),length(Times)))
resuMN <- resuBIAS <- resuSE <- resuMSE <- resuLBCI <- resuUBCI <- array(0,c(length(Betas),length(Times)))


for (idb in 1:length(Betas)) {
  bet = Betas[idb]
  print(paste('beta:',bet))
  
  for (idt in 1:length(Times)) {
    tim = Times[idt]
    print(paste('time (sec):',tim))
    
    
    ## Generate underlying neural activities in RS-fMRI which are correlated & directed
    totaltime <- tim     # time in seconds

    n = totaltime / acc  # number of time points to generate

    para <- list(cop1=PLACKETTcop,  cop2=PLACKETTcop, # Note the singularity
                 para1=5000, para2=5, alpha=1-bet, beta=bet)
    
    # Generate bivariate data on unit interval [0,1]
    for (k in 1:100) {
      print(k)

      D <- simCOP(n=n, cop=composite2COP, para=para, cex=0.5)
      
      #check if any value is infinite
      xtmp = c(D[,1],D[,2])
      if ( !any(is.infinite( xtmp ))  && min(xtmp)>0 && max(xtmp)<1  )
        break
    }
    # k
    # cor(D)

    
    # bivaraite data transformed to standard normal distribution
    Dnorm = matrix(0, n, 2)
    for (i in 1:n) {
      for (j in 1:2) {
        Dnorm[i,j] <- qnorm(D[i,j])
      }
    }
    #print('Correlation before convolution:')
    #print(cor(Dnorm))
    
    
    
    ## Simulate RS-fMRI BOLD signal by convolution with Gaussian HRF
    Dconv = matrix(0, n, 2)
    for (j in 1:2) {
      s.conv <- convolve(gammaHRF(seq(acc, totaltime, acc), fwhm, verbose = FALSE), rev(Dnorm[,j]), type='circular')
      Dconv[,j] <- effectsize * s.conv
    }
    colnames(Dconv) <- c('X','Y')
    Dconv <- data.frame(Dconv)
    #print('Correlation after convolution:')
    #print(cor(Dconv))
    #plot.ts(cbind(Dnorm, Dconv))
    #Correlations[r,"Pearson"] <- cor(Dconv[,1], Dconv[,2])
    
  
    ## Subsampling by fMRI scanning
    Dfmri <- Dconv[seq(TR/acc, n, TR/acc),]
    resuCorDfmri[idb,idt] <- cor(Dfmri[,1],Dfmri[,2]) #***********************#
    tiff(paste0('resu_sim_main_02_scatter_Dfmri_',idb,'_',idt,'.tiff'))
    plot(Dfmri, cex=0.5, cex.lab=1.3)
    dev.off()
    
    
    ## Estimate Copula Directional Dependence (CDD) with bootstrap 
    resu_copuladd <- copuladd_compute(data = Dfmri, bootstrap.number = 100)
    resuCDDutov[idb,idt] <- resu_copuladd$dd_from2to 
    resuCDDvtou[idb,idt] <- resu_copuladd$dd_to2from
    resuCDDdiff[idb,idt] <- resu_copuladd$dd_diff
    
    boot_diff_data = attr(resu_copuladd, 'boot_diff')[[1]]
    resuMN[idb,idt] <- mean(boot_diff_data)   #*** boot_diff_Mean ***#
    resuBIAS[idb,idt] <- ifelse(resu_copuladd$node_from < resu_copuladd$node_to,
                                resuMN[idb,idt] - resuCDDdiff[idb,idt],
                                resuMN[idb,idt] + resuCDDdiff[idb,idt]
    ) #*** boot_diff_Bias ***#
    resuSE[idb,idt] <- sd(boot_diff_data)
    resuMSE[idb,idt] <- resuBIAS[idb,idt]^2 + resuSE[idb,idt]^2
    resuLBCI[idb,idt] <- resu_copuladd$dd_diff_lb #*** boot_diff_lower_bound ***#
    resuUBCI[idb,idt] <- resu_copuladd$dd_diff_ub #*** boot_diff_upper_bound ***#
    
  }
}

save(resuCorDfmri, 
     resuCDDutov, resuCDDvtou, resuCDDdiff, 
     resuMN, resuBIAS, resuSE, resuMSE, resuLBCI, resuUBCI,
     file = 'resu_sim_main_02_sensitivity.RData')


########################################################################
## Draw MSE curves
tiff('resu_sim_main_02_mse_lines.tiff')
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))
plot(Times, resuMSE[1,], type='o', ylim=c(0,1e-3), xaxt='n',
     pch=1, col=1, lty=1, lwd=2, cex=2, 
     xlab='Total time (sec)', ylab='Mean squared error', cex.lab=1.3)
axis(1,at = Times,labels = Times)
for (i in 2:5) {
  lines(Times, resuMSE[i,], type='o', 
        pch=i, col=i, lty=c(1,2,1,2,1)[i], lwd=2, cex=2)
}
legend('topright',legend = paste('beta =',Betas),
       pch=1:5, col=1:5, lty=c(1,2,1,2,1), lwd=2, cex=1.5)
dev.off()
# pdf('resu_sim_main_02_mse_lines.pdf')
# par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))
# plot(Times, resuMSE[1,], type='o', ylim=c(0,1e-3), xaxt='n',
#      pch=1, col=1, lty=1, lwd=2, cex=2, 
#      xlab='Total time (sec)', ylab='Mean squared error', cex.lab=1.3)
# axis(1,at = Times,labels = Times)
# for (i in 2:5) {
#   lines(Times, resuMSE[i,], type='o', 
#         pch=i, col=i, lty=c(1,2,1,2,1)[i], lwd=2, cex=2)
# }
# legend('topright',legend = paste('beta =',Betas),
#        pch=1:5, col=1:5, lty=c(1,2,1,2,1), lwd=2, cex=1.5)
# dev.off()
########################################################################


########################################################################
#Draw RMSE curve
tiff('resu_sim_main_02_rmse_lines.tiff')
par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))
plot(Times, sqrt(resuMSE[1,]), type='o', xaxt='n', ylim=c(0,0.03),
     pch=1, col=1, lty=1, lwd=2, cex=2, 
     xlab='Total time (sec)', ylab='Root Mean Squared Error', cex.lab=1.3)
axis(1,at = Times,labels = Times)
for (i in 2:5) {
  lines(Times, sqrt(resuMSE[i,]), type='o', 
        pch=i, col=i, lty=c(1,2,1,2,1)[i], lwd=2, cex=2)
}
legend('topright',legend = paste('beta =',Betas),
       pch=1:5, col=1:5, lty=c(1,2,1,2,1), lwd=2, cex=1.5)
dev.off()

