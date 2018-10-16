## sim_VAR3D.R
# Simulation for VAR(1) process with d=3 variables
# Evaluate Copula DD with MSEs for different VAR(1) models
#
# models to test:
#   VAR(1) with noise from multivariate normal (1,2,3 are independent)
#   VAR(1) with noise from multivariate normal (1-2 are correlated, 3 independent)
#   VAR(1) with noise from multivariate normal (1<-2 are directionally correlated, 3 independent)
#
# Last modified at Sep. 28, 2018
# by Namgil Lee & Jong-Min Kim

rm(list=ls())

set.seed(111)

library(MASS)
require(copBasic)
library(copula)
library(gcmr)
#library(boot)
#library(igraph)
source('../functions/copuladd_compute.R')
source('../functions/copuladd_compute_2d.R')
source('../functions/empiric_df.R')


## Create VAR(1) coefficient matrix with d=3 variables
d = 3
r = 0.5
A = r * matrix( c(0,1,0, 0,0,1, 1,0,0), d, d, byrow = TRUE)
print('The coefficient matrix A is: ')
print(round(A,2))
# > print(round(A,2))
#      [,1] [,2] [,3]
# [1,]  0.0  0.5  0.0
# [2,]  0.0  0.0  0.5
# [3,]  0.5  0.0  0.0


## Theoretical coefficient matrix of lagged VAR(1) process (lag=5)
lag = 5
print(paste('The lag L is:', lag))
print('The true coef matrix A^L for lagged VAR(1) model:')
A_lag = diag(d); for ( i in 1:lag ) A_lag = A_lag%*%A ; 
print(round(A_lag,5))
# > print(round(A_lag,5))
#         [,1]    [,2]    [,3]
# [1,] 0.00000 0.00000 0.03125
# [2,] 0.03125 0.00000 0.00000
# [3,] 0.00000 0.03125 0.00000


## Generate time series data from VAR(1) model: y = A*x + e
## e_i has mean 0 and variance 1
gen_copula_data <- function(len, beta=0.3) {
  para <- list(cop1=PLACKETTcop,  cop2=PLACKETTcop, 
               para1=5000, para2=5, alpha=1-beta, beta=beta)
  for (k in 1:100) {
    D <- simCOP(n=len, cop=composite2COP, para=para, cex=0.5)
    xtmp = c(D[,1],D[,2])
    if ( !any(is.infinite( xtmp ))  && min(xtmp)>0 && max(xtmp)<1  )
      break
  }
  Dnorm = matrix(0, len, 2)
  for (i in 1:len) {
    for (j in 1:2) {
      Dnorm[i,j] <- qnorm(D[i,j])
    }
  }
  return(Dnorm)
}
gen_noise_demo <- function(len, type, d=1) {
  # noise of length len, dimension d
  
  if (type == 1) {
    # mvrnorm; independent noise
    tsdat = mvrnorm(len, mu = rep(0,d), Sigma = diag(d))  
    
  } else if (type == 2) {
    # mvrnorm; (variables 1-2 dependent, other variables are independent) noise
    S <- diag(d)
    S[1,2] <- S[2,1] <- 0.56
    tsdat = mvrnorm(len, mu = rep(0,d), Sigma = S)  
    
  } else if (type == 3) {
    # (2->1,3)
    Dnorm <- gen_copula_data(len, beta = 0.3) #2->1, 
    tsdat = cbind(Dnorm, mvrnorm(len, rep(0,d-2), diag(1,d-2) ))
  } 
  
  return(tsdat)
}


# Check the correlation of noises from asymmetric Copula
testdata <- gen_copula_data(10000, beta=0.3)
print(cor(testdata)) 
# > print(cor(testdata))
#           [,1]      [,2]
# [1,] 1.0000000 0.5582104
# [2,] 0.5582104 1.0000000


################# SIMULATION #####################################
len_ts0 = 10000   #For original VAR(1) process before lagged
All_Copula_Results = vector('list',3)
names(All_Copula_Results) <- c("Normal_Independent","Normal_Dependent","Normal_Directional")

##### 1. Normal #################
len_ts = len_ts0
tsdat <- gen_noise_demo(len_ts+1, 1, d)

#  generate by VAR process
for (t in 1:len_ts) {
  tsdat[t+1,] = A %*% tsdat[t,] + tsdat[t+1,]
}
tsdat_obs = tsdat[seq(lag+1,len_ts+1,by=lag),]

print(paste('Correlation of y1 and y2 is:', cor(tsdat[,c(1,2)])[2]))
#[1] "Correlation of y1 and y2 is: -0.00728073221419077"

# Compute Copula  
resu_copuladd <- copuladd_compute(data = tsdat_obs, bootstrap.number = 100)
All_Copula_Results$Normal_Independent <- resu_copuladd
print(resu_copuladd)
#   node_from node_to   dd_from2to   dd_to2from      dd_diff   dd_diff_lb   dd_diff_ub
# 1         1       2 2.423731e-04 1.127441e-04 1.296290e-04 5.392030e-05 0.0001763885 
# 2         1       3 8.543149e-05 9.036397e-06 7.639509e-05 2.884985e-06 0.0001233765
# 3         2       3 9.204221e-05 2.396941e-05 6.807280e-05 3.964810e-05 0.0001286255
### <-- The Lower-Bounds are very close to zero (accuracy of <1e-4) ==> not clear directionality
### <-- Directions are all mixed

### investigate edge (1,2): length(x_norm) is 3; The boot results are in x_norm[[1]]
x_norm <- attr(All_Copula_Results$Normal_Independent, 'boot_diff')

mn_norm = mean(x_norm[[1]]); 
sd_norm = sd(x_norm[[1]])/sqrt(length(x_norm[[1]])); 
print(paste("Significance of alpha:", 2*pnorm(0, abs(mn_norm), sd_norm))) #Check if it is smaller than two-sided alpha=0.05, 0.01
#[1] "Significance of alpha: 0.000190388842927035"

bias_norm = mean(x_norm[[1]]) - All_Copula_Results$Normal_Independent$dd_diff[1]  #MSE 
bias_norm
#[1] -1.447456e-05
sd_norm = sd(x_norm[[1]])
sd_norm
#[1] 0.0003086058
mse_norm = bias_norm^2 + sd_norm^2
mse_norm
#[1] 9.544706e-08
rmse_norm = sqrt(mse_norm)
rmse_norm
#[1] 0.0003089451
### <-- All_Copula_Results$Normal_Independent$dd_diff[1] == 1.296290e-04 < rmse_norm == 0.0003089451
rrmse_norm <- rmse_norm / abs(All_Copula_Results$Normal_Independent$dd_diff[1])
print(paste("Relative RMSE:",rrmse_norm))
###     Relative RMSE: rmse_norm / abs(All_Copula_Results$Normal_Independent$dd_diff[1])  == 2.383303
####################################################################################################


##### 2. Normal (dependent) #####
len_ts = len_ts0
tsdat <- gen_noise_demo(len_ts+1, 2, d)

#  generate by VAR process
for (t in 1:len_ts) {
  tsdat[t+1,] = A %*% tsdat[t,] + tsdat[t+1,]
}
tsdat_obs = tsdat[seq(lag+1,len_ts+1,by=lag),]

print(paste('Correlation of y1 and y2 is:', cor(tsdat[,c(1,2)])[2]))
#[1] "Correlation of y1 and y2 is: 0.424721614807866"


# Compute Copula  
resu_copuladd <- copuladd_compute(data = tsdat_obs, bootstrap.number = 100)
All_Copula_Results$Normal_Dependent <- resu_copuladd
print(resu_copuladd)
#   node_from node_to   dd_from2to   dd_to2from      dd_diff   dd_diff_lb   dd_diff_ub
# 1         1       2 0.1598391577 0.1596897437 0.0001494140 -0.0008437488 0.0009079963
# 2         3       1 0.0083513769 0.0080112313 0.0003401457  0.0001090376 0.0005275607  *
# 3         2       3 0.0007102223 0.0002285543 0.0004816680  0.0004172277 0.0005966314  *
### <-- The Lower-Bounds are very close to zero (accuracy <1e-3) ==> not clear directionality
###     Especially for direction bet. 1&2, LB < 0. 
### <-- Directions are all reversed. 

### investigate edge (1,2): length(x_dep) is 3; The boot results are in x_dep[[1]]
x_dep <- attr(All_Copula_Results$Normal_Dependent, 'boot_diff')

mn_dep = mean(x_dep[[1]]); 
sd_dep = sd(x_dep[[1]])/sqrt(length(x_dep[[1]])); 
print(paste("Significance of alpha:", 2*pnorm(0, abs(mn_dep), sd_dep))) #Check if it is smaller than two-sided alpha=0.05, 0.01
#[1] "Significance of alpha: 0.941986173612399

bias_x_dep = mean(x_dep[[1]]) - All_Copula_Results$Normal_Dependent$dd_diff[1]  #MSE 
bias_x_dep
#[1] -0.0001172903
sd_x_dep = sd(x_dep[[1]])
sd_x_dep
#[1] 0.004414197
mse_x_dep = bias_x_dep^2 + sd_x_dep^2
mse_x_dep
#[1] 1.94989e-05
rmse_x_dep = sqrt(mse_x_dep)
rmse_x_dep
#[1] 0.004415755
### <-- All_Copula_Results$Normal_Dependent$dd_diff[1] == 0.000149414 < rmse_x_dep == 0.004415755
rrmse_dep <- rmse_x_dep / abs(All_Copula_Results$Normal_Dependent$dd_diff[1])
print(paste("Relative RMSE:",rrmse_dep))
###     Relative RMSE: rmse_x_dep / abs(All_Copula_Results$Normal_Dependent$dd_diff[1])  == 29.55382
####################################################################################################


##### 3. (2->1,3) direction #####
len_ts = len_ts0
tsdat <- gen_noise_demo(len_ts+1, 3, d)

#  generate by VAR process
for (t in 1:len_ts) {
  tsdat[t+1,] = A %*% tsdat[t,] + tsdat[t+1,]
}
tsdat_obs = tsdat[seq(lag+1,len_ts+1,by=lag),]

print(paste('Correlation of y1 and y2 is:', cor(tsdat[,c(1,2)])[2]))
#[1] "Correlation of y1 and y2 is: 0.41985705953035"

# Compute Copula  
resu_copuladd <- copuladd_compute(data = tsdat_obs, bootstrap.number = 100)
All_Copula_Results$Normal_Directional <- resu_copuladd
print(resu_copuladd)
#   node_from node_to  dd_from2to  dd_to2from      dd_diff  dd_diff_lb   dd_diff_ub
# 1         2       1 0.161447524 0.152748489 0.0086990345 0.007359009 0.0093386327   **
# 2         1       3 0.016391661 0.015461532 0.0009301290 0.000614440 0.0012257971   *
# 3         2       3 0.002573261 0.002060448 0.0005128132 0.000369759 0.0006228495   *
### <-- Edge 2->1 Show significant directionality (LB > 5e-3)
###     Other edges are not clear (LB < 1e-3)

### investigate edge (1,2): length(x_dep) is 3; The boot results are in x_dep[[1]]
x_dir <- attr(All_Copula_Results$Normal_Directional, 'boot_diff')

mn_dir = mean(x_dir[[1]]); 
sd_dir = sd(x_dir[[1]])/sqrt(length(x_dir[[1]])); 
print(paste("Significance of alpha:", 2*pnorm(0, abs(mn_dir), sd_dir) )) #Check if it is smaller than two-sided alpha=0.05, 0.01
#[1] "Significance of alpha: 7.11778930348215e-63

bias_x_dir = mean(x_dir[[1]]) + All_Copula_Results$Normal_Directional$dd_diff[1]  #MSE *** CAREFUL OF THE SIGNS OF ITS TERMS
bias_x_dir
#[1] 0.0003502135
sd_x_dir = sd(x_dir[[1]])
sd_x_dir
#[1] 0.004988425
mse_x_dir = bias_x_dir^2 + sd_x_dir^2
mse_x_dir
#[1] 2.500703e-05
rmse_x_dir = sqrt(mse_x_dir)
rmse_x_dir
#[1] 0.005000703
### <-- All_Copula_Results$Normal_Directional$dd_diff[1] == 0.008699035 > rmse_x_dir == 0.005000703
rrmse_dir <- rmse_x_dir / abs(All_Copula_Results$Normal_Directional$dd_diff[1])
print(paste("Relative RMSE:",rrmse_dir))
###     Relative RMSE: rmse_x_dir / abs(All_Copula_Results$Normal_Directional$dd_diff[1])  == 0.574857221498029 
####################################################################################################


######## Boxplot #######
tiff("sim_VAR3D_boxplot_resampled_diff.tiff")
boxdata <- cbind(x_norm[[1]], x_dep[[1]], x_dir[[1]])
colnames(boxdata) <- c("Independent", "Dependent", "Asymmetric")
boxplot(boxdata, 
        ylab = expression(paste(Delta, rho[U], .[V])), xlab = "Noise type", cex.lab = 1.5, 
        cex.axis = 1.3
)
grid()
dev.off()


##### save

save(All_Copula_Results, 
     file = 'resu_sim_VAR3D_1.RData')

