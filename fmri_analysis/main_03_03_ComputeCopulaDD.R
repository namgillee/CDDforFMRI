## Compute Copula Directional Dependence
#
# 2018.09.27. Namgil Lee & Jong-Min Kim


rm(list=ls())

## Load packages
library(copula)
library(gcmr)
source('../functions/copuladd_compute.R')
source('../functions/copuladd_compute_2d.R')
source('../functions/empiric_df.R')
source('../functions/draw_copuladd_network.R')

nskip = 5 #number of time points to skip, counted from the beginning
  
subj_start = 1
subj_end = 155

#***********************************************************#
## Results will be saved in this directory                  #
#***********************************************************#
resu_ts_dir = './resu_main_03_01/'
resu_dir = './resu_main_03_03/'

sink(file = file(paste0(resu_dir, "c_output", subj_start, ".txt")))


## Iterate over all subjects
for (ids in subj_start:subj_end) {
  strsubj <- sprintf('%03d',ids)
  print(paste('** Subject', strsubj, '**'))
  
 if (!file.exists(paste0(resu_dir, 'resu_main_', strsubj, '_copuladd.Rdata'))) {
  
  #***********************************************************#
  #                  Read time_series data                    #
  #***********************************************************#
  load(paste0(resu_ts_dir, 'resu_main_', strsubj, '_time_series.Rdata'))
  ## tsdat_roi, tsdat_roi_name
   


  #**********************************************************#
  #               Compute Copula DD                          #
  #  By using 'tsdat_roi' and 'tsdat_roi_name'               #
  #  We remove the first '11' points (0~10 TRs)              #
  #  due to artifact.                                        #
  #**********************************************************#
  set.seed(100)

  start_time = Sys.time()
  idxcol_nonnull <- (tsdat_roi_name != "") ##REMOVE MISSING DATA; SEE, "main_03_01...R"
  resu_copuladd <- copuladd_compute(data = tsdat_roi[-(1:nskip),idxcol_nonnull], bootstrap.number = 100)
  end_time = Sys.time()

  time_taken = difftime(end_time, start_time, units='secs')

  print(paste('Total time taken:', round(time_taken,3), 'secs'))
  save(file = paste0(resu_dir, 'resu_main_', strsubj, '_copuladd.Rdata'), resu_copuladd)
  #----------------------------------------------------------#


  #**********************************************************#
  #             Select significant DD                        #
  #**********************************************************#
  resu_copuladd.sig <- resu_copuladd[resu_copuladd$dd_diff_lb>0,,drop=FALSE]

  print(paste('** Number of significant edges:', nrow(resu_copuladd.sig)))

  
 }

}


sink()