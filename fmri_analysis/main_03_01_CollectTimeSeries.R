## Collect time series at selected voxel locations from preprocessed fMRI data
# We will collect data at 14 brain regions
#
# 2018.09.27. Namgil Lee & Jong-Min Kim

rm(list=ls())

## Load packages
library(oro.nifti)   
source('../functions/BA_MNI.R')

subj_start = 1
subj_end = 155

#***********************************************************#
## Results will be saved in this directory                  #
#***********************************************************#
## Preprocessed functional image data
inputdata_dir = './resu_fmridata_preprocessed/'
resu_dir = './resu_main_03_01/'

sink(file = file(paste0(resu_dir, "output", subj_start, ".txt")))


## Iterate over all subjects
for (ids in subj_start:subj_end) {
  strsubj <- sprintf('%03d',ids)
  print(paste('** Subject', strsubj, '**'))
  
 if (!file.exists(paste0(resu_dir, 'resu_main_', strsubj, '_time_series.Rdata'))) {
  #***********************************************************#
  #                  Read and get data                        #
  #***********************************************************#
  
  
  funcfname <- paste0(inputdata_dir, 'swasub-pixar', strsubj, '_task-pixar_bold.nii')
  mydat <- readNIfTI(funcfname, reorient = FALSE)
  # mymask <- (mydat > quantile(mydat, 0.25))   # TOO Large
  thres_mymask <- quantile(mydat, 0.25)
  mydim = dim(mydat)  
  
  
  #**********************************************************#
  #                  Brodmann Areas                          #
  #  http://sprout022.sprout.yale.edu/mni2tal/mni2tal.html   #
  #**********************************************************#
  myROIALL = BA_MNI #3 x 82
  myROIALL = (myROIALL + c(90, 126, 72) + 1) * c(91,109,91) / c(181,217,181)
  myROIALL = round(myROIALL)
  
  
  ## Skip checking maximum distance
  # mydistALL = dist(t(myROIALL), method = 'maximum')
  # hist(mydistALL, 100, main='Histogram of Maxumum Distances Between BAs', xlab = 'Maximum distance')
  
  
  ## Select time series based on BA's
  # BA1(PARIETAL),8(FRONTAL),10(FRONTAL), BA18(OCCIPITAL), 
  # 20(TEMPORAL), BA23(MEDIAL), 38(MEDIAL)
  selectedindex = sapply(colnames(myROIALL), substring, 3, 6) %in% 
    c('BA04', 'BA08', 'BA10', 'BA18', 'BA20', 'BA23', 'BA38')
  myROI = myROIALL[,selectedindex]
  # print('Selected BAs: '); names(myROI)
  
  
  ## Skip verifying distances between ROIs
  ## It is OK if the minimum of the 'maximum distances' is at least 5
  # mydist_R = dist(t(myROI[,seq(1,ncol(myROI),by=2)]), method = 'maximum')
  # mydist_L = dist(t(myROI[,seq(2,ncol(myROI),by=2)]), method = 'maximum')
  # print(paste('The minimum distance bet. ROIs: ', min(mydist_R, mydist_L)))
  # hist(c(mydist_R, mydist_L), 100, main='Histogram of Maxumum Distances Between selected BAs', xlab = 'Maximum distance')

    
  ## We are using preprocessed (but non-smoothed) data
  ## Here, we consider neighbors to be averaged
  nbsize = 2
  # print(paste('We will use neighbors within plus-minus ', nbsize))
  # if (nbsize > (min(mydist_R, mydist_L)-1)/2) {
  #   print('Caution: Some neighbors overlap.')
  # }
  
  
  #**********************************************************#
  #        Extract time series on ROI                        #
  #**********************************************************#
  tsdat_roi = matrix(0,mydim[4],ncol(myROI))
  tsdat_roi_name = character(ncol(myROI))
  cur_roi <- 1
  while (cur_roi <= ncol(myROI)) {
    ## Sweep over ROIs
      
    #print(paste('In voxels:', names(myROI)[cur_roi]))

    ## Collect time series on a neighbor in 3D cube
    tsdat_cube = NULL   
    aidR = myROI[,cur_roi]  # arrayindex at the center, in 3D coordinate
    mydat_ptR_num = mydat[(aidR[1]-nbsize):(aidR[1]+nbsize),
                          (aidR[2]-nbsize):(aidR[2]+nbsize), 
                          (aidR[3]-nbsize):(aidR[3]+nbsize), 1:mydim[4]]

    
    for (id1 in 1:(1+2*nbsize)) {
      for (id2 in 1:(1+2*nbsize)) {
        for (id3 in 1:(1+2*nbsize)) {

          if ( mean(tsmydat <- mydat_ptR_num[id1,id2,id3,]) > thres_mymask )
            tsdat_cube = cbind(tsdat_cube, tsmydat)
        }
      }
    }
    
    
    # time series - averaged
    if (is.null(tsdat_cube)) {
      warning(paste('No data available at:', 
                    names(myROI)[cur_roi]))
    } else {
      print(paste('Number of averaged timeseries:',ncol(tsdat_cube)) )
      
      idx_half <- cur_roi%%2
      if (idx_half == 1) {
        #odd => Right Hemisphere
        tsdat_roi[ ,(cur_roi+1)/2] <- rowMeans(tsdat_cube)
        tsdat_roi_name[(cur_roi+1)/2] <- paste0("R.",substring(names(myROI)[cur_roi], 8))
        
      } else {
        #even =>Left Hemisphere
        tsdat_roi[ ,15-(cur_roi/2)] <- rowMeans(tsdat_cube)
        tsdat_roi_name[15-(cur_roi/2)] <- paste0("L.",substring(names(myROI)[cur_roi], 8))
        
      }
        
    }
    
  
    cur_roi = cur_roi + 1
  }
  
  
  #rm(mydat, mydat_ptR_num, mydat_ptL_num)
  
  
  colnames(tsdat_roi) <- tsdat_roi_name
  
  save(file = paste0(resu_dir, 'resu_main_', strsubj, '_time_series.Rdata'), tsdat_roi, tsdat_roi_name)
  
 }
}


sink()