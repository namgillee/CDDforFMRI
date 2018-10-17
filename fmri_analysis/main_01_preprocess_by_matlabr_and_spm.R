## Preprocessing of fMRI data by using the Matlab package SPM12 and R package spm12r
# Before running, install Matlab and SPM12
# https://cran.r-project.org/web/packages/spm12r/vignettes/fmri_processing_spm12r.html

rm(list=ls())

## Load packages
library(oro.nifti)   
library(neurobase)   
library(matlabr)
library(spm12r)
# library(fmri)
source('../functions/BA_MNI.R')
library(copula)
library(gcmr)
library(boot)
library(igraph)
source('../functions/copuladd_compute.R')
source('../functions/copuladd_compute_2d.R')
source('../functions/empiric_df.R')
source('../functions/draw_copuladd_network.R')



if (!have_matlab()) {
  print('Cannot find Matlab..')
}

#***********************************************************#
## Results will be saved in this directory                  #
#***********************************************************#
fmridata_dir = './openfmri_ds000228_animated_film/ds000228_R1.0.1/'
outputdata_dir = './resu_fmridata_preprocessed/'
sink(file = file(paste0(outputdata_dir, "output2.txt")))



## Iterate over all subjects
for (ids in 1:155) {
#ids = 100  
  
  print(paste('** Subject', sprintf('%03d',ids), '**'))
  
  
  #***********************************************************#
  #            Read and get data                              #
  #***********************************************************#
  strsubj <- sprintf('%03d',ids)
  anatfname <- paste0(fmridata_dir,'sub-pixar', 
                      strsubj, '/anat/sub-pixar', strsubj, '_T1w.nii.gz')
  funcfname <- paste0(fmridata_dir,'/sub-pixar', 
                      strsubj, '/func/sub-pixar', strsubj, '_task-pixar_bold.nii.gz')
  mydat_oro = readNIfTI(funcfname, reorient = FALSE)
  tr <- 2
  DROP <- 5 #10 seconds for stabilization
  times <- (DROP+1):dim(mydat_oro)[4]
  run_fmri = copyNIfTIHeader(mydat_oro, mydat_oro[,,,times], drop = TRUE)
  rm(mydat_oro)
  
  
  #************************************************************#
  #           Preprocessing: 
  # 1. Image Realignment 
  #************************************************************#
  # realign_batch = build_spm12_realign( 
  #   filename = run_fmri, 
  #   register_to = "mean",
  #   reslice = "mean"
  # )
  # 
  # print(names(realign_batch))
  # print(names(unlist(realign_batch$spm)))

  if (have_matlab()) {
    realigned = spm12_realign(
      filename = run_fmri,
      register_to = "mean",
      reslice = "mean",
      clean = FALSE
    )
    #print(names(realigned))
    ##[1] "spm"      "script"   "outfiles" "rp"       "mean"     "mat"     "result"
    #print(class(realigned$outfiles))
    ##[1] "C:/Users/namgi/AppData/Local/Temp/RtmpSurVwV/file1c90686f3583.nii"
  }
  
  #************************************************************#
  #           Preprocessing: 
  # 2. Slice Time Correction
  #************************************************************#
  nslices = oro.nifti::nsli(run_fmri)
  slice_order = 1:nslices
  ref_slice = slice_order[median(seq(nslices))] ##Middel slice
  ta = tr - tr/nslices
  n_time_points = ntim(run_fmri)
  
  # # Using run_fmri because if you do not have
  # # matlab, then realigned not available
  # st_batch = build_spm12_slice_timing(
  #   filename = run_fmri,
  #   time_points = seq(n_time_points),
  #   nslices = nslices,
  #   tr = tr,
  #   ref_slice = ref_slice,
  #   prefix = "a")
  # #print(names(st_batch))
  # # [1] "spm"                 "script"              "orig_filename"      
  # # [4] "base_name"           "temporary_directory" "outfile" 
  # 
  
  if (have_matlab()) {
    
    st_results = spm12_slice_timing(
      filename = realigned[['outfiles']],  ##use realigned
      time_points = seq(n_time_points),
      nslices = nslices,
      tr = tr, 
      slice_order = slice_order,
      ta = ta, 
      ref_slice = ref_slice,
      prefix = "a", 
      clean = FALSE, 
      retimg = FALSE)
    aimg = st_results$outfile
    #print(aimg)  ##"C:\\Users\\...."  ; x = readnii(aimg); 64 64 32 163
    #x = readnii(aimg); print(dim(x)); oro.nifti::orthographic(x)
    mean_img = realigned[["mean"]]
    #mean_nifti = readnii(mean_img)##dim: 64 64 32
    
  }
  
  
  #************************************************************#
  #           Preprocessing: 
  # 3. Spatial Normalization
  #   a) Direct Normalization
  #************************************************************#
  ## Direct normalization to MNI
  # if (have_matlab()) {
  #   bbox = matrix(
  #     c(-90, -126, -72, 
  #       90, 90, 108), 
  #     nrow = 2, byrow = TRUE)
  #   direct_norm = spm12_normalize(
  #     filename = mean_img,
  #     other.files = c(mean_img, aimg),
  #     bounding_box = bbox,
  #     clean = FALSE
  #   )
  #   # print(names(direct_norm)) ##[1] "spm"      "script"   "outfiles" "result"  
  #   dnorm_files = direct_norm$outfiles
  #   dnorm_mean_img = readnii(dnorm_files[1])  #To compare with templage
  #   #oro.nifti::orthographic(dnorm_mean_img) #skip visualization
  #   #x = readnii(readnii(dnorm_files[2])); print(dim(x)); #Error: cannot allocate vector of size 1.1 Gb
  # }
  
  
  #************************************************************#
  #           Preprocessing: 
  # 3. Spatial Normalization
  #   b) Indirect Normalization (==co-registration)
  #     Normalize anatomical image to the mean image
  #************************************************************#
  if (have_matlab()) {
    anatomical = anatfname   #files["anatomical"]
    anat_img = checknii(anatomical)
    #print(anat_img)
    acpc_reorient(
      infiles = anat_img,
      modality = "T1")

    coreg = spm12_coregister(
      fixed = mean_img,
      moving = anat_img,
      prefix = "r")

    coreg_anat = coreg$outfile
    #coreg_img = readnii(coreg_anat)
    #double_ortho(coreg_img, mean_nifti) #skip visualization??
  }
  
  
  #************************************************************#
  #   Anatomical MR Segmentation 
  #************************************************************#
  if (have_matlab()) {
    seg_res = spm12_segment(
      filename = coreg_anat,
      set_origin = FALSE,
      retimg = FALSE,
      clean = FALSE)
    #print(names(seg_res)) #[1]  "spm"  "script"  "result"
    #    "outfiles"            "outmat"              "deformation"
    #    "inverse_deformation"
  }


  ## Hard Segmentation:
  alpha = function(col, alpha = 1) {
    cols = t(col2rgb(col, alpha = FALSE)/255)
    rgb(cols, alpha = alpha)
  }
  # if (have_matlab()) {
  #   seg_files = check_nifti(seg_res$outfiles)
  #   hard_seg = spm_probs_to_seg(seg_files)
  #   hard_seg[ hard_seg > 3] = 0
  # 
  #   ortho2(coreg_img, hard_seg,
  #          col.y = alpha(c("red", "green", "blue"), 0.5))
  # }
  
  
  ## Applying Spatial Normalization Transform
  bbox = matrix(
    c(-90, -126, -72,
      90, 90, 108),
    nrow = 2, byrow = TRUE)
  if (have_matlab()) {
    norm = spm12_normalize_write(
      deformation = seg_res$deformation,
      other.files = c(coreg_anat, mean_img, aimg),
      bounding_box = bbox,
      retimg = FALSE)
    #print(names(norm))
    # [1] "spm"          "script"       "deformation"  "other_fnames"
    # [5] "outfiles"     "result"
    norm_data = norm$outfiles
    names(norm_data) = c("anat", "mean", "fmri")
    # norm_mean_img = readnii(norm_data["mean"])
    # norm_anat_img = readnii(norm_data["anat"])
  }
  # #Now we have the indirect spatially normalized data in MNI template space
  
  
  ## Comparison of Direct and Indirect Normalization
  # if (have_matlab()) {
  #   template_path = file.path(spm_dir(), 
  #                             "canonical", "avg152T1.nii")
  #   
  #   template = readnii(template_path)
  #   
  #   dnorm_mask = dnorm_mean_img > quantile(
  #     dnorm_mean_img[dnorm_mean_img > 0],
  #     probs = 0.6)
  #   norm_mask = norm_mean_img > quantile(
  #    norm_mean_img[norm_mean_img > 0],
  #    probs = 0.6)
  # 
  #   double_ortho(template, norm_anat_img)
  #   double_ortho(template, norm_mean_img)
  #   double_ortho(norm_mean_img, norm_anat_img)
  #   ortho2(template, norm_mask, col.y = alpha("red", 0.5))
  # 
  #   double_ortho(template, dnorm_mean_img)
  #   ortho2(template, dnorm_mask, col.y = alpha("red", 0.5))
  #   ##  UP: dnorm=direct norm is fitting well !!
  #   
  #   ##norm vs dnorm
  #   double_ortho(norm_mean_img, dnorm_mean_img)
  # }
  
  
  #************************************************************#
  #   Spatial Smoothing
  #  WE WILL AVOID SPATIAL SMOOTHING 
  # ref: https://arxiv.org/pdf/1705.02141.pdf
  #************************************************************#
  # if (have_matlab()) {
  #   smoothed = spm12_smooth(
  #     filename = norm_data["fmri"], ##4D fMRI data
  #     fwhm = 8,
  #     prefix = "s",
  #     retimg = FALSE
  #   )
  #   smoothed_data = smoothed$outfiles
  #   #In many applications, this is the data you will 
  #   #use for post-processing and analysis. Motion
  #   #correction has usually been applied above, but 
  #   #some motion correct this data as well.
  # }
  
  # if (have_matlab()) {
  #   smoothed_mean = spm12_smooth(
  #     filename = norm_data["mean"],  ##3D mean data
  #     prefix = "s",
  #     retimg = FALSE
  #   )
  #   smoothed_mean_data = smoothed_mean$outfiles
  #   #Here we can smooth the mean image in MNI space.
  #   
  #   #smooth_mean_img = readnii(smoothed_mean_data)
  #   #ortho2(smooth_mean_img)
  # }
  
  ## Copy the smoothed fMRI data ##
  
  source_fname = norm_data["fmri"] #not smoothed #moothed == smoothed_data
  target_fname = paste0(outputdata_dir, 'swasub-pixar', strsubj, '_task-pixar_bold.nii')
  x = file.copy(source_fname, target_fname)
  if(isTRUE(x)) {
    file.remove(source_fname)
  } else {
    warning(paste('Subject', strsubj,': Final nii file could not be copied.'))
  }
  
  
  rm(run_fmri)
  
  
}


sink()