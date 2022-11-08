#!/bin/bash

# set locations of the Roemmich and Gilson product (for fronts)
export RG_location='~/Downloads/RG_ArgoClim_33pfit_2019_mean.nc'
# SOCCOM data
export soccom_location='~/Downloads/SOCCOM_HiResQC_MLR_05May2021_netcdf/'
export soccom_date='05_05_21'

# location of Argo data mirror
export argo_location='~/Downloads/202106-ArgoData/'
export out_file='test_out.out'
export out_file2='test_out2.out'
export out_file3='test_out3.out'

R CMD BATCH --vanilla paper/code/data_preparation/get_fronts.R $out_file
R CMD BATCH --vanilla paper/code/data_preparation/get_soccom_data.R $out_file2
R CMD BATCH --vanilla paper/code/data_preparation/get_core_data.R $out_file3
  
  
  
  
  