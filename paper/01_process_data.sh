#!/bin/bash
export soccom_date='05_05_21'
export variable='oxy'
export out_file='test_out.out'
export out_file2='test_out2.out'
export out_file3='test_out3.out'
export out_file4='test_out4.out'

# subset data for oxygen and nitrate and core
# The application(s) to execute along with its input arguments and options:
R CMD BATCH --vanilla paper/code/data_preparation/subset_soccom_data.R $out_file

export variable='nitrate'
R CMD BATCH --vanilla paper/code/data_preparation/subset_soccom_data.R $out_file2

export argo_date='06_21'
R CMD BATCH --vanilla paper/code/data_preparation/subset_core_data.R $out_file3
# Make data plots for paper
R CMD BATCH --vanilla paper/code/data_preparation/data_plots.R $out_file4

  