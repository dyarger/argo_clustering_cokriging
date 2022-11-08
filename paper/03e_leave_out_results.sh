#!/bin/bash
export variable='oxy'
export base_folder_end='_5_13_11_15_15/'
export source_folder='paper/results/'
export out_file='test_out.out'

# The application(s) to execute along with its input arguments and options:
R CMD BATCH --vanilla paper/code/data_analysis_results/argo_leave_out_results.R $out_file
  
