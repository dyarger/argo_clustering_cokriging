#!/bin/bash
export G=5
export source_folder='paper/results/'
export base_folder_beginning='oxy_none_NA_'
export base_folder_end='_13_12_15_15/'
export out_file='test_out.out'

# The application(s) to execute along with its input arguments and options:
R CMD BATCH --vanilla paper/code/data_analysis_results/argo_estimate_model_results.R $out_file
  