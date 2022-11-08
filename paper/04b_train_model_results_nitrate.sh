#!/bin/bash
export source_folder='paper/results/'
export G=5
export base_folder_beginning='nitrate_none_NA_'
export base_folder_end='_12_12_15_15/'
export out_file='test_out.out'
R CMD BATCH --vanilla paper/code/data_analysis_results/argo_estimate_model_results.R $out_file

  