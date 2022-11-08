#!/bin/bash
export G=5
export source_folder='paper/results/'
export base_folder_beginning='oxy_none_NA_'
export base_folder_end='_13_12_15_15/'
export m_pred_predictors='5'
export m_pred_response='5'
export grid_size='1'
export date=2020-01-15
export new_response_pcs2='11'
export MC_max_clusters=100
export MC_max_pred=50
export out_file='test_out.out'

# The application(s) to execute along with its input arguments and options:
R CMD BATCH --vanilla paper/code/data_analysis_results/argo_predict_results_new.R $out_file
  
