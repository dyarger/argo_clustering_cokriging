#!/bin/sh
 # computing info
export out_file=test_out.out

# The application(s) to execute along with its input arguments and options:
R CMD BATCH paper/code/data_analysis_results/argo_AIC_results.R $out_file
