#!/bin/bash
#“#SBATCH” directives that convey submission options:

export out_file='test_out.out'
# The application(s) to execute along with its input arguments and options:
R CMD BATCH --vanilla code/simulation/simu_one_var_results.R $out_file