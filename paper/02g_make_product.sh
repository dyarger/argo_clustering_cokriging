#!/bin/sh
# computing info
#SBATCH --job-name=product
#SBATCH --mail-user=dyarger@umich.edu
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=14000m 
#SBATCH --time=7:00:00
#SBATCH --account=stats_dept1
#SBATCH --partition=standard
#SBATCH --output=/home/%u/jobmaterials/%x-%j.log

export G=5
export source_folder='/scratch/stats_dept_root/stats_dept1/shared_data/argo/'
export base_folder_beginning='oxy_none_NA_'
export base_folder_end='_13_12_15_15/'
export m_pred_predictors='5'
export m_pred_response='5'
export grid_size='1'
export new_response_pcs2='11'
#export out_file='test_out.out'
export out_file=/scratch/stats_dept_root/stats_dept1/shared_data/argo/product.out

module load R/4.1.0
module load Rtidyverse/4.1.0
module load Rgeospatial
module load gsl
cd /home/dyarger/argo_functional_regression/
R_LIBS_SITE="/sw/arcts/centos7/stacks/gcc/8.2.0/Rtidyverse/4.1.0/2021-06-28"
# The application(s) to execute along with its input arguments and options:
/sw/arcts/centos7/stacks/gcc/8.2.0/R/4.1.0/bin/R CMD BATCH --vanilla paper/code/data_analysis_results/argo_make_product.R $out_file
  
