#!/bin/sh
 # computing info
#SBATCH --job-name=nitrate
#SBATCH --mail-user=dyarger@umich.edu
#SBATCH --mail-type=END,ARRAY_TASKS
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=41000m 
#SBATCH --time=80:00:00
#SBATCH --account=stats_dept1
#SBATCH --partition=largemem
#SBATCH --output=/home/%u/jobmaterials/%x-%j.log

# number of clusters and pcs
export G='5'
export pc_predictors='12'
export pc_response='12'

# variable and whether to add any core profiles/data subset
export var='nitrate'
export n_core_profiles='0'
export leave_out_type='none'

# number of EM and MC iterations
export EM_iter='40'
export EM_iter_init='12'
export MC_max='15'
export MC_max_AIC='3'

# number of vecchia neighbors for training and model selection
export m_train_predictors='15'
export m_train_response='15'
export m_AIC_predictors='3'
export m_AIC_response='3'
export maxit='20' # number of optim iterations for updating spatial parameters

# MRF info
export nn_strategy='nn'
export nn_type='spacewarp_day'
export nn_time_param='12'
export nn_lat_param='3'
export nn='15'

export out_file=/scratch/stats_dept_root/stats_dept1/shared_data/argo/nitrate_full_$SLURM_ARRAY_TASK_ID.out


# The application(s) to execute along with its input arguments and options:
cd /home/dyarger/argo_functional_regression/
module load R/4.1.0
module load Rtidyverse/4.1.0
module load Rgeospatial
module load gsl
module load gcc/8.2.0
R_LIBS_SITE="/sw/arcts/centos7/stacks/gcc/8.2.0/Rtidyverse/4.1.0/2021-06-28"
/sw/arcts/centos7/stacks/gcc/8.2.0/R/4.1.0/bin/R CMD BATCH --vanilla paper/code/data_analysis/argo_estimate_model.R $out_file



