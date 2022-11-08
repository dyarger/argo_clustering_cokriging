#!/bin/sh
 # computing info
#SBATCH --job-name=floats
#SBATCH --mail-user=dyarger@umich.edu
#SBATCH --mail-type=END,ARRAY_TASKS
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=82000m 
#SBATCH --time=24:00:00
#SBATCH --account=stats_dept1
#SBATCH --partition=largemem
#SBATCH --output=/home/%u/jobmaterials/%x-%j.log
#SBATCH --array=1-10

# number of clusters and pcs
export G='5'
export pc_predictors='13'
export pc_response='11'

# variable and whether to add any core profiles/data subset
export var='oxy'
export n_core_profiles='0'
export leave_out_type='floats'

# number of EM and MC iterations
export EM_iter='10'
export EM_iter_init='12'
export MC_max='15'
export MC_max_AIC='3'

export n_folds=10
export fold=$SLURM_ARRAY_TASK_ID

# number of vecchia neighbors for training and model selection, prediction
export m_train_predictors='15'
export m_train_response='15'
export m_AIC_predictors='3'
export m_AIC_response='3'
export m_pred_predictors='35'
export m_pred_response='35'
export maxit='20' # number of optim iterations for updating spatial parameters

# MRF info
export nn_strategy='nn'
export nn_type='spacewarp_day'
export nn_time_param='12'
export nn_lat_param='3'
export nn='15'

export out_file=/scratch/stats_dept_root/stats_dept1/shared_data/argo/float_$SLURM_ARRAY_TASK_ID.out

# The application(s) to execute along with its input arguments and options:
module load R/4.1.0
module load Rtidyverse/4.1.0
module load Rgeospatial
module load gsl
cd /home/dyarger/argo_functional_regression/
R_LIBS_SITE="/sw/arcts/centos7/stacks/gcc/8.2.0/Rtidyverse/4.1.0/2021-06-28"
/sw/arcts/centos7/stacks/gcc/8.2.0/R/4.1.0/bin/R CMD BATCH --vanilla paper/code/data_analysis/argo_estimate_model.R $out_file 
