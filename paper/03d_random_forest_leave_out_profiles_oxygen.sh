#!/bin/sh
 # computing info
#SBATCH --job-name=rf_profiles
#SBATCH --mail-user=dyarger@umich.edu
#SBATCH --mail-type=END,ARRAY_TASKS
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=14000m 
#SBATCH --time=00:20:00
#SBATCH --account=stats_dept1
#SBATCH --partition=standard
#SBATCH --output=/home/%u/jobmaterials/%x-%j.log
#SBATCH --array=1-5

export leave_out_type='none'
export leave_out_type='profiles'
export n_folds='5'
export soccom_date='05_05_21'
export var='oxy'
export fold=$SLURM_ARRAY_TASK_ID
export out_file=/scratch/stats_dept_root/stats_dept1/shared_data/argo/rf_profiles_$fold_TSO.out  

# The application(s) to execute along with its input arguments and options:
cd /home/dyarger/argo_functional_regression/
module load R/4.1.0
module load Rtidyverse/4.1.0
module load Rgeospatial
module load gsl
module load gcc/8.2.0
R_LIBS_SITE="/sw/arcts/centos7/stacks/gcc/8.2.0/Rtidyverse/4.1.0/2021-06-28"
/sw/arcts/centos7/stacks/gcc/8.2.0/R/4.1.0/bin/R CMD BATCH --vanilla paper/code/data_analysis/argo_random_forest.R $out_file