#!/bin/sh
 # computing info
#SBATCH --job-name=AIC
#SBATCH --mail-user=dyarger@umich.edu
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=41000m 
#SBATCH --time=0:45:00
#SBATCH --account=stats_dept1
#SBATCH --partition=standard
#SBATCH --output=/home/%u/jobmaterials/%x-%j.log
#SBATCH --array=1-56

export G='5'
export var='oxy'
export pc_response='12'
export m_train_predictors='15'
export m_train_response='15'
export m_AIC_predictors='15'
export m_AIC_response='15'
export MC_max_AIC='50'
export array_id=$SLURM_ARRAY_TASK_ID
export base_folder='/scratch/stats_dept_root/stats_dept1/shared_data/argo/'
export out_file=/scratch/stats_dept_root/stats_dept1/shared_data/argo/AIC_TSO_$SLURM_ARRAY_TASK_ID.out

# The application(s) to execute along with its input arguments and options:
cd /home/dyarger/argo_functional_regression/
module load R/4.1.0
module load Rtidyverse/4.1.0
module load Rgeospatial
module load gsl
module load gcc/8.2.0
R_LIBS_SITE="/sw/arcts/centos7/stacks/gcc/8.2.0/Rtidyverse/4.1.0/2021-06-28"
/sw/arcts/centos7/stacks/gcc/8.2.0/R/4.1.0/bin/R CMD BATCH paper/code/data_analysis/argo_AIC.R $out_file
