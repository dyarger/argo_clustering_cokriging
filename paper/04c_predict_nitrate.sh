#!/bin/sh
# computing info
#SBATCH --job-name=predict_nitrate
#SBATCH --mail-user=dyarger@umich.edu
#SBATCH --mail-type=END,ARRAY_TASKS,ERROR
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100000m 
#SBATCH --time=40:00:00
#SBATCH --account=stats_dept1
#SBATCH --partition=largemem
#SBATCH --output=/home/%u/jobmaterials/%x-%j.log
#SBATCH --array=0

export m_pred_predictors='5'
export m_pred_response='5'
export grid_size='1'
export new_response_pcs2='10'
export dates=(2020-01-15)
export MC_max='50'
export date=${dates[$SLURM_ARRAY_TASK_ID]}
export source_folder='/scratch/stats_dept_root/stats_dept1/shared_data/argo/'
export base_folder='nitrate_none_NA_5_12_12_15_15/'
export out_file=/scratch/stats_dept_root/stats_dept1/shared_data/argo/nitrate_pred_TSO_$SLURM_ARRAY_TASK_ID.out

# The application(s) to execute along with its input arguments and options:
module load R/4.1.0
module load Rtidyverse/4.1.0
module load Rgeospatial
module load gsl
cd /home/dyarger/argo_functional_regression/
R_LIBS_SITE="/sw/arcts/centos7/stacks/gcc/8.2.0/Rtidyverse/4.1.0/2021-06-28"
/sw/arcts/centos7/stacks/gcc/8.2.0/R/4.1.0/bin/R CMD BATCH --vanilla paper/code/data_analysis/argo_predict.R $out_file

