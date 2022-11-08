#!/bin/bash
#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=simu
#SBATCH --mail-user=dyarger@umich.edu
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=7000m 
#SBATCH --time=00:25:00
#SBATCH --account=stats_dept2
#SBATCH --partition=standard
#SBATCH --output=/home/%u/jobmaterials/%x-%j.log
#SBATCH --array=1-1000

# The application(s) to execute along with its input arguments and options:
module load R/4.1.0
module load Rtidyverse/4.1.0
module load gsl
module load jags
cd /home/dyarger/argo_functional_regression/
export out_file=/scratch/stats_dept_root/stats_dept1/shared_data/argo/simulation/simu_$SLURM_ARRAY_TASK_ID.out

R CMD BATCH --vanilla paper/code/simulation/simu_one_var.R $out_file


