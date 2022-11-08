#!/bin/sh
# computing info
#SBATCH --job-name=predict
#SBATCH --mail-user=dyarger@umich.edu
#SBATCH --mail-type=END,ARRAY_TASKS,ERROR
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100000m 
#SBATCH --time=50:00:00
#SBATCH --account=stats_dept1
#SBATCH --partition=largemem
#SBATCH --output=/home/%u/jobmaterials/%x-%j.log
#SBATCH --array=0-83

export source_folder='/scratch/stats_dept_root/stats_dept1/shared_data/argo/'
export base_folder='oxy_none_NA_5_13_12_15_15/'
export out_file=/scratch/stats_dept_root/stats_dept1/shared_data/argo/pred_$SLURM_ARRAY_TASK_ID.out
export m_pred_predictors='5'
export m_pred_response='5'
export new_response_pcs2='11'
export MC_max='50'
export grid_size='1'
export dates=(2020-01-15 2020-02-15 2020-03-15 2020-04-15 2020-05-15 2020-06-15 2020-07-15 2020-08-15 2020-09-15 2020-10-15 2020-11-15 2020-12-15 2014-07-15 2014-08-15 2014-09-15 2014-10-15 2014-11-15 2014-12-15 2015-01-15 2015-02-15 2015-03-15 2015-04-15 2015-05-15 2015-06-15 2015-07-15 2015-08-15 2015-09-15 2015-10-15 2015-11-15 2015-12-15 2016-01-15 2016-02-15 2016-03-15 2016-04-15 2016-05-15 2016-06-15 2016-07-15 2016-08-15 2016-09-15 2016-10-15 2016-11-15 2016-12-15 2017-01-15 2017-02-15 2017-03-15 2017-04-15 2017-05-15 2017-06-15 2017-07-15 2017-08-15 2017-09-15 2017-10-15 2017-11-15 2017-12-15 2018-01-15 2018-02-15 2018-03-15 2018-04-15 2018-05-15 2018-06-15 2018-07-15 2018-08-15 2018-09-15 2018-10-15 2018-11-15 2018-12-15 2019-01-15 2019-02-15 2019-03-15 2019-04-15 2019-05-15 2019-06-15 2019-07-15 2019-08-15 2019-09-15 2019-10-15 2019-11-15 2019-12-15 2021-01-15 2021-02-15 2021-03-15 2021-04-15 2021-05-15 2021-06-15)
export date=${dates[$SLURM_ARRAY_TASK_ID]}

# The application(s) to execute along with its input arguments and options:
module load R/4.1.0
module load Rtidyverse/4.1.0
module load Rgeospatial
module load gsl
cd /home/dyarger/argo_functional_regression/
R_LIBS_SITE="/sw/arcts/centos7/stacks/gcc/8.2.0/Rtidyverse/4.1.0/2021-06-28"
/sw/arcts/centos7/stacks/gcc/8.2.0/R/4.1.0/bin/R CMD BATCH --vanilla paper/code/data_analysis/argo_predict.R $out_file

