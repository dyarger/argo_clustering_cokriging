Sys.setenv(R_LIBS_SITE = '/sw/arcts/centos7/stacks/gcc/8.2.0/Rtidyverse/4.1.0/2021-06-28')
library(fstmr)
library(dplyr)
library(fda)
library(Matrix)
# variables and basic setup
set.seed(56)
params = list()
params[['results_folder']] <- "/scratch/stats_dept_root/stats_dept1/shared_data/argo/"

# type of leave out data
params[['leave_out_type']] = Sys.getenv('leave_out_type')
params[['n_folds']] = as.integer(Sys.getenv('n_folds'))
params[['folds']] <-  as.integer(Sys.getenv('fold'))
params[['leave_out_prediction_core']] = F

# set up data variables
params[['n_core_profiles']] <- ifelse(Sys.getenv('n_core_profiles') == 'all', 
                                      'all', as.integer(Sys.getenv('n_core_profiles')))
params[['id']] = 'profile_unique'
params[['ind_name']] = 'pressure'
params[['loc_names']] = c('longitude', 'latitude')
params[['var']] = params[['Y']] = Sys.getenv('var')
params[['X']] = c('temp', 'psal')
params[['domain_Y']] = c(0,2000)
params[['domain_X']] <- list(c(0,2000), c(0,2000))

# initialization
params[['init_strategy']] = 'kmeans'
params[['levels']] = list('nitrate' = seq(0,2000,length.out = 15),
                          'oxy' = seq(0,2000,length.out = 15),
                          'temp'= seq(0,2000,length.out = 15),
                          'psal'= seq(0,2000,length.out = 15))
# MRF parameters
params[['nn_strategy']] = Sys.getenv('nn_strategy')
params[['nn_type']] = Sys.getenv('nn_type')
params[['nn_range_time']] = as.numeric(Sys.getenv('nn_range_time'))
params[['nn_time_param']] = as.numeric(Sys.getenv('nn_time_param'))
params[['nn_lat_param']] = as.numeric(Sys.getenv('nn_lat_param'))
params[['nn']] = as.integer(Sys.getenv('nn'))
params[['remove_between_land']] <- T

# basis parameters
params[['knots']] = c(seq(0,100,5), seq(110, 290, by = 10), seq(300,1000,50), 1100, 1200, 1300, 1400, 1500, 1600, 2000)
params[['basis_response']] = create.bspline.basis(params[['knots']])
params[['basis_predictors']] = list(create.bspline.basis(params[['knots']]),
                                    create.bspline.basis(params[['knots']])) 

# 4 spatial covariance parameters
params[['m_train_predictors']] <- as.integer(Sys.getenv('m_train_predictors'))
params[['m_train_response']] <- as.integer(Sys.getenv('m_train_response'))
params[['m_AIC_predictors']] <- as.integer(Sys.getenv('m_AIC_predictors'))
params[['m_AIC_response']] <- as.integer(Sys.getenv('m_AIC_response'))
params[['m_pred_predictors']] <- as.integer(Sys.getenv('m_pred_predictors'))
params[['m_pred_response']] <-  as.integer(Sys.getenv('m_pred_response'))
params[['maxit']] = as.integer(Sys.getenv('maxit'))
params[['time']] <- 'day'
params[['dayofyear']] <- 'dayofyear'
params[['covariance_function']] <- 'exponential_spheretime_warp'

# number of clusters and pcs
params[['G']] = as.integer(Sys.getenv('G'))
params[['pc_response1']] <- params[['pc_predictors']] <- as.integer(Sys.getenv('pc_predictors'))
params[['pc_response2']] <- as.integer(Sys.getenv('pc_response'))

# smoothing parameters
params[['lambda_mean_response']] <- 10^(3:5)
params[['lambda_mean_predictors']] <- 10^(3:5)
params[['lambda_pcs_response']] <- 10^(7:10)
params[['lambda_lt']] <- 10^(5:8)
params[['lambda_pcs_predictors']] <- 10^(4:7)
params[['cv_skip']] <- 10


# number of MCEM and EM iterations
params[['MC_max']] =  as.integer(Sys.getenv('MC_max'))
params[['MC_max_AIC']] =  as.integer(Sys.getenv('MC_max_AIC'))
params[['EM_iter']] =  as.integer(Sys.getenv('EM_iter'))
params[['EM_iter_init']] <-  as.integer(Sys.getenv('EM_iter_init'))

##### Load and subset data #####
load(paste0('paper/data/soccom_processed_', params[['Y']], '_05_05_21.RData'))
source('paper/code/src/prepare_argo_data_src.R')
data_leave_out <- leave_out_data(df_list, params[['leave_out_type']], 
                                 params[['folds']], params[['n_folds']], params[['leave_out_prediction_core']],
                                 params[['n_core_profiles']], params[['core_sampling_strategy']],
                                 reduce_data = F
)
df_list <- data_leave_out[[1]]; df_list_pred <- data_leave_out[[2]]
held_out_profs <- df_list_pred[[1]]
rm(data_leave_out)

# things to print each MC iteration
params[['dir_name']] <- paste0(params[['results_folder']], params[['Y']], '_', 
                               params[['leave_out_type']], '_', params[['folds']], '_', params[['G']], '_',
                               params[['pc_predictors']], '_', params[['pc_response2']], '_',
                               params[['m_train_predictors']], '_', params[['m_train_response']], '/')
system(paste('mkdir', params[['dir_name']]))

printf <- function(model, diagnostics, parameters, data_inputs, params) {
  print(parameters[['measurement_error_response']])
  print(c(parameters[['lambda_mean_response']], range(params[['lambda_mean_response']])))
  print(c(parameters[['lambda_mean_predictors']], range(params[['lambda_mean_predictors']])))
  print(c(parameters[['lambda_pcs_response']], range(params[['lambda_pcs_response']])))
  print(c(parameters[['lambda_pcs_predictors']], range(params[['lambda_pcs_predictors']])))
  print(c(parameters[['lambda_lt']], range(params[['lambda_lt']])))
  if (params[['leave_out_type']] == 'none') {
    print(diagnostics[['AIC']])
    iteration <- sum(diagnostics[['measurement_error_response']] != 0)
    save(parameters, file = paste0(params[['dir_name']], iteration, '_parameters.RData'))
    gc()
  }
}

# estimate model
estimated_model <- fstmr::fstmr(data = df_list, params,
                                compute_diagnostics = ifelse(params[['leave_out_type']] == 'none', T, F),
                                verbose = 50, printf = printf)

save(params, estimated_model, diagnostics,
     file = paste0(params[['dir_name']], 'final_results.RData'))

# Prediction for leave out
if (params[['leave_out_type']] != 'none') {
  print(Sys.time())
  pred_set_up <- fstmr:::set_up_prediction(estimated_model, params)
  print(Sys.time())
  cluster_mat <- fstmr:::rPotts_prediction(params[['G']], pred_set_up[['neighbors_all']], pred_set_up[['has_data']],
                                           pred_set_up[['p_mat']], pred_set_up[['like_mat']],
                                           estimated_model[['parameters']][['theta']], 
                                           n_samples = params[['MC_max']], skip = 4)
  print(summary(apply(cluster_mat, 1, var)))
  print(Sys.time())
  
  pred_score_results <- fstmr:::pred_E_step(pred_set_up, params, cluster_mat)
  
  print(Sys.time())
  pred_at_p_levels <- fstmr:::pressure_level_prediction(held_out_profs, params, pred_set_up, pred_score_results,
                                                        cluster_mat)
  
  pred_set_up[['prediction_model']][['data_inputs']] <- NULL
  save(pred_set_up, pred_at_p_levels, pred_score_results, held_out_profs, cluster_mat,
       file = paste0(params[['dir_name']], 'prediction.RData'))
}