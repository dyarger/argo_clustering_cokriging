Sys.setenv(R_LIBS_SITE = '/sw/arcts/centos7/stacks/gcc/8.2.0/Rtidyverse/4.1.0/2021-06-28')
library(fstmr)
library(dplyr)
# load data for which we want to predict core memberships
source_folder <-Sys.getenv('source_folder')
base_folder <-Sys.getenv('base_folder')
# source_folder='paper/results/'
# base_folder='oxy_none_NA_5_13_12_15_15/'
load(paste0(source_folder, base_folder, 'final_results.RData'))
new_number <- as.integer(Sys.getenv('new_response_pcs2'))
estimated_model[['parameters']][['Omega_response2']] <- 
  lapply(estimated_model[['parameters']][['Omega_response2']], function(x) {
    x[,1:new_number]
  })
estimated_model[['parameters']][['variances_response']] <- 
  lapply(estimated_model[['parameters']][['variances_response']], function(x) {
    x[1:new_number]
  })
estimated_model[['parameters']][['range_params_r']] <- 
  lapply(estimated_model[['parameters']][['range_params_r']], function(x) {
    x[,1:new_number]
  })
estimated_model[['params']][['pc_response2']] <- new_number
estimated_model[['params']][['m_pred_response']] <- as.integer(Sys.getenv('m_pred_response'))
estimated_model[['params']][['m_pred_predictors']] <- as.integer(Sys.getenv('m_pred_predictors'))
date_to_predict <- as.Date(Sys.getenv('date'), format = '%Y-%m-%d')

MC_max <- as.integer(Sys.getenv('MC_max'))
estimated_model[['params']][['MC_max']] <- MC_max
grid_size <- as.integer(Sys.getenv('grid_size'))

params <- estimated_model[['params']]                  
parameters <- estimated_model[['parameters']]
data_inputs <- estimated_model[['data_inputs']]
G <- params[['G']]

# load core data
load('paper/data/core_processed_06_21.RData')
#core_data_subset <- core_data[1:45000,]
core_data_subset <- core_data
rm(core_data)
core_data_subset[['day']] <- core_data_subset[['date']]
core_data_subset[['dayofyear']] <- julian(as.Date(core_data_subset[['date']]), origin = as.Date('2000-01-01')) %% 365.25

not_dup <- !duplicated(core_data_subset[['profile_unique']])
core_data_subset_locations <- core_data_subset[not_dup,]
locations_combined <- rbind(data_inputs[['locs']][,c('longitude', 'latitude')],
                            core_data_subset_locations[,c('longitude', 'latitude')])

# create grid
source('paper/code/src/pred_at_pressure_levels_src.R')
grid_locs <- get_grid_info(grid_size)
grid_locs_expanded <- cbind(grid_locs, profile_unique = 1:nrow(grid_locs),
                            time = as.vector(julian(as.Date(date_to_predict, format = "%Y-%m-%d"), origin = as.Date("2000-01-01"))))
grid_locs_expanded[['dayofyear']] <- grid_locs_expanded[['time']] %% 365.25
grid_locs_expanded[['day']] <- grid_locs_expanded[['time']] %% 365.25

rdist_comparison <- spam::spam_rdist.earth(as.matrix(grid_locs_expanded[,c('longitude', 'latitude')]), 
                                           as.matrix(locations_combined), 
                                           miles = F, delta = 325*360/(6378.388*2*pi))
rdist_comparison_mat <- spam::as.dgCMatrix.spam(rdist_comparison)
nn_within <- table(rdist_comparison_mat@i)
not_in <- 1:nrow(rdist_comparison_mat) %in% (as.numeric(names(nn_within)) + 1)
grid_use <- grid_locs_expanded[not_in, ] %>% 
  filter(depth > 600, latitude > -75) %>%
  mutate(longitude = ifelse(longitude > 180, longitude - 360, longitude))

Sys.time()
# begin to set things up
pred_set_up <- fstmr:::set_up_prediction(estimated_model, params,core_data = core_data_subset,
                                         grid = grid_use)
Sys.time()
MC_max_cluster <- max(c(100, MC_max))
cluster_mat_full <- fstmr:::rPotts_prediction(params[['G']], pred_set_up[['neighbors_all']], pred_set_up[['has_data']],
                                              pred_set_up[['p_mat']], pred_set_up[['like_mat']],
                                              estimated_model[['parameters']][['theta']], 
                                              n_samples = MC_max_cluster, skip = 4, init_pred_cluster = pred_set_up[['memberships']])

# save a bunch of cluster predictions
cluster_mat_pred <- cluster_mat_full[pred_set_up[['is_prediction']], ]
save(cluster_mat_pred, grid_use,
     file =  paste0(source_folder, base_folder, 'predict_clusters_',
                    date_to_predict, '_', grid_size, '_', MC_max_cluster, '.RData'))

# predict
cluster_mat <- cluster_mat_full[,1:MC_max]
Sys.time()
pred_score_results <- fstmr:::pred_E_step(pred_set_up, params, cluster_mat, 
                                          prediction_only = T)
Sys.time()

data_inputs <- pred_set_up[['prediction_model']][['data_inputs']]
unchanged_file <- paste0(source_folder, base_folder, 'data_inputs_pred_',
                         grid_size, '.RData')
if (!file.exists(unchanged_file)) {
  save(data_inputs, file = unchanged_file)
}

pred_set_up[['prediction_model']][['data_inputs']] <- NULL
save(pred_set_up, cluster_mat, pred_score_results, grid_use, 
     file = paste0(source_folder, base_folder, 'predict_IS_',
                   date_to_predict, '_', grid_size, '_', params[['m_pred_predictors']], '_',
                   params[['m_pred_response']], '_', MC_max, '.RData'))
