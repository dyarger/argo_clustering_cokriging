library(fda)
library(ncdf4)
library(spam)
.libPaths("/Library/Frameworks/R.framework/Versions/4.1/Resources/miworkspace-library")
library(tidyverse)
source('paper/code/src/pred_at_pressure_levels_src.R')
theme_set(theme_bw() +theme(legend.position = 'bottom'))

base_folder_beginning <- Sys.getenv('base_folder_beginning')
base_folder_end <- Sys.getenv('base_folder_end')
G <- as.integer(Sys.getenv('G'))
source_folder <- Sys.getenv('source_folder')
m_pred_predictors <-  as.integer(Sys.getenv('m_pred_predictors'))
m_pred_response <- as.integer(Sys.getenv('m_pred_response'))
grid_size <- as.integer(Sys.getenv('grid_size'))
new_response_pcs2 <- as.integer(Sys.getenv('new_response_pcs2'))


save_folder <- paste0(source_folder,
                      base_folder_beginning,  G,
                      base_folder_end)
variable <- strsplit(base_folder_beginning, '_')[[1]][1]

load(paste0(save_folder, 'final_results.RData'))
estimated_model[['parameters']][['Omega_response2']] <- 
  lapply(estimated_model[['parameters']][['Omega_response2']], function(x) {
    x[,1:new_response_pcs2]
  })
estimated_model[['parameters']][['variances_response']] <- 
  lapply(estimated_model[['parameters']][['variances_response']], function(x) {
    x[1:new_response_pcs2]
  })
estimated_model[['parameters']][['range_params_r']] <- 
  lapply(estimated_model[['parameters']][['range_params_r']], function(x) {
    x[,1:new_response_pcs2]
  })
estimated_model[['params']][['pc_response2']] <- new_response_pcs2

params <- estimated_model[['params']]

# reorder clusters
parameters <- estimated_model[['parameters']]
Phi_pred <- bdiag(eval.basis(params[['basis_predictors']][[1]], evalarg = 1000),
                  eval.basis(params[['basis_predictors']][[1]], evalarg = 1000))
t_vals <- sapply(1:params[['G']], function(g){ (Phi_pred %*% parameters[['means_predictors']][[g]])[1,]})
new_order <- order(order(t_vals))
memberships_use <- new_order[parameters[['cluster_membership']]]

pressures <- c(2.5,   10.0,   20.0,   30.0,   40.0,   50.0,   60.0,   70.0,   80.0,   90.0,  100.0,  110.0,  120.0,
               130.0,  140.0,  150.0,  160.0,  170.0,  182.5,  200.0,  220.0,  240.0,  260.0,  280.0,  300.0,  320.0,
               340.0,  360.0,  380.0,  400.0,  420.0,  440.0,  462.5,  500.0,  550.0,  600.0,  650.0,  700.0,  750.0,
               800.0,  850.0,  900.0,  950.0, 1000.0, 1050.0, 1100.0, 1150.0, 1200.0, 1250.0, 1300.0, 1350.0, 1412.5,
               1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 1975.0)
grid_df <- get_grid_info(grid_size) %>%
  mutate(longitude = ifelse(longitude > 180, longitude - 360, longitude))
n_p <- length(pressures)
n_x <- length(unique(grid_df$longitude))
n_y <- length(unique(grid_df$latitude))

cluster_probabilities <- array(dim = c(n_x, n_y, G, 12), NA)

prediction_cluster_files <- list.files(save_folder, pattern = '^predict_clusters')
prediction_cluster_files_months <- sapply(strsplit(prediction_cluster_files, '-'),
                                          function(x) as.numeric(x[2]))


for (z in 1:12) {
  print(prediction_cluster_files[which(prediction_cluster_files_months == z)[1]])
  load(paste0(save_folder, prediction_cluster_files[which(prediction_cluster_files_months == z)[1]]))
  cluster_assigned <- apply(cluster_mat_pred, 1, function(x) which.max(tabulate(x)))
  cluster_assigned <- new_order[cluster_assigned]
  cluster_probability <- apply(cluster_mat_pred, 
                               1, function(x) tabulate(x)[which.max(tabulate(x))])/
    ncol(cluster_mat_pred)
  
  cluster_probs <- t(apply(cluster_mat_pred, 1, function(x) tabulate(x, nbins = G),
                         simplify = T)/ncol(cluster_mat_pred))
  
  probs_df <- cbind(grid_use, prob =cluster_probs) %>%
    right_join(grid_df) %>%
    arrange(latitude, longitude)
  cond_probs <- as.matrix(probs_df[, paste0('prob.', 1:G)])
  cluster_probabilities[,,,z] <- as.matrix(cond_probs)
}

dimx <- ncdim_def("Longitude", "Degrees", unique(grid_df$longitude)[order(unique(grid_df$longitude))]) 
dimy <- ncdim_def("Latitude", "Degrees",unique(grid_df$latitude)[order(unique(grid_df$latitude))]) 
dimg <- ncdim_def("Cluster", "Probability", 1:G) 
dimm <- ncdim_def("Month", "month", 1:12)
cluster_ncdf <- ncvar_def('CLUSTER_PROB' , units = 'PROBABILITY', dim = list(dimx, dimy, dimg, dimm), 
                          prec = 'double',
                          longname = 'CLUSTER_PROBABILITIES', verbose = T)

cluster_file <- nc_create(paste0(save_folder, 'cluster_product.nc'), 
                          vars = list(cluster_ncdf),
                          verbose = T, force_v4 = T)
ncvar_put(cluster_file, 'CLUSTER_PROB', vals = cluster_probabilities)
nc_close(cluster_file)

prediction_files <- list.files(save_folder, pattern = '^predict_IS')
n_dates <- length(prediction_files)
prediction_files_months <- sapply(strsplit(prediction_cluster_files, '-'),
                                          function(x) as.numeric(x[2]))
print(n_dates)

oxygen_predictions <- array(dim = c(n_x, n_y, n_p, n_dates), NA)
temperature_predictions <- array(dim = c(n_x, n_y, n_p, n_dates), NA)
salinity_predictions <- array(dim = c(n_x, n_y, n_p, n_dates), NA)
oxygen_var <- array(dim = c(n_x, n_y, n_p, n_dates), NA)
oxygen_var_p1 <- array(dim = c(n_x, n_y, n_p, n_dates), NA)

Phi <- eval.basis(params[['basis_response']], pressures)
Phi_pred <- bdiag(eval.basis(params[['basis_predictors']][[1]], evalarg= pressures),
                  eval.basis(params[['basis_predictors']][[2]], evalarg= pressures))

Omega_responseg <- lapply(1:params[['G']], function(g) {
  cbind(Phi %*% parameters[['Omega_response1']][[g]],
        Phi %*% parameters[['Omega_response2']][[g]])
})


# check if grid points are close to any profile
for (z in 1:n_dates) {
  print(prediction_files[z])
  load(paste0(save_folder, prediction_files[z]))
  print(pred_score_results[['wt']])
  pcs_array_predictors <- pred_score_results[['pcs_array_predictors']]
  pcs_array_response <- pred_score_results[['pcs_array_response']]
  MC_max_pred <- dim(pred_score_results[['pcs_array_response']])[2]
  cluster_mat_pred_use <- cluster_mat[pred_set_up[['is_prediction']],]
  # predictions
  values <- array(dim = c(nrow(pcs_array_predictors), MC_max_pred, length(pressures)))
  temp_values <- array(dim = c(nrow(pcs_array_predictors), MC_max_pred, length(pressures)))
  psal_values <- array(dim = c(nrow(pcs_array_predictors), MC_max_pred, length(pressures)))
  for (i in 1:nrow(pcs_array_predictors)) {
    for (mc in 1:MC_max_pred) { 
      z_i_tau <- cluster_mat_pred_use[i,mc]
      values[i,mc, ] <- 
        as.double(Phi %*% (parameters[['means_response']][[z_i_tau]] + 
                             parameters[['Omega_response1']][[z_i_tau]] %*% pcs_array_predictors[i,mc,] + 
                             parameters[['Omega_response2']][[z_i_tau]] %*% pcs_array_response[i,mc,]))
      
      pred_values <- as.double(Phi_pred %*% (parameters[['means_predictors']][[z_i_tau]] + 
                                               parameters[['Omega_predictors']][[z_i_tau]] %*% pcs_array_predictors[i,mc,]))
      
      temp_values[i,mc, ] <- pred_values[1:length(pressures)]
      psal_values[i,mc, ] <- pred_values[-c(1:length(pressures))]
    }
  }
  avg_values <- apply(values, c(1,3), function(x) mean(x * pred_score_results[['wt']]))
  avg_temp_values <- apply(temp_values, c(1,3), function(x) mean(x * pred_score_results[['wt']]))
  avg_psal_values <- apply(psal_values, c(1,3), function(x) mean(x * pred_score_results[['wt']]))
  
  var_values1 <- array(dim = c(nrow(pcs_array_predictors), MC_max_pred, length(pressures)))
  for (i in 1:nrow(pcs_array_predictors)) {
    for (mc in 1:MC_max_pred) { 
      z_i_tau <- cluster_mat_pred_use[i,mc]
      var_values1[i, mc, ] <- 
        rowSums((Omega_responseg[[z_i_tau]] %*% pred_score_results[['variances']][[mc]][[i]])*
                  Omega_responseg[[z_i_tau]])
    }
  }
  avg_var_values <- apply(var_values1, c(1,3), function(x) mean(x* pred_score_results[['wt']]))
  var_avg_values <- apply(values, c(1,3), function(x) {
    mean((x * pred_score_results[['wt']] - mean(x * pred_score_results[['wt']]))^2)
  })

  
  vals_df <- cbind(grid_use, oxy = avg_values,
                   temp = avg_temp_values,
                   psal = avg_psal_values,
                   var_p1 = avg_var_values,
                   var_full = var_avg_values + avg_var_values) %>%
    right_join(grid_df, by = c("latitude", "longitude", "depth", "width", "height")) %>%
    arrange(latitude, longitude)
  oxygen_predictions[,,,z] <- as.matrix(vals_df %>% dplyr::select(starts_with('oxy')))
  temperature_predictions[,,,z] <- as.matrix(vals_df %>% dplyr::select(starts_with('temp')))
  salinity_predictions[,,,z] <- as.matrix(vals_df %>% dplyr::select(starts_with('psal')))
  
  oxygen_var[,,,z] <- as.matrix(vals_df %>% dplyr::select(starts_with('var_full')))
  oxygen_var_p1[,,,z] <- as.matrix(vals_df %>% dplyr::select(starts_with('var_p1')))
}
dates <- sapply(strsplit(prediction_files, '_'), function(x) x[3])
julian_date <- julian(as.Date(dates), origin = as.Date('2000-01-01'))
dimp <- ncdim_def("Pressure", "Decibars", pressures) 
dimt <- ncdim_def("Time", "days", julian_date) 
oxygen_ncdf <- ncvar_def('OXYGEN', units = 'μmol/kg', dim = list(dimx, dimy, dimp, dimt), 
                         prec = 'double',
                         longname = 'OXYGEN PREDICTION', verbose = T)
oxygen_sd_ncdf <- ncvar_def('OXYGEN_SD' , units = 'μmol/kg', dim = list(dimx, dimy, dimp, dimt), 
                         prec = 'double',
                         longname = 'OXYGEN STANDARD DEVIATION', verbose = T)
oxygen_sd1_ncdf <- ncvar_def('OXYGEN_SD1' , units = 'μmol/kg', dim = list(dimx, dimy, dimp, dimt), 
                            prec = 'double',
                            longname = 'OXYGEN STANDARD DEVIATION CONDITIONAL ON CLUSTER', verbose = T)
temperature_ncdf <- ncvar_def('TEMPERATURE' , units = '°C', dim = list(dimx, dimy, dimp, dimt), 
                            prec = 'double',
                            longname = 'TEMPERATURE PREDICTION', verbose = T)
salinity_ncdf <- ncvar_def('SALINITY' , units = 'PRACTICAL SALINITY UNITS', dim = list(dimx, dimy, dimp, dimt), 
                              prec = 'double',
                              longname = 'SALINITY PREDICTION', verbose = T)

oxygen_mean_ncdf <- ncvar_def('OXYGEN_MEAN', units = 'μmol/kg', dim = list(dimx, dimy, dimp), 
                         prec = 'double',
                         longname = 'OXYGEN MEAN PREDICTION', verbose = T)

oxygen_cycle_ncdf <- ncvar_def('OXYGEN_MEAN_CYCLE', units = 'μmol/kg', dim = list(dimx, dimy, dimp, dimm), 
                              prec = 'double',
                              longname = 'OXYGEN MEAN ANNUAL CYCLE', verbose = T)

temperature_mean_ncdf <- ncvar_def('TEMPERATURE_MEAN', units = '°C', dim = list(dimx, dimy, dimp), 
                              prec = 'double',
                              longname = 'TEMPERATURE MEAN PREDICTION', verbose = T)

temperature_cycle_ncdf <- ncvar_def('TEMPERATURE_MEAN_CYCLE', units = '°C', dim = list(dimx, dimy, dimp, dimm), 
                               prec = 'double',
                               longname = 'TEMPERATURE MEAN ANNUAL CYCLE', verbose = T)

salinity_mean_ncdf <- ncvar_def('SALINITY_MEAN', units = 'PSU', dim = list(dimx, dimy, dimp), 
                                   prec = 'double',
                                   longname = 'SALINITY MEAN PREDICTION', verbose = T)

salinity_cycle_ncdf <- ncvar_def('SALINITY_MEAN_CYCLE', units = 'PSU', dim = list(dimx, dimy, dimp, dimm), 
                                    prec = 'double',
                                    longname = 'SALINITY MEAN ANNUAL CYCLE', verbose = T)

oxy_file <- nc_create(paste0(save_folder, 'oxy_product.nc'), 
          vars = list(oxygen_ncdf),
          verbose = T, force_v4 = T)
temp_file <- nc_create(paste0(save_folder, 'temp_product.nc'), 
                     vars = list(temperature_ncdf),
                     verbose = T, force_v4 = T)
salinity_file <- nc_create(paste0(save_folder, 'salinity_product.nc'), 
                           vars = list(salinity_ncdf),
                           verbose = T, force_v4 = T)
cluster_file <- nc_create(paste0(save_folder, 'cluster_product.nc'), 
                      vars = list(cluster_ncdf),
                      verbose = T, force_v4 = T)
oxy_sd_file <- nc_create(paste0(save_folder, 'oxy_sd_product.nc'), 
                          vars = list(oxygen_sd_ncdf),
                          verbose = T, force_v4 = T)
oxy_sd1_file <- nc_create(paste0(save_folder, 'oxy_sd1_product.nc'), 
                         vars = list(oxygen_sd1_ncdf),
                         verbose = T, force_v4 = T)

oxy_mean_file <- nc_create(paste0(save_folder, 'oxy_mean_product.nc'), 
                      vars = list(oxygen_mean_ncdf),
                      verbose = T, force_v4 = T)
oxy_cycle_file <- nc_create(paste0(save_folder, 'oxy_cycle_product.nc'), 
                      vars = list(oxygen_cycle_ncdf),
                      verbose = T, force_v4 = T)

temp_mean_file <- nc_create(paste0(save_folder, 'temp_mean_product.nc'), 
                           vars = list(temperature_mean_ncdf),
                           verbose = T, force_v4 = T)
temp_cycle_file <- nc_create(paste0(save_folder, 'temp_cycle_product.nc'), 
                            vars = list(temperature_cycle_ncdf),
                            verbose = T, force_v4 = T)

salinity_mean_file <- nc_create(paste0(save_folder, 'salinity_mean_product.nc'), 
                            vars = list(salinity_mean_ncdf),
                            verbose = T, force_v4 = T)
salinity_cycle_file <- nc_create(paste0(save_folder, 'salinity_cycle_product.nc'), 
                             vars = list(salinity_cycle_ncdf),
                             verbose = T, force_v4 = T)

ncvar_put(oxy_file, 'OXYGEN', vals = oxygen_predictions)
ncvar_put(oxy_sd_file, 'OXYGEN_SD', vals = oxygen_var)
ncvar_put(oxy_sd1_file, 'OXYGEN_SD1', vals = oxygen_var_p1)
ncvar_put(temp_file, 'TEMPERATURE', vals = temperature_predictions)
ncvar_put(salinity_file, 'SALINITY', vals = salinity_predictions)
ncvar_put(cluster_file, 'CLUSTER_PROB', vals = cluster_probabilities)

oxygen_mean_prediction <- apply(oxygen_predictions, c(1,2,3), mean)
oxygen_cycle_prediction <- sapply(1:12, function(x) {
  apply(oxygen_predictions[, ,, prediction_files_months == x,
                           drop = F], c(1,2,3), mean)
}, simplify = "array")
ncvar_put(oxy_mean_file, 'OXYGEN_MEAN', vals = oxygen_mean_prediction)
ncvar_put(oxy_cycle_file, 'OXYGEN_MEAN_CYCLE', vals = oxygen_cycle_prediction)

temperature_mean_prediction <- apply(temperature_predictions, c(1,2,3), mean)
temperature_cycle_prediction <- sapply(1:12, function(x) {
  apply(temperature_predictions[, ,, prediction_files_months == x,
                           drop = F], c(1,2,3), mean)
}, simplify = "array")
ncvar_put(temp_mean_file, 'TEMPERATURE_MEAN', vals = temperature_mean_prediction)
ncvar_put(temp_cycle_file, 'TEMPERATURE_MEAN_CYCLE', vals = temperature_cycle_prediction)

salinity_mean_prediction <- apply(salinity_predictions, c(1,2,3), mean)
salinity_cycle_prediction <- sapply(1:12, function(x) {
  apply(salinity_predictions[, ,, prediction_files_months == x,
                                drop = F], c(1,2,3), mean)
}, simplify = "array")
ncvar_put(salinity_mean_file, 'SALINITY_MEAN', vals = salinity_mean_prediction)
ncvar_put(salinity_cycle_file, 'SALINITY_MEAN_CYCLE', vals = salinity_cycle_prediction)

nc_close(oxy_file)
nc_close(oxy_mean_file)
nc_close(oxy_cycle_file)
nc_close(oxy_sd_file)
nc_close(oxy_sd1_file)
nc_close(temp_file)
nc_close(salinity_file)
nc_close(cluster_file)
nc_close(salinity_mean_file)
nc_close(salinity_cycle_file)
nc_close(temp_cycle_file)
nc_close(temp_mean_file)
q()