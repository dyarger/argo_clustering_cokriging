library(tidyverse)
source('paper/code/src/plot_src.R')
library(Matrix)
library(fda)
theme_set(theme_bw() +theme(legend.position = 'bottom'))

base_folder_beginning <- Sys.getenv('base_folder_beginning')
base_folder_end <- Sys.getenv('base_folder_end')
G <- as.integer(Sys.getenv('G'))
source_folder <- Sys.getenv('source_folder')
m_pred_predictors <-  as.integer(Sys.getenv('m_pred_predictors'))
m_pred_response <- as.integer(Sys.getenv('m_pred_response'))
grid_size <- as.integer(Sys.getenv('grid_size'))
date_to_predict <- Sys.getenv('date')
new_response_pcs2 <- as.integer(Sys.getenv('new_response_pcs2'))
MC_max_clusters <- as.integer(Sys.getenv('MC_max_clusters'))
MC_max_pred <- as.integer(Sys.getenv('MC_max_pred'))

save_folder <- paste0(source_folder,
                      base_folder_beginning,  G,
                      base_folder_end)
variable <- strsplit(base_folder_beginning, '_')[[1]][1]

if (variable == 'oxy') {
  types <- c('oxygen', 'temperature', 'salinity')
  type_label <- c('Oxygen\n(μmol/kg)', 'Temperature\n(°C)', 'Salinity\n(PSU)')
  type_label_var <- c('Oxygen SD\n(μmol/kg)', 'Temp SD\n(°C)', 'Salinity SD\n(PSU)')
  type_label_no_units <- c('Oxygen', 'Temperature', 'Salinity')
  type_label_name <- c('Oxygen', 'Predictors')
} else {
  types <- c('nitrate', 'temperature', 'salinity')
  type_label <- c('Nitrate\n(μmol/kg)', 'Temperature\n(°C)', 'Salinity\n(PSU)')
  type_label_var <- c('Nitrate SD\n(μmol/kg)', 'Temp SD\n(°C)', 'Salinity SD\n(PSU)')
  type_label_no_units <- c('Nitrate', 'Temperature', 'Salinity')
  type_label_name <- c('Nitrate', 'Predictors')
}
month <- as.integer(substr(date_to_predict, 6, 7))
year <- as.integer(substr(date_to_predict, 1, 4))
day <- as.integer(substr(date_to_predict, 9, 10))
date_label <- paste0(month.name[month], ' ', day, 'th, ', year)
date_label_no_year <- paste0(month.name[month], ' ', day, 'th')

tile_score <- geom_tile(aes(x = longitude, y = latitude, color = score,fill = score, 
                            height = height, width = width))

load(paste0(
  save_folder, 'final_results.RData'))

# reorder clusters
parameters <- estimated_model[['parameters']]
Phi_pred <- bdiag(eval.basis(params[['basis_predictors']][[1]], evalarg= 1000),
                  eval.basis(params[['basis_predictors']][[2]], evalarg=1000))
t_vals <- sapply(1:params[['G']], function(g){ (Phi_pred %*% parameters[['means_predictors']][[g]])[1,]})
new_order <- order(order(t_vals))
memberships_use <- new_order[parameters[['cluster_membership']]]

load(paste0(
  save_folder, 'predict_clusters_', 
  date_to_predict, '_', grid_size, '_', MC_max_clusters, '.RData'))

cluster_assigned <- apply(cluster_mat_pred, 1, function(x) which.max(tabulate(x)))
cluster_assigned <- new_order[cluster_assigned]
cluster_probability <- apply(cluster_mat_pred, 
                          1, function(x) tabulate(x)[which.max(tabulate(x))])/
  ncol(cluster_mat_pred)
plt <- ggplot(data =  data.frame(grid_use, membership_use = cluster_assigned,
                                 max_prob = cluster_probability))+
  geom_tile(aes(x = longitude, y = latitude, 
                fill = factor(membership_use), alpha = max_prob,
                height = height + .022, width = width +.022)) +
  SO_coord + SO_theme +
  fronts_dark + continents +
  latitude_lines + longitude_lines +
  labs(color ='Cluster', fill = 'Cluster', alpha = 'Probability',
       title = 'Predicted cluster membership',
       subtitle = date_label_no_year)+
  theme(panel.grid = element_blank()) +
  guides(alpha = guide_legend(ncol = 1), fill = guide_legend(ncol = 2))
ggsave(plot = plt, filename = paste0(save_folder, date_to_predict, '_pred_prob.png'), 
       width = 4.5, 
       height = 5.11, dpi = 550)

plt <- ggplot(data =  data.frame(grid_use, membership_use = cluster_probability))+
  geom_tile(aes(x = longitude, y = latitude, fill =membership_use,
                height = height + .022, width = width +.022)) +
  SO_coord + SO_theme +
  fronts_dark + continents + 
  latitude_lines + longitude_lines + 
  labs(color ='Probability', fill = 'Probability', 
       title = 'Probability of the most likely cluster')+
  theme(panel.grid = element_blank()) +
  scale_color_viridis_c()+  scale_fill_viridis_c()
ggsave(filename = paste0(save_folder, date_to_predict, '_clus_uncertainty.png'), plt, width = 4.5, 
       height = 5.11, dpi = 550)


# Predictions
#   load in prediction/parameters
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
parameters <- estimated_model[['parameters']]
params <- estimated_model[['params']]
load(paste0(
  save_folder, 'predict_IS_', 
  date_to_predict, '_', grid_size, '_', m_pred_predictors, '_',
  m_pred_response, '_', MC_max_pred, '.RData'))

pcs_array_predictors <- pred_score_results[['pcs_array_predictors']]
pcs_array_response <- pred_score_results[['pcs_array_response']]
pressures <- c(50, 150, 300, 500, 1000, 1500)
Phi <- eval.basis(params[['basis_response']], pressures)
Phi_pred <- bdiag(eval.basis(params[['basis_predictors']][[1]], evalarg= pressures),
                  eval.basis(params[['basis_predictors']][[2]], evalarg= pressures))
cluster_mat_pred_use <- cluster_mat[!pred_set_up[['has_data']],]


Omega_responseg <- lapply(1:params[['G']], function(g) {
  cbind(Phi %*% parameters[['Omega_response1']][[g]],
        Phi %*% parameters[['Omega_response2']][[g]])
})

Omega_predictorsg <- lapply(1:params[['G']], function(g) {
  Phi_pred %*% parameters[['Omega_predictors']][[g]]
})

var_values <- pred_score_results[['variances']]
n_pcs_predictors <- pred_set_up[['prediction_model']][['params']][['pc_predictors']]

for (type in 1:length(types)) {
  
  # predictions
  values <- array(dim = c(nrow(pcs_array_predictors), MC_max_pred, length(pressures)))
  for (i in 1:nrow(pcs_array_predictors)) {
    for (mc in 1:MC_max_pred) { 
      z_i_tau <- cluster_mat_pred_use[i,mc]
      if (types[type] == 'oxygen' | types[type] == 'nitrate') {
        values[i,mc, ] <- 
          as.double(Phi %*% (parameters[['means_response']][[z_i_tau]] + 
                               parameters[['Omega_response1']][[z_i_tau]] %*% pcs_array_predictors[i,mc,] + 
                               parameters[['Omega_response2']][[z_i_tau]] %*% pcs_array_response[i,mc,]))
      } else if (types[type] == 'temperature') {
        values[i,mc, ] <- 
          as.double(Phi_pred %*% (parameters[['means_predictors']][[z_i_tau]] + 
                               parameters[['Omega_predictors']][[z_i_tau]] %*% pcs_array_predictors[i,mc,]))[1:length(pressures)]
      } else if (types[type] == 'salinity') {
        values[i,mc, ] <- 
          as.double(Phi_pred %*% (parameters[['means_predictors']][[z_i_tau]] + 
                               parameters[['Omega_predictors']][[z_i_tau]] %*% pcs_array_predictors[i,mc,]))[-c(1:length(pressures))]
      }
    }
  }
  avg_values <- apply(values,c(1,3), function(x) mean(x * pred_score_results[['wt']]))
  
  # first part of variance
  var_values1 <- array(dim = c(nrow(pcs_array_predictors), MC_max_pred, length(pressures)))
  for (i in 1:nrow(pcs_array_predictors)) {
    for (mc in 1:MC_max_pred) { 
      z_i_tau <- cluster_mat_pred_use[i,mc]
      if (types[type] == 'oxygen'| types[type] == 'nitrate') {
        var_values1[i, mc, ] <- 
          rowSums((Omega_responseg[[z_i_tau]] %*% var_values[[mc]][[i]])*Omega_responseg[[z_i_tau]])
      } else if (types[type] == 'temperature') {
        var_values1[i, mc, ] <- 
          rowSums((Omega_predictorsg[[z_i_tau]][1:length(pressures), ] %*% 
                     var_values[[mc]][[i]][1:n_pcs_predictors, 1:n_pcs_predictors])*
                    Omega_predictorsg[[z_i_tau]][1:length(pressures), ]) 
      } else if (types[type] == 'salinity') {
        var_values1[i, mc, ] <- 
          rowSums((Omega_predictorsg[[z_i_tau]][-c(1:length(pressures)), ] %*% 
                     var_values[[mc]][[i]][1:n_pcs_predictors, 1:n_pcs_predictors])*
                    Omega_predictorsg[[z_i_tau]][-c(1:length(pressures)), ])
      }
    }
    if (i %% 1000 == 0) {
      print(i)
    }
  }
  avg_var_values <- apply(var_values1, c(1,3), function(x) mean(x* pred_score_results[['wt']]))
  # second part of variance
  var_avg_values <- apply(values, c(1,3), function(x) mean((x * pred_score_results[['wt']] -
                                                             mean(x * pred_score_results[['wt']]))^2))
  
  for (p in 1:length(pressures)) {
    a <- ggplot(data = data.frame(grid_use, score = avg_values[,p]))+
      tile_score + scale_color_viridis_c() + scale_fill_viridis_c() + 
      fronts_light + SO_coord + SO_theme + continents + latitude_lines + longitude_lines + 
      labs(fill = type_label[type], color = type_label[type],
           title = paste0(type_label_no_units[type], ' prediction'),
           subtitle = paste0(pressures[p], ' decibars, ', date_label))+ 
      theme(text = element_text(size = 16),
            legend.key.width = unit(1.2, "cm"), legend.text = element_text(size = 14)) 
    ggsave(filename = paste0(save_folder, '/', date_to_predict, '_', pressures[p], '_', types[type], '.png'), a, 
           width = 4.5, height = 5.11)
    
    a <- ggplot(data = data.frame(grid_use, score = sqrt(avg_var_values[,p])))+
      tile_score + 
      scale_color_viridis_c() + scale_fill_viridis_c() +   SO_coord + SO_theme +
      fronts_light + continents + latitude_lines + longitude_lines + 
      labs(fill = type_label_var[type], color = type_label_var[type],
           subtitle = paste0(pressures[p], ' decibars, ', date_label),
           title = paste0(type_label_no_units[type], ' avg conditional SD')) + 
      theme(text = element_text(size = 16),  legend.key.width = unit(1, "cm"))
    ggsave(filename =  paste0(save_folder, date_to_predict, '_var_p1_', pressures[p], '_', types[type], '.png'), a,
           width = 4.5, height = 5.11)
    a <- ggplot(data = data.frame(grid_use, score = sqrt(var_avg_values[,p])))+
      tile_score + 
      scale_color_viridis_c() + scale_fill_viridis_c() +   SO_coord + SO_theme +
      fronts_light + continents + latitude_lines + longitude_lines + 
      labs(fill = type_label_var[type], color = type_label_var[type],
           subtitle = paste0(pressures[p], ' decibars, ', date_label),
           title = paste0(type_label_no_units[type], ' uncertainty from clusters')) + 
      theme(text = element_text(size = 16),  legend.key.width = unit(1, "cm"))
    ggsave(filename =  paste0(save_folder, date_to_predict, '_var_p2_',  pressures[p], '_', types[type], '.png'), a,
           width = 4.5, height = 5.11)
    
    a <- ggplot(data = data.frame(grid_use, score = sqrt(var_avg_values[,p]+avg_var_values[,p])))+
      tile_score + scale_color_viridis_c() + scale_fill_viridis_c() +   SO_coord + SO_theme +
      fronts_light + continents + latitude_lines + longitude_lines + 
      labs(fill =type_label_var[type], color =type_label_var[type],
           subtitle = paste0(pressures[p], ' decibars, ', date_label),
           title = paste0(type_label_no_units[type], ' standard deviation')) + 
      theme(text = element_text(size = 16),  legend.key.width = unit(1, "cm"))
    ggsave(filename =  paste0(save_folder, date_to_predict, '_var_', pressures[p], '_', types[type], '.png'), a,
           width = 4.5, height = 5.11)
  }
  
  if (types[type] %in% c('temperature', 'salinity')) {
    next
  }
  
  # compute unconditional predictions/variances
  # predictions
  uncon_values <- array(dim = c(nrow(pcs_array_predictors), MC_max_pred, length(pressures)))
  for (i in 1:nrow(pcs_array_predictors)) {
    for (mc in 1:MC_max_pred) { 
      z_i_tau <- cluster_mat_pred_use[i,mc]
      uncon_values[i,mc, ] <- 
        as.double(Phi %*% (parameters[['means_response']][[z_i_tau]]))
    }
  }
  uncon_avg_values <- apply(uncon_values,c(1,3), function(x) mean(x* pred_score_results[['wt']]))
  
  # first part of variance
  uncon_vars_response <- lapply(1:length(Omega_responseg), function(g){
    diag((Omega_responseg[[g]] %*% diag(c(parameters[['variances_predictors']][[g]], 
                                          parameters[['variances_response']][[g]]))) %*%
            t(Omega_responseg[[g]]))
  })
  uncon_var_values1 <- array(dim = c(nrow(pcs_array_predictors), MC_max_pred, length(pressures)))
  for (i in 1:nrow(pcs_array_predictors)) {
    for (mc in 1:MC_max_pred) { 
      z_i_tau <- cluster_mat_pred_use[i,mc]
      uncon_var_values1[i, mc, ] <- 
        uncon_vars_response[[z_i_tau]]
    }
  }
  uncon_avg_var_values <- apply(uncon_var_values1, c(1,3), function(x) mean(x* pred_score_results[['wt']]))
  # second part of variance
  uncon_var_avg_values <- apply(uncon_values, c(1,3), function(x) mean((x * pred_score_results[['wt']] -
                                                               mean(x * pred_score_results[['wt']]))^2))
  
  for (p in 1:length(pressures)) {
    a <- ggplot(data = data.frame(grid_use, score = 1- (avg_var_values[,p] + var_avg_values[,p])/
                                    (uncon_avg_var_values[,p] + uncon_var_avg_values[,p])) %>%
                  mutate(score = ifelse(score < 0, 0, score), 
                         score = ifelse(score > 1, 1, score)))+
      tile_score + scale_color_viridis_c() + scale_fill_viridis_c() +   SO_coord + SO_theme +
      fronts_light + continents + latitude_lines + longitude_lines + 
      labs(fill = substitute(paste(nn, ~R^2), list(nn = type_label_name[1])), 
           color = substitute(paste(nn, ~R^2), list(nn = type_label_name[1])),
           title = substitute(paste(nn, ~R^2), list(nn = type_label_name[1])),
           subtitle = paste0('Pressure ', pressures[p], ', ', date_label)) + 
      theme(text = element_text(size = 16),
            legend.key.width = unit(1, "cm"))
    ggsave(filename = paste0(save_folder, date_to_predict, '_r2_', pressures[p], '_', types[type], '.png'), a,
           width = 4.5, height = 5.11)
  }
}


############# Section plots

pressures <- c(seq(0, 2000, by = 5))
final_prediction <- list()
longitude_use <- -90.0833358764648
tol <- .5
grid_use_subset <- grid_use %>%
  filter(abs(longitude -longitude_use) < tol)
var_values_subset <- lapply(var_values, function(x) {x[abs(grid_use$longitude   -longitude_use) < tol]})

Phi <- eval.basis(params[['basis_response']], pressures)
Omega_responseg_subset <- lapply(1:params[['G']], function(g) {
  cbind(Phi %*% parameters[['Omega_response1']][[g]],
        Phi %*% parameters[['Omega_response2']][[g]])
})
pcs_array_predictors_subset <- pcs_array_predictors[abs(grid_use$longitude 
                                                       -longitude_use) < tol,,]
pcs_array_response_subset <- pcs_array_response[abs(grid_use$longitude 
                                                        -longitude_use) < tol,,]
values_reduced <- array(dim = c(nrow(pcs_array_predictors_subset), MC_max_pred, length(pressures)))
cluster_mat_subset <- cluster_mat_pred_use[abs(grid_use$longitude 
                                               -longitude_use) < tol,]
var_values1_reduced <- array(dim = c(nrow(pcs_array_predictors_subset), MC_max_pred, length(pressures)))
values_df <- list()
for (i in 1:nrow(pcs_array_predictors_subset)) {
  for (mc in 1:MC_max_pred) { 
    z_i_tau <- cluster_mat_subset[i,mc]
    values_reduced[i,mc, ] <- 
      as.double(Phi %*% (parameters[['means_response']][[z_i_tau]] + 
                           parameters[['Omega_response1']][[z_i_tau]] %*% pcs_array_predictors[i,mc,] + 
                           parameters[['Omega_response2']][[z_i_tau]] %*% pcs_array_response[i,mc,]))
    var_values1_reduced[i, mc, ] <- 
      rowSums((Omega_responseg_subset[[z_i_tau]] %*% var_values_subset[[mc]][[i]])*
                Omega_responseg_subset[[z_i_tau]])
    # first part of variance
  }
  
  values_df[[i]] <- data.frame('oxygen' = apply(values_reduced[i,,],2, function(x) mean(x* pred_score_results[['wt']])),
                               'v_p1' = apply(var_values1_reduced[i,,],2, function(x) mean(x* pred_score_results[['wt']])),
                               'v_p2' = apply(values_reduced[i,,], 2, function(x) mean((x * pred_score_results[['wt']] -
                                                                              mean(x * pred_score_results[['wt']]))^2)),
                               'oxygen_var' = apply(var_values1_reduced[i,,],2, function(x) mean(x* pred_score_results[['wt']])) +
                                 apply(values_reduced[i,,], 2, function(x) mean((x * pred_score_results[['wt']] -
                                                                                 mean(x * pred_score_results[['wt']]))^2)),
                               'pressure' = pressures,
                               'latitude' = grid_use_subset[i,'latitude'])
  
}
values_df <- bind_rows(values_df) %>% left_join(grid_use_subset)

ggplot(data = values_df, 
       aes(y = pressure, x = latitude, color = oxygen, width = height, fill = oxygen))+
  geom_tile() + 
  scale_color_viridis_c(name = type_label) +
  scale_fill_viridis_c(name = type_label)  +
  scale_y_continuous(trans = 'reverse') + 
  labs(x = 'Latitude', y = 'Pressure (decibars)',
       title = paste0(type_label_no_units[1], ', longitude ', round(longitude_use)),
       subtitle = date_label)+ 
  theme(text = element_text(size = 16), 
        legend.key.width = unit(1, "cm"))
ggsave(filename = paste0(save_folder, date_to_predict,  '_', round(longitude_use), '_section.png'),  
       width = 8.5, 
       height = 5.11)

ggplot(data = values_df, aes(y = pressure, x = latitude, color = sqrt(oxygen_var),
                                    width = height, fill = sqrt(oxygen_var)))+
  geom_tile() + 
  scale_color_viridis_c(name = type_label_var) +
  scale_fill_viridis_c(name = type_label_var)  +
  scale_y_continuous(trans = 'reverse') + 
  labs(x = 'Latitude', y = 'Pressure (decibars)',
       title = paste0(type_label_no_units[1], ' standard deviation, longitude ', round(longitude_use)),
       subtitle = date_label)+ 
  theme(text = element_text(size = 16), 
        legend.key.width = unit(1, "cm"))
ggsave(filename = paste0(save_folder, date_to_predict,  '_', round(longitude_use), '_section_var.png'),  
       width = 8.5, 
       height = 5.11)

latitudes_use <- unique(values_df[['latitude']])[seq(1, length(unique(values_df[['latitude']])), 3)]

ggplot(data = values_df %>% filter(latitude %in% latitudes_use), aes(color = latitude, y = pressure, x = oxygen,
                                           group = latitude))+
  geom_path() + 
  scale_y_continuous(trans = 'reverse') + 
  labs(x = type_label[1], y = 'Pressure (decibars)', color = 'Latitude',
       title = paste0(type_label_no_units[1], ', longitude ', round(longitude_use)),
       subtitle = date_label)+ 
  theme(text = element_text(size = 16), 
        legend.key.width = unit(1, "cm"))+
  scale_color_viridis_c()
ggsave(filename = paste0(save_folder, date_to_predict,  '_', round(longitude_use), '_section_oxy_lines.png'),  
       width = 4, 
       height = 7.11)

ggplot(data = values_df %>% filter(latitude %in% latitudes_use), aes(color = latitude, y = pressure, x = sqrt(oxygen_var),
                                           group = latitude))+
  geom_path() + 
  scale_y_continuous(trans = 'reverse') + 
  labs(x = type_label_var[1], y = 'Pressure (decibars)', color = 'Latitude',
       title = paste0(type_label_no_units[1], ', longitude ', round(longitude_use)),
       subtitle = date_label)+ 
  theme(text = element_text(size = 16), 
        legend.key.width = unit(1, "cm"))+
  scale_color_viridis_c()
ggsave(filename = paste0(save_folder, date_to_predict,  '_', round(longitude_use), '_section_oxy_var.png'),  
       width = 4, 
       height = 7.11)

