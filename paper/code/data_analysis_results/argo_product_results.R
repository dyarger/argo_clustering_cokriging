library(Matrix)
library(fda)
library(ncdf4)
library(tidyverse)
library(reshape2)
library(fields)
theme_set(theme_bw())
source('paper/code/src/pred_at_pressure_levels_src.R')
grid_info <- get_grid_info(1) %>%
  mutate(longitude = ifelse(longitude > 180, longitude - 360, longitude))
source('paper/code/src/plot_src.R')

cluster_info <- nc_open('paper/results/oxy_none_NA_5_13_12_15_15/cluster_product.nc')
cluster_probs <- ncvar_get(cluster_info, varid = 'CLUSTER_PROB')

latitude <- cluster_info[['dim']][['Latitude']][['vals']]
longitude <- cluster_info[['dim']][['Longitude']][['vals']]
cluster <- cluster_info[['dim']][['Cluster']][['vals']]
month <- cluster_info[['dim']][['month']][['vals']]

clust_probs <- array(cluster_probs, dimnames = list(longitude, latitude, cluster, month),
                     dim = dim(cluster_probs))
prob_df <- reshape2::melt(clust_probs)
colnames(prob_df) <- c('longitude', 'latitude', 'cluster', 'month', 'value')

prob_df_sum <- prob_df %>%
  group_by(longitude, latitude, month) %>%
  arrange(cluster) %>%
  summarize(cluster = which.max(value),
            probability =value[which.max(value)])
prob_df_sum <- merge(prob_df_sum, grid_info, by = c('longitude', 'latitude'))

time_df <- data.frame(month = 1:12, 
                      month_name = factor(month.name, levels = month.name))


load('paper/results/oxy_none_NA_5_13_12_15_15/final_results.RData')
parameters <- estimated_model[['parameters']]
Phi_pred <- bdiag(eval.basis(params[['basis_predictors']][[1]], evalarg= 1000),
                  eval.basis(params[['basis_predictors']][[1]], evalarg= 1000))
t_vals <- sapply(1:params[['G']], function(g){ (Phi_pred %*% parameters[['means_predictors']][[g]])[1,]})
new_order <- order(order(t_vals))
memberships_use <- new_order[parameters[['cluster_membership']]]

prob_df_sum <- mutate(prob_df_sum, cluster_use = new_order[cluster])

a <- ggplot()+
  geom_tile(data = prob_df_sum %>% 
              left_join(time_df),
            aes(x = ifelse(longitude > 180, longitude - 360, longitude), y = latitude,
                fill = probability,
                height = height, width = width))+
  SO_coord + SO_theme +
  scale_fill_viridis_c() + 
  fronts_dark + continents + latitude_lines + longitude_lines + 
  facet_wrap(~month_name, ncol = 4)+
  labs(color = 'Maximum Probability', fill = 'Maximum Probability') + 
  theme(legend.position = 'bottom')
ggsave(filename = paste0('paper/results/oxy_none_NA_5_13_12_15_15/', 'product_probabilities.png'), 
       plot = a, width = 8, 
       height = 7)

a <- ggplot()+
  geom_tile(data = prob_df_sum %>% 
              left_join(time_df),
            aes(x = longitude, y = latitude,
                fill = factor(cluster_use),
                height = height, width = width, alpha = probability))+
  SO_coord + SO_theme +
  fronts_dark + continents + latitude_lines + longitude_lines + 
  facet_wrap(~month_name, ncol = 4)+
  labs(alpha = 'Probability', color = 'Cluster', fill = 'Cluster')
ggsave(filename = paste0('paper/results/oxy_none_NA_5_13_12_15_15/', 'product_cluster.png'), 
       plot = a, width = 8, 
       height = 7)

  
  
oxy_mean_nc <- nc_open('paper/results/oxy_none_NA_5_13_12_15_15/oxy_mean_product.nc')
oxy_mean <- ncvar_get(oxy_mean_nc, varid = 'OXYGEN_MEAN')
pressures <- oxy_mean_nc[['dim']][['Pressure']][['vals']]
oxy_mean <- array(oxy_mean, dimnames = list(longitude, latitude, pressures),
                     dim = dim(oxy_mean))
oxy_mean_df <- reshape2::melt(oxy_mean)
colnames(oxy_mean_df) <- c('longitude', 'latitude', 'pressure', 'value')
oxy_mean_df <- merge(oxy_mean_df, grid_info, by = c('longitude', 'latitude'))

a <- ggplot()+
  geom_tile(data = oxy_mean_df %>% filter(pressure == 150) %>%
              filter(!is.na(value)),
            aes(x = longitude, y = latitude, color = value,
                fill = value,
                height = height, width = width)) + 
  scale_fill_viridis_c() + scale_color_viridis_c() + 
  SO_coord + SO_theme +
  fronts_dark + continents + latitude_lines + longitude_lines + 
  labs(color ='Oxygen', fill = 'Oxygen',
       title = 'Oxygen Averaged Mean')+
  theme(panel.grid = element_blank())
ggsave(filename = paste0('paper/results/oxy_none_NA_5_13_12_15_15/', 'product_mean_150.png'), 
       plot = a, width = 4.5, 
       height = 5.11)

a <- ggplot()+
  geom_tile(data = oxy_mean_df %>% filter(pressure %in% c(150, 400, 1000)) %>%
              filter(!is.na(value)),
            aes(x = longitude, y = latitude, color = value,
                fill = value,
                height = height, width = width)) + 
  scale_fill_viridis_c() + scale_color_viridis_c() + 
  SO_coord + SO_theme +
  fronts_dark + continents + latitude_lines + longitude_lines + 
  facet_wrap(~pressure) + 
  labs(color ='Oxygen', fill = 'Oxygen',
       title = 'Oxygen Averaged Mean')+
  theme(panel.grid = element_blank())
ggsave(filename = paste0('paper/results/oxy_none_NA_5_13_12_15_15/', 'product_mean_150_400_1000.png'), 
       plot = a, width = 4.5, 
       height = 5.11)


oxy_cycle_nc <- nc_open('paper/results/oxy_none_NA_5_13_12_15_15/oxy_cycle_product.nc')
oxy_cycle <- ncvar_get(oxy_cycle_nc, varid = 'OXYGEN_MEAN_CYCLE')
months <- oxy_cycle_nc[['dim']][['Month']][['vals']]
oxy_cycle <- array(oxy_cycle, dimnames = list(longitude, latitude, pressures,
                                              factor(month.name[months], levels = month.name)),
                  dim = dim(oxy_cycle))
oxy_cycle_df <- reshape2::melt(oxy_cycle)
colnames(oxy_cycle_df) <- c('longitude', 'latitude', 'pressure', 'month', 'value')
oxy_cycle_df[['month']] <- factor(oxy_cycle_df[['month']],
                                  levels = month.name)

temp_cycle_nc <- nc_open('paper/results/oxy_none_NA_5_13_12_15_15/temp_cycle_product.nc')
temp_cycle <- ncvar_get(temp_cycle_nc, varid = 'TEMPERATURE_MEAN_CYCLE')
temp_cycle <- array(temp_cycle, dimnames = list(longitude, latitude, pressures,
                                              factor(month.name[months], levels = month.name)),
                   dim = dim(temp_cycle))
temp_cycle_df <- reshape2::melt(temp_cycle)
colnames(temp_cycle_df) <- c('longitude', 'latitude', 'pressure', 'month', 'value')
temp_cycle_df[['month']] <- factor(temp_cycle_df[['month']],
                                  levels = month.name)

salinity_cycle_nc <- nc_open('paper/results/oxy_none_NA_5_13_12_15_15/salinity_cycle_product.nc')
salinity_cycle <- ncvar_get(salinity_cycle_nc, varid = 'SALINITY_MEAN_CYCLE')
salinity_cycle <- array(salinity_cycle, dimnames = list(longitude, latitude, pressures,
                                                factor(month.name[months], levels = month.name)),
                    dim = dim(salinity_cycle))
salinity_cycle_df <- reshape2::melt(salinity_cycle)
colnames(salinity_cycle_df) <- c('longitude', 'latitude', 'pressure', 'month', 'value')
salinity_cycle_df[['month']] <- factor(salinity_cycle_df[['month']],
                                   levels = month.name)


pressures <- c(50, 150, 300, 500, 1000, 1500)
for (p in pressures) {
  oxy_cycle_df_one <- merge(oxy_cycle_df %>%filter(pressure == p), grid_info,
                            by = c('longitude', 'latitude'))
  temp_cycle_df_one <- merge(temp_cycle_df %>%filter(pressure == p), grid_info,
                            by = c('longitude', 'latitude'))
  salinity_cycle_df_one <- merge(salinity_cycle_df %>%filter(pressure == p), grid_info,
                             by = c('longitude', 'latitude'))
  a <- ggplot()+
    geom_tile(data = oxy_cycle_df_one %>%
                filter(!is.na(value)),
              aes(x = longitude, y = latitude, color = value,
                  fill = value,
                  height = height, width = width)) + 
    scale_fill_viridis_c() + scale_color_viridis_c() + 
    SO_coord + SO_theme +
    fronts_dark + continents +
    latitude_lines + longitude_lines +
    facet_wrap(~month, nrow = 3, ncol = 4) + 
    labs(color ='Oxygen\n(μmol/kg)', fill = 'Oxygen\n(μmol/kg)',
         title = 'Oxygen Annual Cycle')+
    theme(panel.grid = element_blank())
  ggsave(filename = paste0('paper/results/oxy_none_NA_5_13_12_15_15/', 'product_cycle_',
                           p, '.png'), 
         plot = a, width = 8, 
         height = 7)
  
  a <- ggplot()+
    geom_tile(data = temp_cycle_df_one %>%
                filter(!is.na(value)),
              aes(x = longitude, y = latitude, color = value,
                  fill = value,
                  height = height, width = width)) + 
    scale_fill_viridis_c() + scale_color_viridis_c() + 
    SO_coord + SO_theme +
    fronts_dark + continents +
    latitude_lines + longitude_lines +
    facet_wrap(~month, nrow = 3, ncol = 4) + 
    labs(color ='Temperature\n(°C)', fill = 'Temperature\n(°C)',
         title = 'Temperature Annual Cycle')+
    theme(panel.grid = element_blank())
  ggsave(filename = paste0('paper/results/oxy_none_NA_5_13_12_15_15/', 'product_cycle_temp_',
                           p, '.png'), 
         plot = a, width = 8, 
         height = 7)
  
  a <- ggplot()+
    geom_tile(data = salinity_cycle_df_one %>%
                filter(!is.na(value)),
              aes(x = longitude, y = latitude, color = value,
                  fill = value,
                  height = height, width = width)) + 
    scale_fill_viridis_c() + scale_color_viridis_c() + 
    SO_coord + SO_theme +
    fronts_dark + continents +
    latitude_lines + longitude_lines +
    facet_wrap(~month, nrow = 3, ncol = 4) + 
    labs(color ='PSU', fill = 'PSU',
         title = 'Salinity Annual Cycle')+
    theme(panel.grid = element_blank())
  ggsave(filename = paste0('paper/results/oxy_none_NA_5_13_12_15_15/', 'product_cycle_salinity_',
                           p, '.png'), 
         plot = a, width = 8, 
         height = 7)
  
}
