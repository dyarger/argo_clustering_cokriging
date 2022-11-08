.libPaths("/Library/Frameworks/R.framework/Versions/4.1/Resources/miworkspace-library")
library(dplyr)
library(ggplot2)
library(patchwork)
library(fda)
library(Matrix)
theme_set(theme_bw() + theme(legend.position = 'bottom'))
source('paper/code/src/plot_src.R')
G <- as.integer(Sys.getenv('G'))
source_folder <- Sys.getenv('source_folder')
base_folder_beginning <- Sys.getenv('base_folder_beginning')
variable <- strsplit(base_folder_beginning, '_')[[1]][1]
base_folder_end <- Sys.getenv('base_folder_end')
load(paste0(source_folder, 
            base_folder_beginning, 
            G, 
            base_folder_end,
            'final_results.RData'))
parameters <- estimated_model[['parameters']]
data_inputs <- estimated_model[['data_inputs']]
params <- estimated_model[['params']]
# reorder clusters based on temperature at 1000 dbar
Phi_pred <- bdiag(eval.basis(params[['basis_predictors']][[1]],evalarg = 1000),
                  eval.basis(params[['basis_predictors']][[2]],evalarg = 1000))
t_vals <- sapply(1:params[['G']], function(g){ (Phi_pred %*% parameters[['means_predictors']][[g]])[1,]})
new_order <- order(order(t_vals))
memberships_use <- new_order[parameters[['cluster_membership']]]

p_seq <- seq(0, 2000, by = 2.5)
Phi_pred <- eval.basis(params[['basis_response']], evalarg = p_seq)

nbasis_pred <- nrow(parameters[['Omega_predictors']][[1]])/length(params[['X']])

df_means <- data.frame('pressure' = p_seq, 
                       'mean' = as.double(sapply((1:params[['G']]), function(g){ Phi_pred %*% parameters[['means_response']][[g]]})),
                       'g' = factor(rep(new_order[(1:params[['G']])], each = length(p_seq))))

if (params[['Y']] == 'oxy') {
  var_name <- 'Oxygen'
} else if(params[['Y']] == 'nitrate') {
  var_name <- 'Nitrate'
}
label_use <- paste(var_name, '(μmol/kg)')
ggplot(data = df_means, aes(x = pressure, y = mean, group = g, color = g, linetype = g)) + 
  geom_line(size = 1.1) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = label_use,
       linetype = 'Cluster', color = 'Cluster', 
       title = paste(var_name, 'mean functions'))+
  theme(legend.position = 'bottom', text = element_text(size = 18))
ggsave(paste0(source_folder, 
              base_folder_beginning, 
              G, 
              base_folder_end, 'spatial_var_means.png'),
       height = 8.86, width = 5.11)

df_means <- data.frame('pressure' = p_seq, 
                       'mean' = as.double(sapply(1:params[['G']], function(g){ Phi_pred %*% parameters[['means_predictors']][[g]][1:nbasis_pred]})),
                       'g' = factor(rep(new_order[1:params[['G']]], each = length(p_seq))))

ggplot(data = df_means, aes(x = pressure, y = mean, group = g, color = g, linetype = g)) + 
  geom_line(size = 1.1) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = 'Temperature (°C)',
       linetype = 'Cluster', color = 'Cluster',
       title = 'Temperature mean functions')+
  theme(legend.position = 'bottom', text = element_text(size = 18))
ggsave(paste0(source_folder, 
              base_folder_beginning, 
              G, 
              base_folder_end, 'spatial_temp_means.png'),
       height = 8.86, width = 5.11)

df_means <- data.frame('pressure' = p_seq, 
                       'mean' = as.double(sapply(1:params[['G']], function(g){ Phi_pred %*% parameters[['means_predictors']][[g]][-(1:nbasis_pred)]})),
                       'g' = factor(rep(new_order[1:params[['G']]], each = length(p_seq))))

ggplot(data = df_means, aes(x = pressure, y = mean, group = g, color = g, linetype = g)) + 
  geom_line(size = 1.1) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = 'Salinity (PSU)',
       linetype = 'Cluster', color = 'Cluster',
       title = 'Salinity mean functions')+
  theme(legend.position = 'bottom', text = element_text(size = 18))
ggsave(paste0(source_folder, 
              base_folder_beginning, 
              G, 
              base_folder_end, 'spatial_psal_means.png'),
       height = 8.86, width = 5.11)

var_TS <- lapply(1:params[['G']], function(g) {
  diag(Phi_pred %*% parameters[['Omega_response1']][[g]] %*% diag(parameters[['variances_predictors']][[g]]) %*%
         t(Phi_pred %*% parameters[['Omega_response1']][[g]]))
})

var_resp <- lapply(1:params[['G']], function(g) {
  diag(Phi_pred %*% parameters[['Omega_response2']][[g]] %*% diag(parameters[['variances_response']][[g]]) %*%
         t(Phi_pred %*% parameters[['Omega_response2']][[g]]))
})
df_var <- data.frame(cluster = rep(new_order[1:params[['G']]], each = length(var_TS[[1]])),
                     pressure = p_seq,
                     var_TS = sqrt(unlist(var_TS) + unlist(var_resp) + 
                                     parameters[['measurement_error_response']]), var_resp = sqrt(
                       unlist(var_resp)),
                     var_me = sqrt(
                         parameters[['measurement_error_response']]))

df_var_labels <- data.frame('Type' = factor(c('Total SD', 'SD if \nTS Known', 'Measurement\nError SD'),
                                            levels = c('Measurement\nError SD', 'Total SD', 'SD if \nTS Known')),
                            'type' = c('var_TS', 'var_resp', 'var_me'))
df_var_long <- tidyr::pivot_longer(df_var, cols = starts_with('var'),  names_to = 'type', values_to = 'SD')
ggplot(data = df_var_long %>% left_join(df_var_labels), aes(pressure, SD, color = Type)) + 
  geom_line() + 
  facet_wrap(~cluster, ncol = 1) + 
  theme(legend.position = 'bottom') + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse', name = 'Pressure (decibars)') + 
  labs(title = 'Decomposition of variation in oxygen')
ggsave(paste0(source_folder, 
              base_folder_beginning, 
              G, 
              base_folder_end, 'spatial_standard_deviation_cluster.png'),height = 8, width = 4)

r2_by_cluster <- df_var %>%
  mutate(var_proportion = 1 - var_resp / var_TS)  %>%
  group_by(cluster) %>%
  summarise(mean(var_proportion))
plot(r2_by_cluster)

ggplot(data = data.frame(time = data_inputs[['time']] %% 365.25, mem = memberships_use)) + 
  geom_histogram(aes(x = time), breaks = seq(0, 365.25, by = 14))+
  facet_wrap(~mem, ncol = 5) +
  labs(x = 'Day of Year', y ='Number of Profiles',
       title = 'Temporal distribution of cluster assignment of profiles') 
ggsave(paste0(source_folder, 
              base_folder_beginning, 
              G, 
              base_folder_end, 'spatial_cluster_time.png'), height = 2, width = 6)

label_use <- paste(var_name, 'function (μmol/kg)')
df_pcs <- data.frame('pressure' = p_seq, 
                       'pc' = as.double(sapply((1:params[['G']]), function(g){ Phi_pred %*% parameters[['Omega_response1']][[g]][,1]})),
                       'g' = factor(rep(new_order[(1:params[['G']])], each = length(p_seq))))
ggplot(data = df_pcs, aes(x = pressure, y = pc, group = g, color = g, linetype = g)) + 
  geom_line(size = 1.1) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = label_use,
       linetype = 'Cluster', color = 'Cluster',
       title = 'First transformation functions')+
  theme(legend.position = 'bottom', text = element_text(size = 18))
ggsave(paste0(source_folder, 
              base_folder_beginning, 
              G, 
              base_folder_end, 'spatial_var_pcs_response1.png'),
       height = 8.86, width = 5.11)
label_use <- paste(var_name, 'PC (μmol/kg)')

df_pcs <- data.frame('pressure' = p_seq, 
                     'pc' = as.double(sapply((1:params[['G']]), function(g){ Phi_pred %*% parameters[['Omega_response2']][[g]][,1]})),
                     'g' = factor(rep(new_order[(1:params[['G']])], each = length(p_seq))))
ggplot(data = df_pcs, aes(x = pressure, y = pc, group = g, color = g, linetype = g)) + 
  geom_line(size = 1.1) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = label_use,
       linetype = 'Cluster', color = 'Cluster',
       title = paste('First', tolower(var_name),  'PC functions'))+
  theme(legend.position = 'bottom', text = element_text(size = 18))
ggsave(paste0(source_folder, 
              base_folder_beginning, 
              G, 
              base_folder_end, 'spatial_var_pcs_response2.png'),
       height = 8.86, width = 5.11)

df_pcs <- data.frame('pressure' = p_seq, 
                     'pc' = as.double(sapply((1:params[['G']]), function(g){ Phi_pred %*% parameters[['Omega_predictors']][[g]][1:nbasis_pred,1]})),
                     'g' = factor(rep(new_order[(1:params[['G']])], each = length(p_seq))))
ggplot(data = df_pcs, aes(x = pressure, y = pc, group = g, color = g, linetype = g)) + 
  geom_line(size = 1.1) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = 'Temperature PC (°C)',
       linetype = 'Cluster', color = 'Cluster',
       title = 'First temperature PC functions')+
  theme(legend.position = 'bottom', text = element_text(size = 18))
ggsave(paste0(source_folder, 
              base_folder_beginning, 
              G, 
              base_folder_end, 'spatial_var_pcs_predictors1.png'),
       height = 8.86, width = 5.11)

df_pcs <- data.frame('pressure' = p_seq, 
                     'pc' = as.double(sapply((1:params[['G']]), function(g){ Phi_pred %*% 
                         parameters[['Omega_predictors']][[g]][(nbasis_pred + 1):(2*nbasis_pred),1]})),
                     'g' = factor(rep(new_order[(1:params[['G']])], each = length(p_seq))))
ggplot(data = df_pcs, aes(x = pressure, y = pc, group = g, color = g, linetype = g)) + 
  geom_line(size = 1.1) + 
  coord_flip() + 
  scale_x_continuous(trans = 'reverse')+
  labs(x = 'Pressure (decibars)', y = 'Salinity PC (PSU)',
       linetype = 'Cluster', color = 'Cluster',
       title = 'First salinity PC functions')+
  theme(legend.position = 'bottom', text = element_text(size = 18))
ggsave(paste0(source_folder, 
              base_folder_beginning, 
              G, 
              base_folder_end, 'spatial_var_pcs_predictors2.png'),
       height = 8.86, width = 5.11)
likelihood_use <- estimated_model[['diagnostics']][['likelihood']]
likelihood_use[likelihood_use == Inf] <- (lead(likelihood_use) + lag(likelihood_use))[likelihood_use == Inf]/2
likelihood_df <- data.frame(iter = 0:(length(likelihood_use) - 1),
                            likelihood = likelihood_use)
ggplot(data = likelihood_df, aes(iter, likelihood))+
  geom_line() +
  labs(x = 'Iteration', y =  'Log Likelihood',
       title = 'Log likelihood trace plot by iteration')+
  theme(text = element_text(size = 18))
ggsave(paste0(source_folder, 
              base_folder_beginning, 
              G, 
              base_folder_end, 'log_likelihood.png'),
       height = 5.11, width = 8.86)


estimated_model[['diagnostics']][['measurement_error_predictors']]
estimated_model[['diagnostics']][['measurement_error_response']]
estimated_model[['parameters']][['variances_response']]
estimated_model[['parameters']][['range_params_r']]
estimated_model[['parameters']][['range_params_r']]

ggplot(data = rbind(data.frame(data_inputs[['locs']], mems = memberships_use, type=  data_inputs[['BGC']]))) +
  geom_point(aes(x = longitude, y = latitude), size = .02, alpha = .6) + facet_grid(~mems)  + 
  SO_coord + SO_theme + continents +
  labs(title = 'BGC profile locations, by cluster')
ggsave(paste0(source_folder, 
              base_folder_beginning, 
              G, 
              base_folder_end, 'spatial_clustering.png'),height = 4, width = 8)




