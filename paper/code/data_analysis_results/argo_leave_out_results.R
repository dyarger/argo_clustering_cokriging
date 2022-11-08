library(dplyr)
library(ggplot2)
library(Matrix)
library(patchwork)
library(fda)
variable <- Sys.getenv('variable')
base_folder_end <- Sys.getenv('base_folder_end')
source_folder <- Sys.getenv('source_folder')
var_label <- 'Oxygen'
# new leave out results
df_list <- list()
for (i in 1:5) {
  load(paste0('paper/results/', variable, '_profiles_', i, base_folder_end, 'prediction.RData'))
  pred_df <- pred_at_p_levels$pred_df
  df_list[[i]] <- pred_df %>%
    group_by(index) %>%
    summarise(avg_full = mean(pred_full*weight),
              avg_TS = mean(pred_TS*weight),
              avg_mean = mean(pred_mean*weight),
              avg_var = mean(var_at_pressure*weight),
              var_avg = mean((pred_full - mean(pred_full*weight))^2 * weight),
              profile_unique = profile_unique[1], 
              type = 'profile') %>%
    left_join(cbind(held_out_profs, index = 1:nrow(held_out_profs)))
}

df_list_floats <- list()
for (i in 1:10) {
  load(paste0('paper/results/', variable, '_floats_', i, base_folder_end, 'prediction.RData'))
  pred_df <- pred_at_p_levels$pred_df
  df_list_floats[[i]] <- pred_df %>%
    group_by(index) %>%
    summarise(avg_full = mean(pred_full*weight),
              avg_TS = mean(pred_TS*weight),
              avg_mean = mean(pred_mean*weight),
              avg_var = mean(var_at_pressure*weight),
              var_avg = mean((pred_full - mean(pred_full*weight))^2 * weight),
              profile_unique = profile_unique[1], 
              type = 'float') %>%
    left_join(cbind(held_out_profs, index = 1:nrow(held_out_profs)))
}

mean_prediction <- rbind(dplyr::bind_rows(df_list), bind_rows(df_list_floats))
pressure_summary <- mean_prediction %>%
  mutate(pgroup = (seq(0, 2000, by = 100) + 50)[findInterval(pressure, seq(0, 2000, by = 100))]) %>%
  group_by(pgroup, type) %>%
  summarise(med_abs = median(abs(oxy - avg_full)),
            rmse = sqrt(mean((oxy - avg_full)^2)),
            med_abs_mean = median(abs(oxy - avg_mean)),
            rmse_mean = sqrt(mean((oxy - avg_mean)^2)),
            med_abs_TS = median(abs(oxy - avg_TS)),
            rmse_TS = sqrt(mean((oxy - avg_TS)^2)),
            cov_1 = mean(abs(oxy - avg_full) < sqrt(avg_var + var_avg)),
            cov_2 = mean(abs(oxy - avg_full) < 2*sqrt(avg_var + var_avg)),
            cov_3 = mean(abs(oxy - avg_full) < 3*sqrt(avg_var + var_avg)),
            med_length = median(sqrt(avg_var + var_avg)))


df_label <- data.frame('Type' = c('Mean', 'Full\nprediction', 'Temperature\nand salinity'), 
                       'type_var' = c('med_abs_mean', 'med_abs', 'med_abs_TS', 'rmse_mean', 'rmse', 'rmse_TS'))
df_model_label <- data.frame('ModelType' = c('Profiles', 'Floats'), 
                             'type' = c('profile', 'float'))
df_cov_label <- data.frame('type_var' = paste0('cov_', 1:3), 
                           'type_label' = paste0('Within ', 1:3, ' SD', c('', 's', 's')))

theme_set(theme_bw())


a <- ggplot(data = tidyr::pivot_longer(pressure_summary, cols = starts_with('med_abs'),
                                       names_to = 'type_var', values_to = 'value') %>%
              left_join(df_label) %>% left_join(df_model_label) ,
            aes(x = pgroup, shape = Type,color = Type, y = value,
                linetype = Type)) + 
  geom_point() + 
  geom_line() + 
  scale_y_continuous(limits = c(0, 20)) + 
  labs(x = 'Pressure (in 100 dbar bins)', y = paste(var_label, 'Median Absolute Error\n(μmol/kg)'),
       linetype = 'Prediction\nType', shape = 'Prediction\nType', color = 'Prediction\nType',
       title = 'Leave out point prediction performance')+
  guides(linetype = guide_legend(),
         color =guide_legend(), shape = guide_legend()) + 
  theme(text = element_text(size = 14)) + 
  facet_wrap(~ModelType, ncol = 2)
a
ggsave(paste0(source_folder, variable, '_', 
              'profiles', '_',
              1, base_folder_end, 'spatial_error_pressure.png'), plot = a, height = 4,# height = 3, 
       width = 12)


a <- ggplot(data = tidyr::pivot_longer(pressure_summary, cols = starts_with('cov_'),
                                       names_to = 'type_var', values_to = 'value') %>%
              left_join(df_cov_label)%>% left_join(df_model_label),
            aes(x = pgroup, color = type_label, y = value, group = type_label,
                linetype = type_label, shape = type_label
            )
) + 
  geom_point() + 
  geom_line() + 
  geom_hline(data = data.frame('type_var' = paste0('cov_', 1:3), 
                               value = 1 - (1-pnorm(1:3))*2, 
                               amount_label = paste0('Within ', 1:3, 
                                                     ' SD', c('', 's', 's'))) %>%
               left_join(df_cov_label),
             aes(yintercept = value))+
  scale_y_continuous(limits = c(.5, NA)) + 
  facet_wrap(~ModelType, ncol = 2)+
  labs(x = 'Pressure (in 100 dbar bins)', y = 'Coverage\n',
       title = 'Prediction coverages',
       linetype = 'Interval\nLength', shape = 'Interval\nLength', color = 'Interval\nLength')+
  theme(text = element_text(size = 14))
a
b <- ggplot(data = tidyr::pivot_longer(pressure_summary, cols = starts_with('med_length'),
                                       names_to = 'type_var', values_to = 'value') %>%
              left_join(df_cov_label)%>% left_join(df_model_label),
            aes(x = pgroup, y = value, group = type, color = ModelType, linetype = ModelType, 
                shape = ModelType)) + 
  geom_point(size = .7) + 
  geom_line() + 
  scale_y_continuous(limits = c(0, 12)) + 
  labs(x = 'Pressure (in 100 dbar bins)', 
       y = 'Median Prediction Standard\nDeviation (μmol/kg)',
       title = 'Size of predicted uncertainties',
       color = 'Leave out\ntype',
       linetype = 'Leave out\ntype',
       shape = 'Leave out\ntype')
b
patchwork::wrap_plots(list(a,b), ncol = 2, widths = c(2,1))
ggsave(paste0(source_folder, variable, '_', 
              'profiles', '_',
              1, base_folder_end, 'spatial_uncertainty_pressure.png'), 
       height = 3*1.5, width =  5.8*2)

ggsave(plot = a, paste0(source_folder, variable, '_', 
                        'profiles', '_',
                        1, base_folder_end, 'spatial_coverage_pressure.png'), 
       height = 4, width = 12)

ggsave(plot = b, paste0(source_folder, variable, '_', 
                        'profiles', '_',
                        1, base_folder_end, 'spatial_uncertainty_pressure_lens.png'), 
       height = 3*1.5, width = 5.8*2)
