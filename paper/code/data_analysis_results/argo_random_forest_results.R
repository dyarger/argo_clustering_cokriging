
.libPaths("/Library/Frameworks/R.framework/Versions/4.1/Resources/miworkspace-library")
library(tidyverse)
variable <- Sys.getenv('variable')
base_folder_beginning <- Sys.getenv('base_folder_beginning')
base_folder_end <- Sys.getenv('base_folder_end')
source_folder <- Sys.getenv('source_folder')
n_folds <- 5
type = 'profiles'

rf_list <- list()
for (i in 1:n_folds) {
  load(paste0(source_folder, 'randomForest/', type, '_', variable, '_', n_folds, '_',  i, '_rf.RData'))
  load(paste0(source_folder, variable,'_', type, '_',  i, base_folder_end, 'prediction.RData'))
  our_pred <- pred_at_p_levels$pred_df %>%
    group_by(index) %>%
    summarise(avg_full = mean(pred_full*weight)) %>% pull(avg_full) 
  our_results <- cbind(held_out_profs, our_pred = our_pred) %>%
    filter(pressure < 155, pressure > 145) %>%
    mutate(fold = i, our_residual = oxy - our_pred) %>%
    dplyr::select(fold, profile_unique, pressure, oxy, our_pred, our_residual)
  
  results_combined <- full_join(our_results, rf_results)
  rf_list[[i]] <- results_combined
  print(i)
}
rf_results <- bind_rows(rf_list)
summary(rf_results$residual)
sqrt(mean(rf_results$residual^2))
median(abs(rf_results$residual))

sqrt(mean(rf_results$our_residual^2))
median(abs(rf_results$our_residual))

n_folds <- 10
type = 'floats'

rf_list <- list()
for (i in 1:n_folds) {
  load(paste0(source_folder, 'randomForest/', type, '_', variable, '_', n_folds, '_',  i, '_rf.RData'))
  load(paste0(source_folder, variable,'_', type, '_',  i, base_folder_end, 'prediction.RData'))
  our_pred <- pred_at_p_levels$pred_df %>%
    group_by(index) %>%
    summarise(avg_full = mean(pred_full*weight)) %>% pull(avg_full)
  our_results <- cbind(held_out_profs, our_pred = our_pred) %>%
    filter(pressure < 155, pressure > 145) %>%
    mutate(fold = i, our_residual = oxy - our_pred) %>%
    dplyr::select(fold, profile_unique, pressure, oxy, our_pred, our_residual)
  results_combined <- full_join(our_results, rf_results)
  
  rf_list[[i]] <- results_combined
}
rf_results <- bind_rows(rf_list)
summary(rf_results$residual)
sqrt(mean(rf_results$residual^2))
median(abs(rf_results$residual))

sqrt(mean(rf_results$our_residual^2))
median(abs(rf_results$our_residual))

