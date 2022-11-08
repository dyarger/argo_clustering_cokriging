library(Matrix)
library(fstmr)
array_id <- as.integer(Sys.getenv('array_id'))
variable <- Sys.getenv('var')
G <- Sys.getenv('G')
pc_response <- Sys.getenv('pc_response')
m_train_response <- Sys.getenv('m_train_response')
m_train_predictors <- Sys.getenv('m_train_predictors')
m_AIC_response <- Sys.getenv('m_AIC_response')
m_AIC_predictors <- Sys.getenv('m_AIC_predictors')
MC_max_AIC <- Sys.getenv('MC_max_AIC')

combinations <- data.frame(expand.grid('n_pcs_predictors' = c(5, 8:13), 
                                       'n_pcs_response' = c(5:12)))

base_folder <- Sys.getenv('base_folder')
n_pcs_start <-combinations[['n_pcs_predictors']][array_id]
new_number <- combinations[['n_pcs_response']][array_id]

results_folder <- paste0(base_folder, variable, '_none_NA_', G, '_', n_pcs_start, '_',
                         pc_response, '_', m_train_predictors, '_', m_train_response, '/')
load(paste0(results_folder, 'final_results.RData'))

# delete response pcs
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

if (is.na(MC_max_AIC) | is.null(MC_max_AIC)) {
  MC_max_AIC <- estimated_model[['params']][['MC_max']]
}
estimated_model[['params']][['MC_max']] <- MC_max_AIC
model <- fstmr:::compute_information_criteria(model = estimated_model, 
                                              params = estimated_model[['params']], 
                                              data_inputs = estimated_model[['data_inputs']], 
                                              parameters = estimated_model[['parameters']])
model[['data_inputs']] <- NULL
model[['parameters']] <- NULL
save(model, file = paste0(results_folder, 'AIC_', new_number, '.RData'))
