init_diagnostics <- function(model, parameters, params, time) {

  diagnostics <- list()
  
  diagnostics[['measurement_error_predictors']] = matrix(nrow = params[['EM_iter']] + 1,
                                                         ncol = 2)
  
  diagnostics[['measurement_error_response']] = vector(length = params[['EM_iter']] + 1,
                                                       'numeric')
  
  diagnostics[['time']] = vector(length = params[['EM_iter']]+1, 'numeric')
  diagnostics[['likelihood']] = vector(length = params[['EM_iter']]+1, 'numeric')
  diagnostics[['AIC']] = vector(length = params[['EM_iter']]+1, 'numeric')
  diagnostics[['BIC']] = vector(length = params[['EM_iter']]+1, 'numeric')
  diagnostics[['penalty']] = vector(length = params[['EM_iter']]+1, 'numeric')
  diagnostics[['ESS']] = vector(length = params[['EM_iter']]+1, 'numeric')
  
  diagnostics[['lambda_mean_response']] = vector(length = params[['EM_iter']], 'numeric')
  if (class(model) != 'one_var'){
    diagnostics[['lambda_mean_predictors']] = vector(length = params[['EM_iter']], 'numeric')
  }
    
  diagnostics
}

update_diagnostics <- function(model, i_EM, parameters, params, diagnostics, time,
                               MCEM_res, data_inputs) {
  
  diagnostics[['lambda_mean_response']][i_EM+1] = parameters[['lambda_mean_response']]
  if (class(model) != 'one_var'){
    diagnostics[['lambda_mean_predictors']][i_EM+1] = parameters[['lambda_mean_predictors']]
  }
  
  if (class(model) == 'mult_var_spat'){
    diagnostics[['penalty']][i_EM+1] =
      sum(sapply(1:params[['G']], function(g) {
        parameters[['lambda_mean_response']] * as.double(t(parameters[['means_response']][[g]]) %*%
                                                             (data_inputs[['penalty_mat_response']] %*% parameters[['means_response']][[g]]))
      })) + 
      sum(sapply(1:params[['G']], function(g) {
        sapply(1:params[['pc_response2']], function(q) {
          parameters[['lambda_pcs_response']] * as.double(t(parameters[['Omega_response2']][[g]][,q]) %*%
                                                           (data_inputs[['penalty_mat_response']] %*% parameters[['Omega_response2']][[g]][,q]))
        })
      })) + 
      sum(sapply(1:params[['G']], function(g) {
        parameters[['lambda_mean_predictors']] * as.double(t(parameters[['means_predictors']][[g]]) %*%
                                                               (data_inputs[['penalty_mat_predictors']] %*% parameters[['means_predictors']][[g]]))
      })) + 
      sum(sapply(1:params[['G']], function(g) {
        sapply(1:params[['pc_predictors']], function(q) {
          parameters[['lambda_pcs_predictors']] * as.double(t(parameters[['Omega_predictors']][[g]][,q]) %*%
                                                            (data_inputs[['penalty_mat_predictors']] %*% parameters[['Omega_predictors']][[g]][,q]))
        })
      }))
  }
  model = compute_information_criteria(model, params, data_inputs, parameters, MCEM_res=MCEM_res)
  
  diagnostics[['AIC']][i_EM+1] = model[['AIC']]
  diagnostics[['BIC']][i_EM+1] = model[['BIC']]
  diagnostics[['likelihood']][i_EM+1] = model[['likelihood']]
  
  if (params[['EM_iter']] > i_EM){
    diagnostics[['ESS']][i_EM+1] = params[['MC_max']]^2 / sum(exp(2 * log(MCEM_res[['weights']])))
  }

  diagnostics[['measurement_error_response']][i_EM + 1] = parameters[['measurement_error_response']]
  
  if(class(model) != 'one_var'){
    diagnostics[['measurement_error_predictors']][i_EM + 1,] = parameters[['measurement_error_predictors']]
  }
  
  diagnostics[['time']][i_EM+1] = as.numeric(time)
  
  diagnostics
}
