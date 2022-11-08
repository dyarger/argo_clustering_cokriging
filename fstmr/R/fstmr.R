#' This is nice.
#' 
#' @param data List of data.frames or data.frame containing the data
#' @param params Hyperparameters for the EM algorithm
#' @param diagnostics Should diagnostics such as AIC be computed?
#' @param verbose Provide some info during fitting
#' @return A model object
#' @export
fstmr <- function(data, params, compute_diagnostics = F, verbose = F, printf=NULL){
  
  # Set necessary hyper-parameters
  params = init_params(params, compute_diagnostics)

  if(verbose > 2) print("Initialized hyperparameters.")
  
  # Prepare EM algorithm
  data_inputs = prepare_EM_data(data, params)
  if(verbose > 1) print("Data inputs prepared.")
  
  # Check which model
  model = determine_model(params)
  if(verbose > 0) print(paste0("Determined model to be: ", class(model), "."))
  
  parameters = initialize(model, params, data_inputs)
  if(verbose > 0) print("Initialized parameters.")
  
  if(compute_diagnostics) diagnostics = init_diagnostics(model, parameters, params, Sys.time())

  for (i in 1:params[['EM_iter']]){
    MCEM_res = Estep(model, params, data_inputs, parameters)
    if(verbose > 1) print(paste0("Estep done in iteration ", as.character(i), "."))
    
    if(compute_diagnostics) diagnostics = update_diagnostics(model, i-1, parameters, params, diagnostics,
                                                     Sys.time(), MCEM_res, data_inputs)
    
    if(!is.null(printf) & class(printf) == 'function'){
      printf(model, diagnostics, parameters, data_inputs, params)
    }
    
    parameters = Mstep(model, params, data_inputs, parameters, MCEM_res, i)
    if(verbose > 1) print(paste0("Mstep done in iteration ", as.character(i), "."))
  }
  
  if(verbose > 0) print("EM algorithm done.")
  
  if(compute_diagnostics) {
    diagnostics = update_diagnostics(model, params[['EM_iter']], parameters, params, diagnostics,
                                     Sys.time(), NULL, data_inputs)
    model[['diagnostics']] = diagnostics
    model[['AIC']] = diagnostics[['AIC']]
    model[['BIC']] = diagnostics[['BIC']]
    model[['likelihood']] = diagnostics[['likelihood']]
    model[['diagnostics']] = diagnostics
  } else {
    model = compute_information_criteria(model, params, data_inputs, parameters)
  }
  
  if(verbose > 1) print("Computed information criteria.")
  
  model$parameters = parameters
  model$data_inputs = data_inputs
  model$params = params
  return(model)
}

init_params <- function(params, compute_diagnostics=F){
  if(compute_diagnostics | length(params[['lambda_mean_response']]) > 1){
    params[['compute_likelihood']] = T
  }
  if(is.null(params$G)){
    params$G = 1
  }
  if(is.null(params$loc_names)){
    params$loc_names = c('longitude', 'latitude')
  }
  if(is.null(params$init_strategy)){
    params[['init_strategy']] = 'kmeans'
  }
  if(is.null(params[['id']])){
    params[['id']] = 'profile_unique'
  }
  if(is.null(params[['ind_name']])){
    params[['ind_name']] = 'pressure'
  }
  if(is.null(params[['use_MRF']])){
    params[['use_MRF']] = T
  }
  if(is.null(params[['cluster_sampling']])){
    if (is.null(params[['covariance_function']])){
      params[['cluster_sampling']] = 'independence'
    } else {
      params[['cluster_sampling']] = 'importance_sampling'
    }
  }
  if(is.null(params[['cv_skip']])){
    params[['cv_skip']] = 20
  }
  return(params)
}

determine_model <- function(params){
  model = structure(list(), class = 'mult_var_spat')
  if(is.null(params[['X']])){
    class(model) = 'one_var'
  } else if (is.null(params[['covariance_function']])){
    class(model) = 'mult_var_ind'
  }
  return(model)
}




