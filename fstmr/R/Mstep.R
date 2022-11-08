Mstep <- function(model, params, data_inputs, parameters, MCEM_res, iter){
  UseMethod('Mstep')
}

Mstep.one_var <- function(model, params, data_inputs, parameters, MCEM_res, iter){
  
  # To improve the readability of the code
  weights = MCEM_res[['weights']]
  cluster_mat = MCEM_res[['cluster_mat']]
  response_pcs_array = MCEM_res[['pcs_array_response']]
  n_profiles = data_inputs[['n_profiles']]
  G = params[['G']]
  
  cmi = cluster_mat - 1
  mode(cmi) = 'integer'
  Linv_mats <- MCEM_res[['Linv_mat']]
  vecchia_info <- MCEM_res[['vecchia_info']]
  
  # Update all the parameters
  u_mats = compute_u_matrices(Phi_x_Phi_resp = data_inputs[['Phi_x_Phi_response']], 
                              weights = weights,
                              cmi = NULL,
                              cmi_BGC = cmi,
                              n_e_step_samples=ncol(cmi),
                              G = G,
                              pcs_r = response_pcs_array)
  
  if (length(params[['lambda_mean_response']]) == 1 | (iter-1) %% params[['cv_skip']] != 0){
    parameters[['means_response']] = update_means(U = u_mats[['mean_resp']],
                                                  Phi_prof = data_inputs[['Phi_response_prof']],
                                                  profs = data_inputs[['response_profs']],
                                                  cmi = cmi,
                                                  pcs_array1 = response_pcs_array,
                                                  Omegas1 = parameters[['Omega_response']],
                                                  weights = weights,
                                                  pen_mat = data_inputs[['penalty_mat_response']],
                                                  G = G,
                                                  lambda = parameters[['lambda_mean_response']])
  } else {
    AIC = Inf
    for (lambda in params[['lambda_mean_response']]){
      old_means = parameters[['means_response']]; old_lambda = parameters[['lambda_mean_response']]
      parameters[['means_response']] = update_means(U = u_mats[['mean_resp']],
                                                    Phi_prof = data_inputs[['Phi_response_prof']],
                                                    profs = data_inputs[['response_profs']],
                                                    cmi = cmi,
                                                    pcs_array1 = response_pcs_array,
                                                    Omegas1 = parameters[['Omega_response']],
                                                    weights = weights,
                                                    pen_mat = data_inputs[['penalty_mat_response']],
                                                    G = G,
                                                    lambda = lambda)
      parameters[['lambda_mean_response']] = lambda
      model = compute_information_criteria(model, params, data_inputs, parameters, u_mats)
      if(model$AIC > AIC){
        parameters[['lambda_mean_response']] = old_lambda; parameters[['means_response']] = old_means
      }
      AIC = ifelse(model$AIC > AIC, AIC, model$AIC)
    }
  }
  
  if (length(params[['lambda_pcs_response']]) == 1 | (iter-1) %% params[['cv_skip']] != 0){
    parameters[['Omega_response']] = update_principal_components(u_mats = u_mats[['pc_resp']],
                                                                 cmi = cmi,
                                                                 pcs_r1 = response_pcs_array,
                                                                 profs = data_inputs[['response_profs']],
                                                                 Phi_prof = data_inputs[['Phi_response_prof']],
                                                                 Omegas1 = parameters[['Omega_response']],
                                                                 means = parameters[['means_response']], 
                                                                 weights = rep(1, ncol(cluster_mat)),
                                                                 pen_mat = data_inputs[['penalty_mat_response']],
                                                                 G = G,
                                                                 lambda = parameters[['lambda_pcs_response']])
  } else {
    AIC = Inf
    for (lambda in params[['lambda_pcs_response']]){
      old_pcs = parameters[['Omega_response']]; old_lambda = parameters[['lambda_pcs_response']]
      parameters[['Omega_response']] = update_principal_components(u_mats = u_mats[['pc_resp']],
                                                                   cmi = cmi,
                                                                   pcs_r1 = response_pcs_array,
                                                                   profs = data_inputs[['response_profs']],
                                                                   Phi_prof = data_inputs[['Phi_response_prof']],
                                                                   Omegas1 = parameters[['Omega_response']],
                                                                   means = parameters[['means_response']], 
                                                                   weights = rep(1, ncol(cluster_mat)),
                                                                   pen_mat = data_inputs[['penalty_mat_response']],
                                                                   G = G,
                                                                   lambda = lambda)
      parameters[['lambda_pcs_response']] = lambda
      model = compute_information_criteria(model, params, data_inputs, parameters, u_mats)
      if(model$AIC > AIC){
        parameters[['lambda_pcs_response']] = old_lambda; parameters[['Omega_response']] = old_pcs
      }
      AIC = ifelse(model$AIC > AIC, AIC, model$AIC)
    }
  }
  
  parameters[['measurement_error_response']] = update_measurement_error(cluster_mat = cluster_mat,
                                                                        pcs_array1 = response_pcs_array,
                                                                        profs = data_inputs$response_profs,
                                                                        phi_prof =  data_inputs$Phi_response_prof,
                                                                        means = parameters$means_response,
                                                                        Omegas1 = parameters$Omega_response, 
                                                                        n_profiles = data_inputs$n_profiles, 
                                                                        weights = weights, 
                                                                        G = G,
                                                                        prof_ns = data_inputs[['profile_lengths_response']])

  parameters[['variances_response']] = update_variances(cluster_mat, response_pcs_array, parameters[['range_params']],
                                                           params$G, dim(response_pcs_array)[3],
                                                           ncol(cluster_mat), 
                                                           params[['covariance_function']],
                                                           independent_model = is.null(params[['covariance_function']]),
                                                           Linv = Linv_mats, variances = parameters[['variances_response']], 
                                                           vecchia_info =  vecchia_info)
  if (!is.null(params[['covariance_function']])) {
    parameters[['range_params']] = update_range_params(G = params[['G']],
                                                       n_pcs = params[['pc_response']],
                                                       variances = parameters[['variances_response']],
                                                       cluster_mat = cluster_mat,
                                                       pcs_array = response_pcs_array,
                                                       vecchia_info = vecchia_info,
                                                       mc_weights = weights, 
                                                       covfun_name = params[['covariance_function']])
  }
  
  # Sort
  for(g in 1:params$G){
    inds = order(parameters$variances_response[[g]], decreasing=T)
    parameters$Omega_response[[g]] = parameters$Omega_response[[g]][,inds]
    parameters$variances_response[[g]] = parameters$variances_response[[g]][inds]
    parameters$range_params[[g]] = parameters$range_params[[g]][,inds, drop = F]
  }
    
  # Orthogonality constraints
  ortho_info <- orthogonalize(Omegas1 = parameters[['Omega_response']], 
                              variances1 = parameters[['variances_response']],
                              inner_prod_mat1 = data_inputs[['inner_prod_mat_response']],
                              G = G)
  
  parameters[['Omega_response']] <- lapply(1:G, function(g) {ortho_info[[g]][[2]]})
  parameters[['variances_response']] <- lapply(1:G, function(g) {ortho_info[[g]][[1]]})
  
  cluster_list = update_cluster_membership(cluster_mat, n_profiles, G)
  parameters[['cluster_membership']] = cluster_list[[1]]
  parameters[['conditional_cluster_probs']] = cluster_list[[2]]
  parameters[['theta']] = update_theta(c(0,3), data_inputs[['nn_list']], 
                                       cluster_mat = cluster_mat,
                                       n_profiles = n_profiles, G,
                                       weights = weights,
                                       theta = parameters[['theta']])
  
  parameters[['marginal_probabilities']] <- compute_marginal_probabilities(parameters[['cluster_membership']], parameters[['theta']], neighbors = data_inputs[['nn_list']], n_profiles, G)
  
  return(parameters)
}


Mstep.mult_var_ind <- function(model, params, data_inputs, parameters, MCEM_res, iter){
  # To improve the readability of the code
  cluster_mat = MCEM_res[[1]]
  response_pcs_array = MCEM_res[[2]]
  predictor_pcs_array = MCEM_res[[3]]
  n_profiles = data_inputs[['n_profiles']]
  n_profiles_TS = data_inputs[['n_profiles_TS']]
  G = params[['G']]
  cluster_mat_BGC <- cluster_mat[data_inputs[['BGC']], ,drop = F]
  predictor_pcs_array_BGC <- predictor_pcs_array[data_inputs[['BGC']], , , drop = F]

  cmi = cluster_mat -1 
  mode(cmi) = 'integer'
  cmi_BGC = cluster_mat_BGC -1 
  mode(cmi_BGC) = 'integer'
  
  weights = rep(1, ncol(cmi))
  # Update all the parameters
  
  u_mats = compute_u_matrices(Phi_x_Phi_resp = data_inputs[['Phi_x_Phi_response']], 
                              weights = weights,
                              cmi = cmi,
                              cmi_BGC = cmi_BGC,
                              n_e_step_samples=ncol(cmi),
                              G = G,
                              pcs_r = response_pcs_array,
                              Phi_x_Phi_pred = data_inputs[['Phi_x_Phi_predictors']],
                              pcs_p1 = predictor_pcs_array,
                              pcs_p2 = predictor_pcs_array_BGC,
                              ind = T)
  
  if (length(params[['lambda_mean_response']]) == 1| (iter-1) %% params[['cv_skip']] != 0){
    parameters[['means_response']] = update_means(U = u_mats[['mean_resp']],
                                                  Phi_prof = data_inputs[['Phi_response_prof']],
                                                  profs = data_inputs[['response_profs']],
                                                  cmi = cmi_BGC,
                                                  pcs_array1 = response_pcs_array,
                                                  Omegas1 = parameters[['Omega_response']],
                                                  weights = weights,
                                                  pen_mat = data_inputs[['penalty_mat_response']],
                                                  G = G,
                                                  lambda = parameters[['lambda_mean_response']])
  } else {
    AIC = Inf
    for (lambda in params[['lambda_mean_response']]){
      old_means = parameters[['means_response']]; old_lambda = parameters[['lambda_mean_response']]
      parameters[['means_response']] = update_means(U = u_mats[['mean_resp']],
                                                    Phi_prof = data_inputs[['Phi_response_prof']],
                                                    profs = data_inputs[['response_profs']],
                                                    cmi = cmi_BGC,
                                                    pcs_array1 = response_pcs_array,
                                                    Omegas1 = parameters[['Omega_response']],
                                                    weights = weights,
                                                    pen_mat = data_inputs[['penalty_mat_response']],
                                                    G = G,
                                                    lambda = lambda)
      parameters[['lambda_mean_response']] = lambda
      model = compute_information_criteria(model, params, data_inputs, parameters, u_mats)
      if(model$AIC > AIC){
        parameters[['lambda_mean_response']] = old_lambda; parameters[['means_response']] = old_means
      }
      AIC = ifelse(model$AIC > AIC, AIC, model$AIC)
    }
  }
  
  if (length(params[['lambda_mean_predictors']]) == 1| (iter-1) %% params[['cv_skip']] != 0){
    parameters[['means_predictors']] = update_means(U = u_mats[['mean_pred']],
                                                    Phi_prof = data_inputs[['Phi_predictors_prof']],
                                                    profs = data_inputs[['predictor_profs']],
                                                    cmi = cmi,
                                                    pcs_array1 = predictor_pcs_array,
                                                    Omegas1 = parameters[['Omega_predictors']],
                                                    weights = weights,
                                                    pen_mat = data_inputs[['penalty_mat_predictors']],
                                                    G = G,
                                                    lambda = parameters[['lambda_mean_predictors']])
  } else {
    AIC = Inf
    for (lambda in params[['lambda_mean_predictors']]){
      old_means = parameters[['means_predictors']]; old_lambda = parameters[['lambda_mean_predictors']]
      parameters[['means_predictors']] = update_means(U = u_mats[['mean_pred']],
                                                      Phi_prof = data_inputs[['Phi_predictors_prof']],
                                                      profs = data_inputs[['predictor_profs']],
                                                      cmi = cmi,
                                                      pcs_array1 = predictor_pcs_array,
                                                      Omegas1 = parameters[['Omega_predictors']],
                                                      weights = weights,
                                                      pen_mat = data_inputs[['penalty_mat_predictors']],
                                                      G = G,
                                                      lambda = lambda)
      parameters[['lambda_mean_predictors']] = lambda
      model = compute_information_criteria(model, params, data_inputs, parameters, u_mats)
      if(model$AIC > AIC){
        parameters[['lambda_mean_predictors']] = old_lambda; parameters[['means_predictors']] = old_means
      }
      AIC = ifelse(model$AIC > AIC, AIC, model$AIC)
    }
  }
  
  if (length(params[['lambda_pcs_response']]) == 1| (iter-1) %% params[['cv_skip']] != 0){
    parameters[['Omega_response']] = update_principal_components(u_mats = u_mats[['pc_resp']],
                                                                  cmi = cmi_BGC,
                                                                  pcs_r1 = response_pcs_array,
                                                                  profs = data_inputs[['response_profs']],
                                                                  Phi_prof = data_inputs[['Phi_response_prof']],
                                                                  Omegas1 = parameters[['Omega_response']],
                                                                  means = parameters[['means_response']], 
                                                                  weights = rep(1, ncol(cluster_mat_BGC)),
                                                                  pen_mat = data_inputs[['penalty_mat_response']],
                                                                  G = G,
                                                                  lambda = parameters[['lambda_pcs_response']])
  } else {
    AIC = Inf
    for (lambda in params[['lambda_pcs_response']]){
      old_pcs = parameters[['Omega_response']]; old_lambda = parameters[['lambda_pcs_response']]
      parameters[['Omega_response']] = update_principal_components(u_mats = u_mats[['pc_resp']],
                                                                   cmi = cmi_BGC,
                                                                   pcs_r1 = response_pcs_array,
                                                                   profs = data_inputs[['response_profs']],
                                                                   Phi_prof = data_inputs[['Phi_response_prof']],
                                                                   Omegas1 = parameters[['Omega_response']],
                                                                   means = parameters[['means_response']], 
                                                                   weights = rep(1, ncol(cluster_mat_BGC)),
                                                                   pen_mat = data_inputs[['penalty_mat_response']],
                                                                   G = G,
                                                                   lambda = lambda)
      parameters[['lambda_pcs_response']] = lambda
      model = compute_information_criteria(model, params, data_inputs, parameters, u_mats)
      if(model$AIC > AIC){
        parameters[['lambda_pcs_response']] = old_lambda; parameters[['Omega_response']] = old_pcs
      }
      AIC = ifelse(model$AIC > AIC, AIC, model$AIC)
    }
  }
  
  if (length(params[['lambda_pcs_predictors']]) == 1| (iter-1) %% params[['cv_skip']] != 0){
    parameters[['Omega_predictors']] = update_principal_components(u_mats = u_mats[['pc_pred']],
                                                                   cmi = cmi,
                                                                   pcs_r1 = predictor_pcs_array,
                                                                   profs = data_inputs[['predictor_profs']],
                                                                   Phi_prof = data_inputs[['Phi_predictors_prof']],
                                                                   Omegas1 = parameters[['Omega_predictors']],
                                                                   means = parameters[['means_predictors']], 
                                                                   weights = rep(1, ncol(cluster_mat)),
                                                                   pen_mat = data_inputs[['penalty_mat_predictors']],
                                                                   G = G,
                                                                   lambda = parameters[['lambda_pcs_predictors']])
  } else {
    AIC = Inf
    for (lambda in params[['lambda_pcs_predictors']]){
      old_pcs = parameters[['Omega_predictors']]; old_lambda = parameters[['lambda_pcs_predictors']]
      parameters[['Omega_predictors']] = update_principal_components(u_mats = u_mats[['pc_pred']],
                                                                     cmi = cmi,
                                                                     pcs_r1 = predictor_pcs_array,
                                                                     profs = data_inputs[['predictor_profs']],
                                                                     Phi_prof = data_inputs[['Phi_predictors_prof']],
                                                                     Omegas1 = parameters[['Omega_predictors']],
                                                                     means = parameters[['means_predictors']], 
                                                                     weights = rep(1, ncol(cluster_mat)),
                                                                     pen_mat = data_inputs[['penalty_mat_predictors']],
                                                                     G = G,
                                                                     lambda = lambda)
      parameters[['lambda_pcs_predictors']] = lambda
      model = compute_information_criteria(model, params, data_inputs, parameters, u_mats)
      if(model$AIC > AIC){
        parameters[['lambda_pcs_predictors']] = old_lambda; parameters[['Omega_predictors']] = old_pcs
      }
      AIC = ifelse(model$AIC > AIC, AIC, model$AIC)
    }
  }
  
  parameters[['Lambdas']] = update_Lambda(cluster_mat_BGC, response_pcs_array, predictor_pcs_array_BGC, G)
  
  parameters[['measurement_error_response']] = update_measurement_error(cluster_mat = cluster_mat_BGC,
                                                                        pcs_array1 = response_pcs_array,
                                                                        profs = data_inputs[['response_profs']],
                                                                        phi_prof = data_inputs[['Phi_response_prof']],
                                                                        means = parameters[['means_response']],
                                                                        Omegas1 = parameters[['Omega_response']],
                                                                        n_profiles = n_profiles, 
                                                                        G = G,
                                                                        prof_ns = data_inputs[['profile_lengths_response']])
  
  parameters[['measurement_error_predictors']] = update_measurement_error(cluster_mat = cluster_mat,
                                                                          pcs_array1 = predictor_pcs_array,
                                                                          profs = data_inputs[['predictor_profs']],
                                                                          phi_prof = data_inputs[['Phi_predictors_prof']],
                                                                          means = parameters[['means_predictors']],
                                                                          Omegas1 = parameters[['Omega_predictors']],
                                                                          n_profiles = n_profiles_TS,
                                                                          G = G,
                                                                          n_vars = length(params[['X']]),
                                                                          prof_ns = data_inputs[['profile_lengths_predictors']])
  
  
  cluster_list = update_cluster_membership(cluster_mat, n_profiles_TS, G)
  parameters[['cluster_membership']] = cluster_list[[1]]
  parameters[['conditional_cluster_probs']] = cluster_list[[2]]
  parameters[['theta']] = update_theta(c(0,3), data_inputs[['nn_list']], 
                                       cluster_mat = cluster_mat,
                                       n_profiles = n_profiles_TS, G,
                                       weights = weights,
                                       theta = parameters[['theta']])

  orthogonalize_list = orthogonalize(Omegas1 = parameters[['Omega_response']],
                                     inner_prod_mat1 = data_inputs[['inner_prod_mat_response']],
                                     G = G,
                                     Omegas2 =  parameters[['Omega_predictors']],
                                     inner_prod_mat2 = data_inputs[['inner_prod_mat_predictors']],
                                     Lambdas = parameters[['Lambdas']],
                                     cluster_mat = cluster_mat_BGC,
                                     pcs1 = response_pcs_array, 
                                     pcs2 = predictor_pcs_array_BGC)
  
  parameters[['Omega_response']] = orthogonalize_list[[1]]
  parameters[['Omega_predictors']] = orthogonalize_list[[2]]
  parameters[['variances_response']] = orthogonalize_list[[3]]
  parameters[['variances_predictors']] = orthogonalize_list[[4]]
  parameters[['Lambdas']] = orthogonalize_list[[5]]
  parameters[['Sigma_eta_inv']] = orthogonalize_list[[6]]
  
  parameters[['marginal_probabilities']] <- compute_marginal_probabilities(parameters[['cluster_mat']], parameters[['theta']], neighbors = data_inputs[['nn_list']], n_profiles, G)
  
  return(parameters)
}

Mstep.mult_var_spat <- function(model, params, data_inputs, parameters, MCEM_res, iter){
  # To improve the readability of the code
  weights = MCEM_res[['weights']]
  cluster_mat = MCEM_res[['cluster_mat']]
  response_pcs_array = MCEM_res[['pcs_array_response']]
  predictor_pcs_array = MCEM_res[['pcs_array_predictors']]
  n_profiles = data_inputs[['n_profiles']]
  n_profiles_TS = data_inputs[['n_profiles_TS']]
  G = params[['G']]
  cluster_mat_BGC <- cluster_mat[data_inputs[['BGC']], ,drop = F]
  predictor_pcs_array_BGC <- predictor_pcs_array[data_inputs[['BGC']], , , drop = F]
  
  cmi = cluster_mat - 1
  cmi_BGC = cluster_mat_BGC - 1
  mode(cmi) = 'integer'
  mode(cmi_BGC) = 'integer'
  Linv_mats <- MCEM_res[['Linv_mat']]
  vecchia_info <- MCEM_res[['vecchia_info']]

  # Update all the parameters
  
  u_mats = compute_u_matrices(Phi_x_Phi_resp = data_inputs[['Phi_x_Phi_response']], 
                              weights = weights,
                              cmi = cmi,
                              cmi_BGC = cmi_BGC,
                              n_e_step_samples=ncol(cmi),
                              G = G,
                              pcs_r = response_pcs_array,
                              Phi_x_Phi_pred = data_inputs[['Phi_x_Phi_predictors']],
                              pcs_p1 = predictor_pcs_array,
                              pcs_p2 = predictor_pcs_array_BGC)

  if (length(params[['lambda_mean_response']]) == 1 | (iter-1) %% params[['cv_skip']] != 0){
      parameters[['means_response']] = update_means(U = u_mats[['mean_resp']],
                                                    Phi_prof = data_inputs[['Phi_response_prof']],
                                                    profs = data_inputs[['response_profs']],
                                                    cmi = cmi_BGC,
                                                    pcs_array1 = predictor_pcs_array_BGC,
                                                    Omegas1 = parameters[['Omega_response1']],
                                                    weights = weights,
                                                    pen_mat = data_inputs[['penalty_mat_response']],
                                                    G = G,
                                                    Omegas2 = parameters[['Omega_response2']],
                                                    pcs_array2 = response_pcs_array, 
                                                    lambda = parameters[['lambda_mean_response']])
  } else {
    AIC = Inf
    for (lambda in params[['lambda_mean_response']]){
      old_means = parameters[['means_response']]; old_lambda = parameters[['lambda_mean_response']]
      parameters[['means_response']] = update_means(U = u_mats[['mean_resp']],
                                                    Phi_prof = data_inputs[['Phi_response_prof']],
                                                    profs = data_inputs[['response_profs']],
                                                    cmi = cmi_BGC,
                                                    pcs_array1 = predictor_pcs_array_BGC,
                                                    Omegas1 = parameters[['Omega_response1']],
                                                    weights = weights,
                                                    pen_mat = data_inputs[['penalty_mat_response']],
                                                    G = G,
                                                    Omegas2 = parameters[['Omega_response2']],
                                                    pcs_array2 = response_pcs_array, 
                                                    lambda = lambda)
      parameters[['lambda_mean_response']] = lambda
      model = compute_information_criteria(model, params, data_inputs, parameters, u_mats)
      if(model$AIC > AIC){
        parameters[['lambda_mean_response']] = old_lambda; parameters[['means_response']] = old_means
      }
      AIC = ifelse(model$AIC > AIC, AIC, model$AIC)
    }
  }
  
  if (length(params[['lambda_mean_predictors']]) == 1 | (iter-1) %% params[['cv_skip']] != 0){
    parameters[['means_predictors']] = update_means(U = u_mats[['mean_pred']],
                                                    Phi_prof = data_inputs[['Phi_predictors_prof']],
                                                    profs = data_inputs[['predictor_profs']],
                                                    cmi = cmi,
                                                    pcs_array1 = predictor_pcs_array,
                                                    Omegas1 = parameters[['Omega_predictors']],
                                                    weights = weights,
                                                    pen_mat = data_inputs[['penalty_mat_predictors']],
                                                    G = G,
                                                    lambda = parameters[['lambda_mean_predictors']])
  } else {
    AIC = Inf
    for (lambda in params[['lambda_mean_predictors']]){
      old_means = parameters[['means_predictors']]; old_lambda = parameters[['lambda_mean_predictors']]
      parameters[['means_predictors']] = update_means(U = u_mats[['mean_pred']],
                                                      Phi_prof = data_inputs[['Phi_predictors_prof']],
                                                      profs = data_inputs[['predictor_profs']],
                                                      cmi = cmi,
                                                      pcs_array1 = predictor_pcs_array,
                                                      Omegas1 = parameters[['Omega_predictors']],
                                                      weights = weights,
                                                      pen_mat = data_inputs[['penalty_mat_predictors']],
                                                      G = G,
                                                      lambda = lambda)
      parameters[['lambda_mean_predictors']] = lambda
      model = compute_information_criteria(model, params, data_inputs, parameters, u_mats)
      if(model$AIC > AIC){
        parameters[['lambda_mean_predictors']] = old_lambda; parameters[['means_predictors']] = old_means
      }
      AIC = ifelse(model$AIC > AIC, AIC, model$AIC)
    }
  }

  if (length(params[['lambda_pcs_response']]) == 1 | (iter-1) %% params[['cv_skip']] != 0){
    parameters[['Omega_response2']] = update_principal_components(u_mats = u_mats[['pc_resp']],
                                                                  cmi = cmi_BGC,
                                                                  pcs_r1 = response_pcs_array,
                                                                  profs = data_inputs[['response_profs']],
                                                                  Phi_prof = data_inputs[['Phi_response_prof']],
                                                                  Omegas1 = parameters[['Omega_response2']],
                                                                  means = parameters[['means_response']], 
                                                                  weights = weights,
                                                                  pen_mat = data_inputs[['penalty_mat_response']],
                                                                  G = G,
                                                                  pcs_r2 = predictor_pcs_array_BGC,
                                                                  Omegas2 = parameters[['Omega_response1']],
                                                                  lambda = parameters[['lambda_pcs_response']])
  } else {
    AIC = Inf
    for (lambda in params[['lambda_pcs_response']]){
      old_pcs = parameters[['Omega_response2']]; old_lambda = parameters[['lambda_pcs_response']]
      parameters[['Omega_response2']] = update_principal_components(u_mats = u_mats[['pc_resp']],
                                                                    cmi = cmi_BGC,
                                                                    pcs_r1 = response_pcs_array,
                                                                    profs = data_inputs[['response_profs']],
                                                                    Phi_prof = data_inputs[['Phi_response_prof']],
                                                                    Omegas1 = parameters[['Omega_response2']],
                                                                    means = parameters[['means_response']], 
                                                                    weights = weights,
                                                                    pen_mat = data_inputs[['penalty_mat_response']],
                                                                    G = G,
                                                                    pcs_r2 = predictor_pcs_array_BGC,
                                                                    Omegas2 = parameters[['Omega_response1']],
                                                                    lambda = lambda)
      parameters[['lambda_pcs_response']] = lambda
      model = compute_information_criteria(model, params, data_inputs, parameters, u_mats)
      if(model$AIC > AIC){
        parameters[['lambda_pcs_response']] = old_lambda; parameters[['Omega_response2']] = old_pcs
      }
      AIC = ifelse(model$AIC > AIC, AIC, model$AIC)
    }
  }
  
  if (length(params[['lambda_lt']]) == 1 | (iter-1) %% params[['cv_skip']] != 0){
    parameters[['Omega_response1']] = update_principal_components(u_mats = u_mats[['lt']],
                                                                  cmi = cmi_BGC,
                                                                  pcs_r1 = predictor_pcs_array_BGC,
                                                                  profs = data_inputs[['response_profs']],
                                                                  Phi_prof = data_inputs[['Phi_response_prof']],
                                                                  Omegas1 = parameters[['Omega_response1']],
                                                                  means = parameters[['means_response']], 
                                                                  weights = weights,
                                                                  pen_mat = data_inputs[['penalty_mat_response']],
                                                                  G = G,
                                                                  pcs_r2 = response_pcs_array,
                                                                  Omegas2 = parameters[['Omega_response2']],
                                                                  lambda = parameters[['lambda_lt']])
  } else {
    AIC = Inf
    for (lambda in params[['lambda_lt']]){
      old_pcs = parameters[['Omega_response1']]; old_lambda = parameters[['lambda_lt']]
      parameters[['Omega_response1']] = update_principal_components(u_mats = u_mats[['lt']],
                                                                    cmi = cmi_BGC,
                                                                    pcs_r1 = predictor_pcs_array_BGC,
                                                                    profs = data_inputs[['response_profs']],
                                                                    Phi_prof = data_inputs[['Phi_response_prof']],
                                                                    Omegas1 = parameters[['Omega_response1']],
                                                                    means = parameters[['means_response']], 
                                                                    weights = weights,
                                                                    pen_mat = data_inputs[['penalty_mat_response']],
                                                                    G = as.integer(G),
                                                                    pcs_r2 = response_pcs_array,
                                                                    Omegas2 = parameters[['Omega_response2']],
                                                                    lambda = lambda)
      parameters[['lambda_lt']] = lambda
      model = compute_information_criteria(model, params, data_inputs, parameters, u_mats)
      if(model$AIC > AIC){
        parameters[['lambda_lt']] = old_lambda; parameters[['Omega_response1']] = old_pcs
      }
      AIC = ifelse(model$AIC > AIC, AIC, model$AIC)
    }
  }
  
  if (length(params[['lambda_pcs_predictors']]) == 1 | (iter-1) %% params[['cv_skip']] != 0){
    parameters[['Omega_predictors']] = update_principal_components(u_mats = u_mats[['pc_pred']],
                                                                   cmi = cmi,
                                                                   pcs_r1 = predictor_pcs_array,
                                                                   profs = data_inputs[['predictor_profs']],
                                                                   Phi_prof = data_inputs[['Phi_predictors_prof']],
                                                                   Omegas1 = parameters[['Omega_predictors']],
                                                                   means = parameters[['means_predictors']], 
                                                                   weights = weights,
                                                                   pen_mat = data_inputs[['penalty_mat_predictors']],
                                                                   G = G,
                                                                   lambda = parameters[['lambda_pcs_predictors']])
  } else {
    AIC = Inf
    for (lambda in params[['lambda_pcs_predictors']]){
      old_pcs = parameters[['Omega_predictors']]; old_lambda = parameters[['lambda_pcs_predictors']]
      parameters[['Omega_predictors']] = update_principal_components(u_mats = u_mats[['pc_pred']],
                                                                     cmi = cmi,
                                                                     pcs_r1 = predictor_pcs_array,
                                                                     profs = data_inputs[['predictor_profs']],
                                                                     Phi_prof = data_inputs[['Phi_predictors_prof']],
                                                                     Omegas1 = parameters[['Omega_predictors']],
                                                                     means = parameters[['means_predictors']], 
                                                                     weights = weights,
                                                                     pen_mat = data_inputs[['penalty_mat_predictors']],
                                                                     G = G,
                                                                     lambda = lambda)
      parameters[['lambda_pcs_predictors']] = lambda
      model = compute_information_criteria(model, params, data_inputs, parameters, u_mats)
      if(model$AIC > AIC){
        parameters[['lambda_pcs_predictors']] = old_lambda; parameters[['Omega_predictors']] = old_pcs
      }
      AIC = ifelse(model$AIC > AIC, AIC, model$AIC)
    }
  }
  
  parameters[['measurement_error_response']] = update_measurement_error(cluster_mat = cluster_mat_BGC,
                                                                        pcs_array1 = predictor_pcs_array_BGC,
                                                                        profs = data_inputs[['response_profs']],
                                                                        phi_prof = data_inputs[['Phi_response_prof']],
                                                                        means = parameters[['means_response']],
                                                                        Omegas1 = parameters[['Omega_response1']],
                                                                        n_profiles = data_inputs[['n_profiles']],
                                                                        Omegas2 = parameters[['Omega_response2']],
                                                                        pcs_array2 = response_pcs_array,
                                                                        G = G)
  
  parameters[['measurement_error_predictors']] = update_measurement_error(cluster_mat = cluster_mat,
                                                                          pcs_array1 = predictor_pcs_array,
                                                                          profs = data_inputs[['predictor_profs']],
                                                                          phi_prof = data_inputs[['Phi_predictors_prof']],
                                                                          means = parameters[['means_predictors']],
                                                                          Omegas1 = parameters[['Omega_predictors']],
                                                                          n_profiles = data_inputs[['n_profiles_TS']],
                                                                          prof_ns = data_inputs[['profile_lengths_predictors']],
                                                                          n_vars = length(params[['X']]),
                                                                          G = G)

  parameters[['variances_response']] <- update_variances(cluster_mat = cluster_mat_BGC,
                                                         pcs_array = response_pcs_array,
                                                         range_params = parameters[['range_params_r']],
                                                         G = params[['G']],
                                                         n_pcs = dim(response_pcs_array)[3],
                                                         MC = ncol(cluster_mat), 
                                                         cov_fun = params[['covariance_function']],
                                                         independent_model = is.null(params[['covariance_function']]),
                                                         Linv = Linv_mats[['response']],
                                                         variances = parameters[['variances_response']], 
                                                         vecchia_info =  vecchia_info[['response']])
  
  parameters[['variances_predictors']] <- update_variances(cluster_mat = cluster_mat,
                                                           pcs_array = predictor_pcs_array,
                                                           range_params = parameters[['range_params_p']],
                                                           G = params$G,
                                                           n_pcs = dim(predictor_pcs_array)[3],
                                                           MC = ncol(cluster_mat),
                                                           cov_fun = params[['covariance_function']],
                                                           independent_model = is.null(params[['covariance_function']]),
                                                           Linv = Linv_mats[['predictors']],
                                                           variances = parameters[['variances_predictors']],
                                                           vecchia_info =  vecchia_info[['predictors']])
  
  # Sort
  for(g in 1:params$G){
    inds = order(parameters$variances_predictors[[g]], decreasing=T)
    parameters$Omega_predictors[[g]] = parameters$Omega_predictors[[g]][,inds]
    parameters$variances_predictors[[g]] = parameters$variances_predictors[[g]][inds]
    parameters$Omega_response1[[g]] = parameters$Omega_response1[[g]][,inds]
    parameters$range_params_p[[g]] = parameters$range_params_p[[g]][,inds]
  }
  for(g in 1:params$G){
    inds = order(parameters$variances_response[[g]], decreasing=T)
    parameters$Omega_response2[[g]] = parameters$Omega_response2[[g]][,inds]
    parameters$variances_response[[g]] = parameters$variances_response[[g]][inds]
    parameters$range_params_r[[g]] = parameters$range_params_r[[g]][,inds]
  }
  
  # Orthogonality constraints
  ortho_info <- orthogonalize(Omegas1 = parameters[['Omega_predictors']],
                              variances1 = parameters[['variances_predictors']],
                              inner_prod_mat1 = data_inputs[['inner_prod_mat_predictors']],
                              G = G,
                              Omegas2 = parameters[['Omega_response2']], 
                              variances2 = parameters[['variances_response']],
                              inner_prod_mat2 = data_inputs[['inner_prod_mat_response']],
                              Lambdas = parameters[['Omega_response1']])
  
  parameters[['variances_predictors']] <- lapply(1:G, function(g) {ortho_info[[g]][[1]]})
  parameters[['Omega_response1']] <- lapply(1:G, function(g) {ortho_info[[g]][[2]]})
  parameters[['Omega_predictors']] <- lapply(1:G, function(g) {ortho_info[[g]][[3]]})
  parameters[['variances_response']] <- lapply(1:G, function(g) {ortho_info[[g]][[4]]})
  parameters[['Omega_response2']] <- lapply(1:G, function(g) {ortho_info[[g]][[5]]})
  
  if (!is.null(params[['covariance_function']])){
    parameters[['range_params_r']] <- update_range_params(G = params[['G']],
                                                          n_pcs = params[['pc_response2']], 
                                                          variances = parameters[['variances_response']],
                                                          cluster_mat = cluster_mat_BGC, 
                                                          pcs_array = response_pcs_array,
                                                          vecchia_info = vecchia_info[['response']],
                                                          mc_weights = weights, covfun_name = params[['covariance_function']],
                                                          start_params_previous = parameters[['range_params_r']],
                                                          maxit = params[['maxit']])
  } 
  if (!is.null(params[['covariance_function']])) {
    parameters[['range_params_p']] <- update_range_params(G = params[['G']],
                                                          n_pcs = params$pc_predictors, 
                                                          variances = parameters[['variances_predictors']],
                                                          cluster_mat = cluster_mat, 
                                                          pcs_array = predictor_pcs_array,
                                                          vecchia_info = vecchia_info[['predictors']],
                                                          mc_weights = weights, covfun_name = params[['covariance_function']],
                                                          maxit = params[['maxit']],
                                                          start_params_previous = parameters[['range_params_p']])
  }
  
  cluster_list = update_cluster_membership(cluster_mat, n_profiles_TS, G)
  parameters[['cluster_membership']] = cluster_list[[1]]
  parameters[['conditional_cluster_probs']] = cluster_list[[2]]
  parameters[['theta']] = update_theta(c(0,3), nn_list = data_inputs[['nn_list']], 
                                       cluster_mat = cluster_mat,
                                       n_profiles = n_profiles_TS,G = G,
                                       weights = weights,
                                       theta = parameters[['theta']])
  
  parameters[['marginal_probabilities']] <- compute_marginal_probabilities(parameters[['cluster_membership']], parameters[['theta']], neighbors = data_inputs[['nn_list']], n_profiles, G)
  
  return(parameters)
}