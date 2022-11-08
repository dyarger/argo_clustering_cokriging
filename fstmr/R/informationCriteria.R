compute_information_criteria <- function(model, params, data_inputs, parameters, u_mats=NULL, MCEM_res=NULL){
  
  lik_res = eval_likelihood_est(model, params, data_inputs, parameters, MCEM_res)
  likelihood = lik_res[['likelihood']]
  MCEM_res = lik_res[['MCEM_res']]
  
  if (is.null(u_mats)){
    weights = MCEM_res[['weights']]
    cluster_mat = MCEM_res[['cluster_mat']]
    response_pcs_array = MCEM_res[['pcs_array_response']]
    G = params[['G']]
    cmi = cluster_mat - 1
    mode(cmi) = 'integer'
    
    if (class(model) == 'one_var'){
      
      u_mats = compute_u_matrices(Phi_x_Phi_resp = data_inputs[['Phi_x_Phi_response']], 
                                  weights = weights,
                                  cmi = NULL,
                                  cmi_BGC = cmi,
                                  n_e_step_samples=ncol(cmi),
                                  G = G,
                                  pcs_r = response_pcs_array)
    } else if (class(model) == 'mult_var_ind'){
      predictor_pcs_array = MCEM_res[[3]]
      cluster_mat_BGC <- cluster_mat[data_inputs[['BGC']], ,drop = F]
      predictor_pcs_array_BGC <- predictor_pcs_array[data_inputs[['BGC']], , , drop = F]
      
      cmi_BGC = cluster_mat_BGC -1 
      mode(cmi_BGC) = 'integer'
      
      weights = rep(1, ncol(cmi))
      
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
    } else if (class(model) == 'mult_var_spat'){
      predictor_pcs_array = MCEM_res[['pcs_array_predictors']]
      cluster_mat_BGC <- cluster_mat[data_inputs[['BGC']], ,drop = F]
      predictor_pcs_array_BGC <- predictor_pcs_array[data_inputs[['BGC']], , , drop = F]
      
      cmi_BGC = cluster_mat_BGC - 1
      mode(cmi_BGC) = 'integer'
      
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
    }
  }
  
  n_params = compute_n_params(params=params, u_mats=u_mats,
                              lambda_mean_resp=parameters[['lambda_mean_response']],
                              lambda_pcs_resp=parameters[['lambda_pcs_response']],
                              pen_mat_resp=data_inputs[['penalty_mat_response']],
                              lambda_mean_pred=parameters[['lambda_mean_predictors']],
                              lambda_pcs_pred=parameters[['lambda_pcs_predictors']],
                              pen_mat_pred=data_inputs[['penalty_mat_predictors']],
                              lambda_lt=parameters[['lambda_lt']])

  model[['BIC']] = BIC_fun(prof_lengths_resp = data_inputs[['profile_lengths_response']],
                           n_params = n_params,
                           likelihood = likelihood,
                           prof_lengths_pred=data_inputs[['profile_lengths_predictors']])

  model[['AIC']] = AIC_fun(n_params = n_params, likelihood = likelihood) 
  model[['likelihood']] = likelihood
  model[['likelihood_spat']] <- MCEM_res[['likelihood_spat']] 
  model[['likelihood_ind']] <- MCEM_res[['likelihood_ind']]
  model[['likelihood_MRF']] <- lik_res[['like_MRF']]
  return(model)
}


eval_likelihood_est <- function(model, params, data_inputs, parameters, MCEM_res=NULL){
  UseMethod('eval_likelihood_est')
}

eval_likelihood_est.one_var <- function(model, params, data_inputs, parameters, MCEM_res=NULL) {
  use_MRF_old <- params[['use_MRF']]
  if (is.null(params[['m_AIC_response']])) {
    params[['m_AIC_response']] <- params[['m_train_response']]
  } else if (is.na(params[['m_AIC_response']])) {
    params[['m_AIC_response']] <- params[['m_train_response']]
  }
  if (is.null(params[['MC_max_AIC']])) {
    params[['MC_max_AIC']] <- params[['MC_max']]
  } else if (is.na(params[['m_AIC_response']])) {
    params[['MC_max_AIC']] <- params[['MC_max']]
  }
  if (params[['use_MRF']] | is.null(MCEM_res)) {
    params[['use_MRF']] <- F
    params[['m_train_response']] <- params[['m_AIC_response']]
    params[['MC_max']] <- params[['MC_max_AIC']]
    MCEM_res = Estep(model, params, data_inputs, parameters)
  }
  
  G = params[['G']]
  cluster_mat = MCEM_res[['cluster_mat']]
  cluster_mat_i = cluster_mat - 1
  mode(cluster_mat_i) = 'integer'
  weights = MCEM_res[['weights']]
  
  MC = ncol(cluster_mat)
  
  means_mat_resp = as.matrix(do.call(cbind, parameters[['means_response']]))
  sigma2 = parameters[['measurement_error_response']]
  
  s3 <- rep(0, MC)
  for(i in 1:data_inputs[['n_profiles']]){
    
    neighbors <- data_inputs[['nn_list']][[i]]
    
    cluster_neighbors = cluster_mat[neighbors,, drop = F]
    
    for (mc in 1:MC) {
      if (use_MRF_old) {
        s3[mc] = s3[mc] + parameters[['theta']] * sum(
          cluster_neighbors[,mc] == cluster_mat[i, mc] & neighbors < i)
      } 
    }
  }
  likelihood_vec <- s3 + MCEM_res[['likelihood_spat']] - MCEM_res[['likelihood_ind']]
  c_use <- max(likelihood_vec)
  likelihood_return <-   c_use + log(sum(exp(likelihood_vec - c_use)))
  
  return(list(likelihood = likelihood_return, MCEM_res = MCEM_res, like_MRF = s3))
}

eval_likelihood_est.mult_var_ind <- function(model, params, data_inputs, parameters, MCEM_res=NULL) {
  if (is.null(MCEM_res)){
    MCEM_res = Estep(model, params, data_inputs, parameters)
  }
  
  G = params[['G']]
  s1_resp = 0
  s1_pred = 0
  s2 = 0
  s3 = 0
  BGC_index = 1
  cluster_mat = MCEM_res[['cluster_mat']]
  cluster_mat_i = cluster_mat - 1
  mode(cluster_mat_i) = 'integer'
  weights = rep(1, ncol(cluster_mat))
  
  MC = ncol(cluster_mat)
  
  means_mat_pred = as.matrix(do.call(cbind, parameters[['means_predictors']]))
  means_mat_resp = as.matrix(do.call(cbind, parameters[['means_response']]))
  sigma2_resp <- parameters[['measurement_error_response']]
  
  scores_covariance_chols <- lapply(1:params[['G']], function(g) {
    C_mat <- base::diag(parameters[['variances_predictors']][[g]]) %*%
      base::t(parameters[['Lambdas']][[g]])
    tryCatch({
      base::t(base::chol(base::rbind(base::cbind(base::diag(parameters[['variances_response']][[g]]), 
                                                 base::t(C_mat)), 
                                     base::cbind(C_mat, base::diag(parameters[['variances_predictors']][[g]])))))
    }, error = function(z) {
      print(base::rbind(base::cbind(base::diag(parameters[['variances_response']][[g]] + 1), 
                                    base::t(C_mat)), 
                        base::cbind(C_mat, base::diag(parameters[['variances_predictors']][[g]] + 1))))
      base::t(base::chol(base::rbind(base::cbind(base::diag(parameters[['variances_response']][[g]] + 1), 
                                                 base::t(C_mat)), 
                                     base::cbind(C_mat, base::diag(parameters[['variances_predictors']][[g]] + 1)))))
    })
  })
  
  
  for(i in 1:data_inputs[['n_profiles_TS']]){
    is_BGC <- data_inputs[['BGC']][i]
    if (is_BGC) {
      i_resp = BGC_index
      BGC_index = BGC_index + 1
      mi_resp = length(data_inputs[['response_profs']][[i_resp]])
    }
    
    sigma2_pred <- rep(parameters[['measurement_error_predictors']],
                       each = data_inputs[['profile_lengths_predictors']][i,1])
    
    cluster_neighbors = cluster_mat[data_inputs[['nn_list']][[i]],, drop = F]
    
    
    diff = c_compute_squared_sparse(cluster_mat_i = cluster_mat_i[i,],
                                    pcs_mat = as.matrix(MCEM_res[['pcs_array_predictors']][i,,]),
                                    profile = data_inputs[['predictor_profs']][[i]],
                                    G = as.integer(G),
                                    basis_eval = data_inputs[['Phi_predictors_prof']][[i]],
                                    means_mat = means_mat_pred,
                                    Omegas = parameters[['Omega_predictors']],
                                    weights = weights)
    
    s1_pred = s1_pred - 0.5 * MC* sum(log(sigma2_pred)) -
      .5 *  (sum(diff/sigma2_pred)) -
      .5 * length(diff) * MC * log(2 * pi)
    
    if (is_BGC) {
      diff <- c_compute_squared_sparse(cluster_mat_i = cluster_mat_i[i,],
                                       pcs_mat = as.matrix(MCEM_res[['pcs_array_response']][i_resp,,]),
                                       profile = data_inputs[['response_profs']][[i_resp]],
                                       G = as.integer(G),
                                       basis_eval = data_inputs[['Phi_response_prof']][[i_resp]],
                                       means_mat = means_mat_resp,
                                       Omegas = parameters[['Omega_response']],
                                       weights = weights)
      
      s1_resp = s1_resp - 0.5 * MC* mi_resp * sum(log(sigma2_resp)) -
        .5 *  ( sum(diff/sigma2_resp)) -
        .5 * mi_resp * MC * log(2 * pi)
    }
    
    if (params[['use_MRF']]){
      U_ik_clust <- sapply(1:G, function(g) {
        parameters[['theta']] * colSums(cluster_neighbors == g) * weights
      })
      s3 = s3 + sum(U_ik_clust[cbind(1:MC, cluster_mat[i,])] - log(rowSums(exp(U_ik_clust)))) 
    } else {
      s3 = 0
    }
    
    for (mc in 1:MC) {
      c_i <- cluster_mat_i[i, mc] + 1
      if (is_BGC) {
        scores <- c(MCEM_res[['pcs_array_response']][i_resp,mc,],
                    MCEM_res[['pcs_array_predictors']][i,mc,])

        s2 = s2 +  -1/2 *  (2*sum(log(base::diag(scores_covariance_chols[[c_i]])))+ 
                              sum(base::solve(scores_covariance_chols[[c_i]], scores)^2) + 
                              length(scores) * log(2* pi))
        s2 = s2 +  - sum(log(base::diag(scores_covariance_chols[[c_i]])))-
           sum(base::solve(scores_covariance_chols[[c_i]], scores)^2)/2 -
                              length(scores) * log(2* pi)/2
      } else {
        scores <- MCEM_res[['pcs_array_predictors']][i,mc,]
        cov_mat_scores <- diag(parameters[['variances_predictors']][[c_i]])
        s2 = s2 - sum(log(parameters[['variances_predictors']][[c_i]]))/2-
          sum(scores^2 /parameters[['variances_predictors']][[c_i]])/2 - 
          length(scores) * log(2*pi)/2
      }
    }
  }
  
  s1 = (s1_resp + s1_pred)/MC
  s2 = s2/MC
  s3 = s3/MC
  
  return(list(likelihood = s1 + s2 + s3, MCEM_res = MCEM_res))
}

eval_likelihood_est.mult_var_spat <- function(model, params, data_inputs, parameters, MCEM_res=NULL) {
  
  use_MRF_old <- params[['use_MRF']]
  if (is.null(params[['m_AIC_response']])) {
    params[['m_AIC_response']] <- params[['m_train_response']]
  } else if (is.na(params[['m_AIC_response']])) {
    params[['m_AIC_response']] <- params[['m_train_response']]
  }
  if (is.null(params[['m_AIC_predictors']])) {
    params[['m_AIC_predictors']] <- params[['m_train_predictors']]
  } else if (is.na(params[['m_AIC_predictors']])) {
    params[['m_AIC_predictors']] <- params[['m_train_predictors']]
  }
  if (is.null(params[['MC_max_AIC']])) {
    params[['MC_max_AIC']] <- params[['MC_max']]
  } else if (is.na(params[['MC_max_AIC']])) {
    params[['MC_max_AIC']] <- params[['MC_max']]
  }
  
  params[['use_MRF']] <- F
  params[['m_train_response']] <- params[['m_AIC_response']]
  params[['m_train_predictors']] <- params[['m_AIC_predictors']]
  params[['MC_max']] <- params[['MC_max_AIC']]
  if (is.null(params[['use_MRF']])) {
    MCEM_res = Estep(model, params, data_inputs, parameters)
  } else if (params[['use_MRF']] | is.null(MCEM_res)) {
    MCEM_res = Estep(model, params, data_inputs, parameters)
  }
  
  G = params[['G']]
  cluster_mat = MCEM_res[['cluster_mat']]
  cluster_mat_i = cluster_mat - 1
  mode(cluster_mat_i) = 'integer'
  weights = MCEM_res[['weights']]
  
  MC = ncol(cluster_mat)
  
  means_mat_pred = as.matrix(do.call(cbind, parameters[['means_predictors']]))
  means_mat_resp = as.matrix(do.call(cbind, parameters[['means_response']]))
  sigma2_resp <- parameters[['measurement_error_response']]
  
  
  # evaluate P(Z)
  s3 <- rep(0, MC)
  for(i in 1:data_inputs[['n_profiles_TS']]){
    
    neighbors <- data_inputs[['nn_list']][[i]]
    
    cluster_neighbors = cluster_mat[neighbors,, drop = F]
    
    for (mc in 1:MC) {
      if (use_MRF_old) {
        s3[mc] = s3[mc] + parameters[['theta']] * sum(
          cluster_neighbors[,mc] == cluster_mat[i, mc] & neighbors < i)
      } 
    }
  }
  likelihood_vec <- s3 + MCEM_res[['likelihood_spat']] - MCEM_res[['likelihood_ind']]
  c_use <- max(likelihood_vec)
  likelihood_return <-   c_use + log(sum(exp(likelihood_vec - c_use)))

  return(list(likelihood = likelihood_return, MCEM_res = MCEM_res, like_MRF = s3))
}


BIC_fun <- function(prof_lengths_resp, n_params, likelihood, prof_lengths_pred=NULL) {
  n_tilde = sum(prof_lengths_resp) + ifelse(is.null(prof_lengths_pred), 0,
                                     sum(prof_lengths_pred))
  -2* likelihood +  n_params * log(n_tilde)
}

AIC_fun <- function(n_params, likelihood) {
  -2* likelihood + 2 * n_params
}

compute_n_params <- function(params, u_mats, cluster_mat, lambda_mean_resp, lambda_pcs_resp,
                             pen_mat_resp, lambda_mean_pred=NULL, lambda_pcs_pred=NULL,
                             pen_mat_pred=NULL, lambda_lt=NULL){
  
  G = params[['G']]
  
  n_basis_resp = params[['basis_response']][['nbasis']]
  n_basis_pred = ifelse(is.null(params[['basis_response']]), 0,
                        Reduce(`+`,params[['basis_response']][['nbasis']]))
  
  Q_resp1 = ifelse(is.null(params[['pc_response1']]), 0, params[['pc_response1']])
  Q_resp2 = ifelse(is.null(params[['pc_response2']]), 0, params[['pc_response2']])
  Q_resp = ifelse(is.null(params[['pc_response']]), 0, params[['pc_response']])
  Q_pred = ifelse(is.null(params[['pc_predictors']]), 0, params[['pc_predictors']])
  
  n_range = ifelse(is.null(params[['covariance_function']]), 0, 
                   ifelse(params[['covariance_function']] %in%
                            c('exponential_spheretime_warp',
                              'matern_spheretime_warp'),
                          7, 2))
  
  pen_mat_resp = Matrix::bdiag(lapply(1:params[['G']], function(g) pen_mat_resp))
  if (!is.null(pen_mat_pred)) {
    pen_mat_pred = Matrix::bdiag(lapply(1:params[['G']], function(g) pen_mat_pred))
  }
  
  effective_df_mean_resp = tryCatch({sum(Matrix::rowSums(sparseinv::Takahashi_Davis(u_mats[['mean_resp']] + lambda_mean_resp * pen_mat_resp) * u_mats[['mean_resp']]))},
                                    error = function(x) {return(0)})
  
  effective_df_pcs_resp = sum(sapply(1:(Q_resp + Q_resp2), function(q){
    tryCatch({sum(Matrix::rowSums(sparseinv::Takahashi_Davis(u_mats[['pc_resp']][[q]] + lambda_pcs_resp*pen_mat_resp) * u_mats[['pc_resp']][[q]]))},
             error = function(x) {return(0)})
  }))
  
  effective_df_mean_pred = ifelse(
    is.null(lambda_mean_pred), 0,
    tryCatch({sum(Matrix::rowSums(sparseinv::Takahashi_Davis(u_mats[['mean_pred']] + lambda_mean_pred*pen_mat_pred) * u_mats[['mean_pred']]))},
             error = function(x) {return(0)})
  )
  
  effective_df_pcs_pred = ifelse(
    is.null(lambda_pcs_pred), 0,
    sum(sapply(1:Q_pred, function(q){
      tryCatch({sum(Matrix::rowSums(sparseinv::Takahashi_Davis(u_mats[['pc_pred']][[q]] + lambda_pcs_pred*pen_mat_pred) * u_mats[['pc_pred']][[q]]))},
               error = function(x) {return(0)})
    }))
  )
  
  effective_df_lt = ifelse(
    is.null(lambda_lt), 0,
    sum(sapply(1:Q_pred, function(q){
      tryCatch({sum(Matrix::rowSums(sparseinv::Takahashi_Davis(u_mats[['lt']][[q]] + lambda_lt*pen_mat_resp) * u_mats[['lt']][[q]]))},
               error = function(x) {return(0)})
    }))
  )
  
  n_pred = length(params[['preds']])
  n_resp = 1
  
  effective_df_mean_pred + # mean functions
    effective_df_mean_resp + # mean functions
    effective_df_pcs_pred + # Omega
    effective_df_lt + # lambda1
    effective_df_pcs_resp  + # one_var
    G * (Q_resp2 + Q_pred + Q_resp) + # variance of scores
    ifelse(!is.null(params[['covariance_function']]), G * Q_resp2 * n_range, 0) +  # range parameters for response
    ifelse(!is.null(params[['covariance_function']]), G * Q_pred * n_range, 0) + # range parameters for predictors
    ifelse(!is.null(params[['covariance_function']]), G * Q_resp * n_range, 0) + # range parameters for one_var
    ifelse(params[['use_MRF']], 1, 0) +      # MRF parameter
    n_pred + n_resp        # measurement error variances
}

