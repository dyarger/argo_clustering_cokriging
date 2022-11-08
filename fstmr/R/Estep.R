Estep <- function(model, params, data_inputs, parameters){
  UseMethod('Estep')
}

Estep.one_var <- function(model, params, data_inputs, parameters){
  
  MC_max = params[['MC_max']]
  G = params[['G']]
  
  if (params[['cluster_sampling']] == 'gibbs') {
    MC_max <- MC_max + 10
  }
  
  # Frequently used constants
  n_profiles = data_inputs[['n_profiles']]
  
  cluster_mat = matrix(0, ncol=MC_max, nrow=n_profiles)
  # Matrices of principal component samples
  pcs_array_response = array(NA, c(n_profiles, MC_max, params[['pc_response']]))
  
  like_spatial <- rep(0, MC_max)
  like_ind <- rep(0, MC_max)

  probs <- return_probabilities(vars_r = parameters$variances_response,
                                response_profs = data_inputs$response_profs,
                                phi_resp = data_inputs$Phi_response_prof,
                                means_resp = parameters$means_response,
                                me_resp = parameters$measurement_error_response,
                                Omegas_r1 = parameters$Omega_response,
                                G = G,
                                n_profiles = n_profiles)
  
  p_mat <- probs[[1]]; like_mat <- probs[[2]]
  
  if(params[['cluster_sampling']] != 'gibbs' & params[['use_MRF']]) {
    cluster_mat = rPotts(G=G,
                         neighbors = data_inputs[['nn_list']],
                         p_mat=p_mat,
                         likelihoods = like_mat,
                         theta = parameters[['theta']],
                         n_samples = MC_max)
  } 
  if (params[['cluster_sampling']] != 'gibbs' & params[['use_MRF']] == F) {
    cluster_mat <- matrix(nrow = n_profiles, ncol = MC_max)
    for (mc in 1:MC_max) {
      for (x in 1:n_profiles) {
        cluster_mat[x,mc] <- sample(1:G, size = 1, prob = p_mat[x,])
      }
    }
  }
  
  # save vecchia information so that we don't have to redo it in update_variances/range_params
  Linv_mat <- list()
  vecchia_info <- list()
  E_step_list = list()
  
  for (mc in 1:MC_max){
    
    if (params[['cluster_sampling']] == 'gibbs') {
      if (mc == 1 & !is.null(parameters[['cluster_membership']])) {
        cluster_mat[,mc] <- parameters[['cluster_membership']]
      } else if (mc == 1 & is.null(parameters[['cluster_membership']])) {
        cluster_mat[,mc] <- parameters[['cluster_mat']]
      } else {
        probs_gibbs <- compute_probabilities_gibbs(pcs_array_response[,mc-1,],
                                                   parameters[['marginal_probabilities']],
                                                   n_profiles = n_profiles, n_scores = dim(pcs_array_response)[3],
                                                   G = params[['G']], parameters, data_inputs)
        for (x in 1:data_inputs$n_profiles){
          cluster_mat[x, mc] = sample(1:params[['G']], size = 1, prob = probs_gibbs[x,])
        }
      }
    }
    
    cmi = cluster_mat[,mc] -1
    mode(cmi) ='integer'
    
    if (mc == 1 || !all(cluster_mat[,mc] == cluster_mat[,mc - 1])){
      E_step_list[['UTU']] =  c_compute_UTU_single(data_inputs[['Phi_x_Phi_response']],
                                                   parameters[['Omega_response']],
                                                   as.numeric(parameters[['measurement_error_response']]),
                                                   cmi, 
                                                   data_inputs[['n_profiles']], 
                                                   params[['G']])
      
      
      means_mat_r = as.matrix(do.call(cbind, parameters[['means_response']]))
      UTX = c_compute_UTX_single(data_inputs[['Phi_response_prof']], data_inputs[['response_profs']],
                                 parameters[['Omega_response']], means_mat_r,
                                 cmi, as.numeric(parameters[['measurement_error_response']]))
      
      sse_mean <- c_compute_centered_obs_single(data_inputs[['Phi_response_prof']], data_inputs[['response_profs']],
                                                cmi, means_mat_r, as.numeric(parameters[['measurement_error_response']]))
      # Create the spatial covariance matrices
      E_step_list[['UTU']] = UTU = Matrix::forceSymmetric(E_step_list[['UTU']], uplo = 'L')
      
      if (is.null(data_inputs$time)){
        locs = data_inputs[['locs']]
      } else {
        locs = cbind(data_inputs[['locs']], time = data_inputs[['time']])
      }
      
      resp_scores_chol <- spatial_inv_covariance_scores(z = cluster_mat[,mc], 
                                                        n_pcs = params[['pc_response']],
                                                        G = params[['G']], 
                                                        variances = parameters[['variances_response']],
                                                        range_params = parameters[['range_params']],
                                                        model=ifelse(!is.null(params[['covariance_function']]), 'vecchia', 'independence'),
                                                        covfun_name = params[['covariance_function']],
                                                        locs = locs, 
                                                        reorder = !is.null(params[['covariance_function']]),
                                                        m = params[['m_train_response']])
      
      Linv_mat[[mc]] <- resp_scores_chol[['Linv_mats']]
      vecchia_info[[mc]] <- resp_scores_chol[['vecchia_info']]
      
      # reorder things back to the way it was before
      resp_scores_full <- Matrix::crossprod(resp_scores_chol[['Linv_all']])
      if (!is.null(params[['covariance_function']])){
        reorder_resp <- order(vecchia_info[[mc]][['order_overall']])
      } else {
        reorder_resp <- 1:nrow(resp_scores_full)
      }
      
      E_step_list[['A_inv']] <- resp_scores_full[reorder_resp, reorder_resp]
      
      # Compute relevant matrices/cholesky
      E_step_list[['Sigma_inv']] <- E_step_list[['A_inv']] + E_step_list[['UTU']]
      E_step_list[['Sigma_inv_chol']]  <-
        Matrix::Cholesky(E_step_list[['Sigma_inv']], perm = T, super = NA, LDL = F)
      
      # sample
      cond_mean = Matrix::solve(E_step_list[['Sigma_inv_chol']], UTX, system = 'A')
      cond_sample <- cond_mean + 
        Matrix::solve(E_step_list[['Sigma_inv_chol']], 
              Matrix::solve(E_step_list[['Sigma_inv_chol']], rnorm(length(cond_mean)),
                    system = 'Lt'), system = 'Pt')
      
      # compute likelihoods
      data_ll <- compute_likelihood_for_IS(sigma_chol = E_step_list[['Sigma_inv_chol']],
                                           A_inv_chol = resp_scores_chol[['Linv_all']],
                                           me_response = parameters[['measurement_error_response']],
                                           response_profs = data_inputs[['response_profs']],
                                           cond_mean = cond_mean,
                                           UTX = UTX,
                                           sse_mean = sse_mean)
      
    } else {
      # if no change in cluster memberships, only compute what you need
      cond_sample <- cond_mean + 
        Matrix::solve(E_step_list[['Sigma_inv_chol']], 
              Matrix::solve(E_step_list[['Sigma_inv_chol']], rnorm(length(cond_mean)),
                    system = 'Lt'), system = 'Pt')
      Linv_mat[[mc]] <- resp_scores_chol[['Linv_mats']]
      vecchia_info[[mc]] <- resp_scores_chol[['vecchia_info']]
    }
    
    normalizing_sums <- apply(like_mat, 1, function(x) max(x) + log(sum(exp(x - max(x)))))
    like_ind[mc] <- sum(like_mat[cbind(1:nrow(cluster_mat), cluster_mat[,mc])]) - 
      sum(normalizing_sums)
    like_spatial[mc] <- data_ll 
    
    pcs_array_response[,mc,] = matrix(cond_sample,
                                      nrow = data_inputs[['n_profiles']], byrow = T)
  }
  
  weights <- like_spatial - like_ind
  weights_normal <- vector('numeric', length(weights))
  for (i in 1:length(weights)) {
    weights_normal[i] <- 1/(sum(exp(weights - weights[i])))
  }
  wt <- weights_normal * MC_max
  
  if (params[['cluster_sampling']] %in% c('P_I', 'independence', 'gibbs')) {
    wt <- rep(1, params[['MC_max']])
  }
  
  if (params[['cluster_sampling']] %in% c('gibbs')) {
    cluster_mat <- cluster_mat[,-c(1:10)]
    pcs_array_response <- pcs_array_response[,-c(1:10),]
    Linv_mat <- Linv_mat[-c(1:10)]
    vecchia_info <- vecchia_info[-c(1:10)]
    like_ind <- like_ind[-c(1:10)]
    like_spatial <- like_spatial[-c(1:10)]
  }
  
  return(list('cluster_mat' = cluster_mat, 'pcs_array_response' = pcs_array_response,
              'p_mat' = p_mat, 
              'cond_mean' = cond_mean, 'likelihood_ind' = like_ind, 
              'Linv_mat' = Linv_mat, 'vecchia_info' = vecchia_info,
              'weights' = wt, 'likelihood_spat' = like_spatial))
}

Estep.mult_var_ind <- function(model, params, data_inputs, parameters){
  
  # Frequently used constants
  MC_max = params[['MC_max']]
  n_profiles = data_inputs[['n_profiles']]
  n_profiles_TS = data_inputs[['n_profiles_TS']]
  
  cluster_mat = matrix(0, ncol=MC_max, nrow=n_profiles_TS)
  # Matrices of principal component samples
  pcs_array_response = array(NA, c(n_profiles, MC_max, params[['pc_response']]))
  pcs_array_predictors = array(NA, c(n_profiles_TS, MC_max, params[['pc_predictors']]))
  # compute probabilities of being in each class

  p_mat <- return_probabilities(vars_r = parameters[['variances_response']],
                                response_profs = data_inputs[['response_profs']],
                                is_bgc = data_inputs[['BGC']],
                                phi_resp = data_inputs[['Phi_response_prof']],
                                means_resp = parameters[['means_response']],
                                me_resp = parameters[['measurement_error_response']],
                                Omegas_r1 = parameters[['Omega_response']],
                                G = params[['G']],
                                n_profiles = data_inputs[['n_profiles_TS']],
                                vars_p = parameters[['variances_predictors']],
                                predictor_profs = data_inputs[['predictor_profs']],
                                phi_pred = data_inputs[['Phi_predictors_prof']],
                                means_pred = parameters[['means_predictors']],
                                me_pred = parameters[['measurement_error_predictors']],
                                Omegas_pred = parameters[['Omega_predictors']],
                                profile_lengths_p = data_inputs[['profile_lengths_predictors']],
                                n_preds = length(params[['X']]), 
                                Lambdas = parameters[['Lambdas']])
  
  lik_mat = p_mat$lik_mat; p_mat = p_mat$p_mat
  
  # In fact a list of the parameters determining the distributions at each profiles for each cluster
  compute_conditional_distributions(response_prof = data_inputs$response_profs, 
                                    predictors_prof = data_inputs$predictor_profs, 
                                    phi_resp =  data_inputs$Phi_response_prof,
                                    phi_pred = data_inputs$Phi_predictors_prof, 
                                    phi_x_phi_resp = data_inputs$Phi_x_Phi_response, 
                                    phi_x_phi_pred = data_inputs$Phi_x_Phi_predictors,
                                    n_samples = as.integer(data_inputs$n_profiles_TS),
                                    means_resp = parameters$means_response, 
                                    means_pred = parameters$means_predictors,
                                    Omegas_resp = parameters$Omega_response,
                                    Omegas_pred = parameters$Omega_predictors,
                                    Lambdas = parameters$Lambdas,
                                    Sigma_eta_inv = parameters$Sigma_eta_inv,
                                    me_resp = parameters$measurement_error_response,
                                    me_pred = parameters$measurement_error_predictors,
                                    vars_resp = parameters$variances_response,
                                    vars_pred = parameters$variances_predictors,
                                    basis_lengths_p = data_inputs$nbasis_predictors,
                                    cond_probs =  p_mat,
                                    is_bgc = data_inputs$BGC,
                                    conditional_distributions = parameters$conditional_distributions)
  
  for (mc in 1:MC_max){
    cluster_mat[,mc] = sample_memberships(n_profiles_TS, params[['G']], p_mat)
    sample = create_scores_sample(cluster_mat[,mc], n_profiles_TS, params[['pc_response']] + params[['pc_predictors']], parameters$conditional_distributions)
    pcs_array_response[,mc,] = sample[,1:params[['pc_response']]][data_inputs[['BGC']], , drop = F]
    pcs_array_predictors[,mc,] = sample[,(params[['pc_response']]+1):(params[['pc_response']] + params[['pc_predictors']])]
  }
  
  return(list(cluster_mat = cluster_mat,
              pcs_array_response = pcs_array_response,
              pcs_array_predictors = pcs_array_predictors,
              p_mat = p_mat))
}

Estep.mult_var_spat <- function(model, params, data_inputs, parameters){
  
  # Frequently used constants
  G = params[['G']]
  MC_max = params[['MC_max']]
  n_profiles = data_inputs[['n_profiles']]
  n_profiles_TS = data_inputs[['n_profiles_TS']]
  
  cluster_mat = matrix(0, ncol=MC_max, nrow=n_profiles_TS)
  # Matrices of principal component samples
  pcs_array_response = array(NA, c(n_profiles, MC_max, params[['pc_response2']]))
  pcs_array_predictors = array(NA, c(n_profiles_TS, MC_max, params[['pc_predictors']]))
  
  like_spatial <- rep(0, MC_max)
  like_ind <- rep(0, MC_max)

  probs <- return_probabilities(vars_r = parameters[['variances_response']],
                                response_profs = data_inputs[['response_profs']],
                                is_bgc = data_inputs[['BGC']],
                                phi_resp = data_inputs[['Phi_response_prof']],
                                means_resp = parameters[['means_response']],
                                me_resp = parameters[['measurement_error_response']],
                                Omegas_r1 = parameters[['Omega_response1']],
                                G = G,
                                n_profiles = data_inputs[['n_profiles_TS']],
                                Omegas_r2 = parameters[['Omega_response2']],
                                vars_p = parameters[['variances_predictors']],
                                predictor_profs = data_inputs[['predictor_profs']],
                                phi_pred = data_inputs[['Phi_predictors_prof']],
                                means_pred = parameters[['means_predictors']],
                                me_pred = parameters[['measurement_error_predictors']],
                                Omegas_pred = parameters[['Omega_predictors']],
                                profile_lengths_p = data_inputs[['profile_lengths_predictors']],
                                n_preds = length(params[['X']]))
  
  p_mat <- probs[[1]]; like_mat <- probs[[2]]
  
  # save vecchia information so that we don't have to redo it in update_variances/range_params
  Linv_mat <- list()
  Linv_mat[['predictors']] <- list()
  Linv_mat[['response']] <- list()
  vecchia_info <- list()
  vecchia_info[['predictors']] <- list()
  vecchia_info[['response']] <- list()
  E_step_list = list()
  if (is.null(params[['use_MRF']])) {
    cluster_mat = rPotts(G=G,
                         neighbors = data_inputs[['nn_list']],
                         p_mat = p_mat,
                         likelihoods = like_mat,
                         theta = parameters[['theta']],
                         n_samples = MC_max)
  } else if (params[['use_MRF']]) {
    cluster_mat = rPotts(G=G,
                         neighbors = data_inputs[['nn_list']],
                         p_mat = p_mat,
                         likelihoods = like_mat,
                         theta = parameters[['theta']],
                         n_samples = MC_max)
  } else {
    cluster_mat <- matrix(nrow = n_profiles_TS, ncol = MC_max)
    for (mc in 1:MC_max) {
      for (x in 1:n_profiles_TS) {
        cluster_mat[x,mc] <- sample(1:G, size = 1, prob = p_mat[x,])
      }
    }
  }
  
  for (mc in 1:MC_max){

    cmi = cluster_mat[,mc] -1
    mode(cmi) ='integer'
    
    if (mc == 1 || !all(cluster_mat[,mc] == cluster_mat[,mc - 1])){
      if (mc > 1 & length(data_inputs[['Phi_x_Phi_predictors']]) > 50000) {
        rm(E_step_list)
        E_step_list <- list()
        gc()
      }
      
      E_step_list[['UTU']] =  c_compute_UTU(data_inputs[['Phi_x_Phi_response']],
                                            data_inputs[['Phi_x_Phi_predictors']],
                                            parameters[['Omega_response1']],
                                            parameters[['Omega_response2']],
                                            parameters[['Omega_predictors']],
                                            as.numeric(parameters[['measurement_error_response']]),
                                            as.numeric(parameters[['measurement_error_predictors']]),
                                            as.integer(data_inputs[['nbasis_predictors']]), 
                                            cmi, 
                                            data_inputs[['n_profiles']],
                                            data_inputs[['n_profiles_TS']],
                                            params[['G']],
                                            as.integer(data_inputs[['BGC']]))
      
      
      means_mat_r = as.matrix(do.call(cbind, parameters[['means_response']]))
      means_mat_p = as.matrix(do.call(cbind, parameters[['means_predictors']]))
      UTX = c_compute_UTX(data_inputs[['Phi_response_prof']], data_inputs[['Phi_predictors_prof']],
                          data_inputs[['response_profs']], data_inputs[['predictor_profs']],
                          parameters[['Omega_response1']], parameters[['Omega_response2']], parameters[['Omega_predictors']],
                          means_mat_r, means_mat_p,
                          cmi, as.numeric(parameters[['measurement_error_response']]),
                          as.numeric(parameters[['measurement_error_predictors']]),
                          as.integer(data_inputs[['nbasis_predictors']]), 
                          as.integer(data_inputs[['BGC']]))
      sse_mean <- c_compute_centered_obs(data_inputs[['Phi_response_prof']],
                                         data_inputs[['Phi_predictors_prof']],
                                         lapply(data_inputs[['response_profs']], function(x) x/sqrt(parameters[['measurement_error_response']])),
                                         lapply(1:length(data_inputs[['predictor_profs']]), function(x) data_inputs[['predictor_profs']][[x]]/rep(sqrt(parameters[['measurement_error_predictors']]),
                                                                                                                                                  times = data_inputs[['profile_lengths_predictors']][x,])),
                                         means_mat_r/sqrt(parameters[['measurement_error_response']]),
                                         means_mat_p/rep(sqrt(parameters[['measurement_error_predictors']]),
                                                         times = data_inputs[['nbasis_predictors']]), 
                                         cmi, 
                                         as.integer(data_inputs[['nbasis_predictors']]), 
                                         as.integer(data_inputs[['BGC']]))
      
      # Create the spatial covariance matrices
      cluster_mat_BGC <- cluster_mat[data_inputs[['BGC']],mc]
      E_step_list[['UTU']] = UTU = Matrix::forceSymmetric(E_step_list[['UTU']], uplo = 'L')
      
      if (is.null(data_inputs$time)){
        locs = data_inputs[['locs']]
      } else {
        locs = cbind(data_inputs[['locs']], time = data_inputs[['time']])
      }
      
      pred_scores_chol <- spatial_inv_covariance_scores(z = cluster_mat[,mc], n_pcs = params[['pc_predictors']], 
                                                        G = params[['G']], variances = parameters[['variances_predictors']],
                                                        range_params = parameters[['range_params_p']],
                                                        model=ifelse(!is.null(params[['covariance_function']]), 'vecchia', 'independence'),
                                                        covfun_name = params[['covariance_function']], 
                                                        locs = locs,
                                                        m = params[['m_train_predictors']])
      Linv_mat[['predictors']][[mc]] <- pred_scores_chol[['Linv_mats']]
      vecchia_info[['predictors']][[mc]] <- pred_scores_chol[['vecchia_info']]
      
      resp_scores_chol <- spatial_inv_covariance_scores(cluster_mat_BGC, n_pcs = params[['pc_response2']],
                                                        G = params[['G']], variances = parameters[['variances_response']],
                                                        range_params = parameters[['range_params_r']],
                                                        model=ifelse(!is.null(params[['covariance_function']]), 'vecchia', 'independence'),
                                                        covfun_name = params[['covariance_function']],
                                                        locs = locs[data_inputs[['BGC']],],
                                                        m = params[['m_train_response']])
      Linv_mat[['response']][[mc]] <- resp_scores_chol[['Linv_mats']]
      vecchia_info[['response']][[mc]] <- resp_scores_chol[['vecchia_info']]
      
      # reorder things back to the way it was before
      pred_scores_full <- Matrix::crossprod(pred_scores_chol[['Linv_all']])
      if (!is.null(params[['covariance_function']])) {
        reorder_pred <- order(vecchia_info[['predictors']][[mc]][['order_overall']])
      } else {
        reorder_pred <- 1:nrow(pred_scores_full)
      }
      pred_scores_full_reordered <- pred_scores_full[reorder_pred, reorder_pred]
      
      resp_scores_full <- Matrix::crossprod(resp_scores_chol[['Linv_all']])
      if (!is.null(params[['covariance_function']])) {
        reorder_resp <- order(vecchia_info[['response']][[mc]][['order_overall']])
      } else {
        reorder_resp <- 1:nrow(resp_scores_full)
      }
      resp_scores_full_reordered <- resp_scores_full[reorder_resp, reorder_resp]
      
      # create this matrix, we will need its determinant, but it is not reordered
      E_step_list[['A_inv_chol_wrong_order']] <- Matrix::bdiag(pred_scores_chol[['Linv_all']], 
                                                       resp_scores_chol[['Linv_all']])
      
      # Compute relevant matrices/cholesky
      E_step_list[['A_inv']] <- Matrix::bdiag(pred_scores_full_reordered,
                                      resp_scores_full_reordered)
      E_step_list[['Sigma_inv']] <- E_step_list[['A_inv']] + E_step_list[['UTU']]

      E_step_list[['Sigma_inv_chol']]  <-
        Matrix::Cholesky(E_step_list[['Sigma_inv']], perm = T, super = NA, LDL = F)
      
      # sample
      cond_mean = Matrix::solve(E_step_list[['Sigma_inv_chol']], UTX,
                        system = 'A')
      
      cond_sample <- cond_mean + 
        Matrix::solve(E_step_list[['Sigma_inv_chol']], 
              Matrix::solve(E_step_list[['Sigma_inv_chol']], rnorm(length(cond_mean)),
                    system = 'Lt'), system = 'Pt')
      
      # compute likelihoods
      data_ll <- compute_likelihood_for_IS(sigma_chol = E_step_list[['Sigma_inv_chol']],
                                           A_inv_chol = E_step_list[['A_inv_chol_wrong_order']], 
                                           me_response = parameters[['measurement_error_response']],
                                           response_profs = data_inputs[['response_profs']],
                                           cond_mean = cond_mean,
                                           UTX = UTX,
                                           sse_mean = sse_mean,
                                           predictor_profs = data_inputs[['predictor_profs']],
                                           profile_lengths_predictors = data_inputs[['profile_lengths_predictors']],
                                           me_predictors = parameters[['measurement_error_predictors']])
    } else {
      # if no change in cluster memberships, only conpute what you need
      cond_sample <- cond_mean + 
        Matrix::solve(E_step_list[['Sigma_inv_chol']], 
              Matrix::solve(E_step_list[['Sigma_inv_chol']], rnorm(length(cond_mean)),
                    system = 'Lt'), system = 'Pt')
      Linv_mat[['predictors']][[mc]] <- pred_scores_chol[[2]]
      vecchia_info[['predictors']][[mc]] <- pred_scores_chol[[3]]
      Linv_mat[['response']][[mc]] <- resp_scores_chol[[2]]
      vecchia_info[['response']][[mc]] <- resp_scores_chol[[3]]
    }
    normalizing_sums <- apply(like_mat, 1, function(x) max(x) + log(sum(exp(x - max(x)))))
    like_ind[mc] <- sum(like_mat[cbind(1:nrow(cluster_mat), cluster_mat[,mc])])- 
      sum(normalizing_sums)
    like_spatial[mc] <- data_ll
    
    start = (data_inputs[['n_profiles_TS']]*params[['pc_predictors']])
    pcs_array_predictors[,mc,] = matrix(cond_sample[1:start], nrow = data_inputs[['n_profiles_TS']], byrow = T)
    pcs_array_response[,mc,] = matrix(cond_sample[(start+1):(start + data_inputs[['n_profiles']]*params[['pc_response2']])],
                                      nrow = data_inputs[['n_profiles']], byrow = T)
  }
  weights <- like_spatial - like_ind
  weights_normal <- vector('numeric', length(weights))
  for (i in 1:length(weights)) {
    weights_normal[i] <- 1/(sum(exp(weights - weights[i])))
  }
  wt <- weights_normal * MC_max
  
  return(list('cluster_mat' = cluster_mat, 'pcs_array_response' = pcs_array_response,
              'pcs_array_predictors' = pcs_array_predictors, 'p_mat' = p_mat, 
              'cond_mean' = cond_mean, 'likelihood_ind' = like_ind, 
              'Linv_mat' = Linv_mat, 'vecchia_info' = vecchia_info, 'weights' = wt, 
              'likelihood_spat' = like_spatial))
}




