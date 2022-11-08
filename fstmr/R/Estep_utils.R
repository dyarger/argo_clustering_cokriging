return_probabilities <- function(vars_r, response_profs, is_bgc, phi_resp, means_resp,
                                 me_resp, Omegas_r1, G, n_profiles, Omegas_r2=NULL, vars_p=NULL,
                                 predictor_profs=NULL, phi_pred=NULL, means_pred=NULL, me_pred=NULL,
                                 Omegas_pred=NULL, profile_lengths_p=NULL, n_preds=NULL, Lambdas=NULL){
  
  means_mat_resp = as.matrix(do.call(cbind, means_resp))
  variances =  me_resp
  
  if (!is.null(Omegas_r2)){
    vars_r = as.matrix(do.call(cbind, vars_r))
    means_mat_pred = as.matrix(do.call(cbind, means_pred))
    variances = c(me_pred, variances)
    vars_p = as.matrix(do.call(cbind, vars_p))
    mode(profile_lengths_p) = 'integer'

    lik_mat = c_compute_E_step_likelihoods(response_profs, predictor_profs, as.integer(is_bgc),
                                           as.integer(n_profiles), phi_resp, phi_pred,
                                           Omegas_r1, Omegas_r2, Omegas_pred, means_mat_resp,
                                           means_mat_pred, variances, vars_r, vars_p, profile_lengths_p,
                                           as.integer(G), as.integer(n_preds))
  } else if (is.null(Omegas_pred)) {
    vars_r = as.matrix(do.call(cbind, vars_r))
    lik_mat = c_compute_E_step_likelihoods_single(response_profs, as.integer(n_profiles), phi_resp,
                                                  Omegas_r1, means_mat_resp, variances, vars_r, as.integer(G))
  } else if (!is.null(Lambdas)){
    
    Gamma_list = form_joint_cov_scores(vars_r, vars_p, Lambdas, G)
    Gamma_list = lapply(Gamma_list, function(g) as.matrix(g))
    
    means_mat_pred = as.matrix(do.call(cbind, means_pred))
    variances = c(variances, me_pred)
    vars_p = as.matrix(do.call(cbind, vars_p))
    mode(profile_lengths_p) = 'integer'
    vars_r = as.matrix(do.call(cbind, vars_r))
    
    lik_mat = c_compute_E_step_likelihoods_ind(profs_resp = response_profs,
                                               profs_pred = predictor_profs,
                                               is_bgc = as.integer(is_bgc),
                                               n_profiles = n_profiles,
                                               basis_evals_resp = phi_resp,
                                               basis_evals_pred = phi_pred,
                                               Gammas = Gamma_list,
                                               Omegas_resp = Omegas_r1,
                                               Omegas_preds = Omegas_pred,
                                               means_resp = means_mat_resp,
                                               means_pred = means_mat_pred,
                                               variances = variances,
                                               vars_pred = vars_p,
                                               profile_lengths_p = profile_lengths_p,
                                               G = as.integer(G),
                                               n_preds = n_preds)
  }

  p_mat <- lik_mat
  
  for (x in 1:nrow(p_mat)){
    p_mat[x,] = log_like_to_probs(p_mat[x,], G)
  }
  return(list('p_mat' = p_mat, 'lik_mat' =  lik_mat))
}

spatial_inv_covariance_scores <- function(z, n_pcs, G, variances, range_params,
                                          model='dependence', locs,covfun_name = NULL,
                                          m = 25, reorder = T) {
  if (is.null(covfun_name)){
    covfun_name <- 'exponential_spheretime'
  }
  n <- length(z)
  if (model == 'independence') {
    reorder <- F
  }
  NN_list <- list()
  gamma_list <- list()
  gamma_list_mat <- list()
  locs_list <- list()
  order_list <- list()
  order_list_df <- list()
  for (g in 1:G) { 
    indexes <- which(z == g)
    if (length(indexes) == 0)  {
      next
    }
    locs_use <- locs[indexes,,drop = F]
    if (reorder & length(indexes) > 1) {
      if (covfun_name=='exponential_isotropic'){
        ord = GpGp::order_maxmin(locs_use, lonlat = T)
      } else {
        ord <- GpGp::order_maxmin(locs_use, lonlat = T, space_time = T, st_scale = NULL)
      }
    } else {
      ord <- 1:nrow(locs_use)
    }
    order_list[[g]] <- ord
    order_list_df[[g]] <- list()
    locs_ord <- locs_use[ord,,drop = F]
    if (length(indexes) != 1) {
      NN_list[[g]] <- GpGp::find_ordered_nn(locs_ord, m = m, lonlat = T)
    } else {
      NN_list[[g]] <- matrix(nrow = 1,ncol = m+1, c(1, rep(NA, m)))
    }
    locs_list[[g]] <- locs_ord
    gamma_list[[g]] <- list()
    gamma_list_mat[[g]] <- list()
    for (q in 1:n_pcs) {
      reorder_vec <- seq(q, n_pcs*n, by=n_pcs)[indexes]
      if (model == 'independence' | length(indexes) == 1) {
        gamma_list[[g]][[q]] <- list('i' = reorder_vec,
                                     'j' = reorder_vec,
                                     'x' = rep(1/sqrt(variances[[g]][q]),
                                               length(reorder_vec)))
        gamma_list_mat <- NN_list <- locs_list <- order_list_df <- order_list <-  NULL
      } else if (model == 'vecchia') {
        rp_one <- range_params[[g]][,q]
        if (covfun_name == 'exponential_spheretime') {
          cov_params_use <-  c(variances[[g]][q], rp_one[1], rp_one[2], .00001)
        } else if (covfun_name == 'exponential_spheretime_warp') {
          cov_params_use <-  c(variances[[g]][q], rp_one[1], rp_one[2], .00001, rp_one[3:7])
        } else if (covfun_name == 'matern_spheretime') {
          cov_params_use <-  c(variances[[g]][q], rp_one[1], rp_one[2], .5, .00001)
        } else if (covfun_name == 'exponential_isotropic'){
          cov_params_use <-  c(variances[[g]][q], rp_one[1], .00001)
        }
        # form inverse cholesky entries and then it as a sparse Matrix
        L_inv <- GpGp::vecchia_Linv(covparms = cov_params_use,
                                     covfun_name = covfun_name, 
                                     locs  = locs_ord, NNarray = NN_list[[g]])
        i_vec <- rep(1:length(indexes), each = ncol(NN_list[[g]]))
        j_vec <- as.double(t(NN_list[[g]]))
        x_vec <- as.double(t(L_inv))
        gamma_list[[g]][[q]] <- list('i' = reorder_vec[i_vec[!is.na(j_vec)]],
                                     'j' = reorder_vec[j_vec[!is.na(j_vec)]],
                                     'x' = x_vec[!is.na(j_vec)])
        gamma_list_mat[[g]][[q]] <- Matrix::sparseMatrix(i = i_vec[!is.na(j_vec)],
                                                 j = j_vec[!is.na(j_vec)],
                                                 x = gamma_list[[g]][[q]][['x']])
        order_list_df[[g]][[q]] <- data.frame(matrix_indexes = reorder_vec,
                                              original_index = 1:length(reorder_vec),
                                              ord,
                                              matrix_indexes_ord = reorder_vec[ord])
      }
    }
  }
  
  if (sum(sapply(gamma_list, length)) > 1 | G == 1) {
    gamma_df <- dplyr::bind_rows(gamma_list)
    gamma_tilde <- Matrix::sparseMatrix(i = gamma_df$i,
                                        j = gamma_df$j,
                                        x = gamma_df$x, dims = c(n_pcs *n, n_pcs*n))
  }
  if (model == 'vecchia') {
    order_list_df <- dplyr::bind_rows(order_list_df)
    order_list_df <- order_list_df[['matrix_indexes_ord']][order(order_list_df[['matrix_indexes']])]
  }
  vecchia_info <- list('NN_list' = NN_list, 'locs_list' = locs_list, 'order' = order_list,
                       'order_overall' = order_list_df)
  list('Linv_all' = gamma_tilde, 
       'Linv_mats' = gamma_list_mat, 
       'vecchia_info' = vecchia_info)
}

compute_likelihood_for_IS <- function(sigma_chol, A_inv_chol, me_response, response_profs, cond_mean, UTX, sse_mean,
                                      predictor_profs=NULL, profile_lengths_predictors=NULL, me_predictors){
  n_length <- sum(sapply(response_profs, length))
  rss <- sse_mean - sum(as.double(cond_mean) * UTX)
  
  if(!is.null(predictor_profs)){
    pred_length <- sum(sapply(predictor_profs, length))/ncol(profile_lengths_predictors)
    
    log_det = as.double(2*Matrix::determinant(sigma_chol)[['modulus']] - 2 * sum(log(Matrix::diag(A_inv_chol))) + 
                          n_length * log(me_response) + 
                          pred_length * sum(log(me_predictors)))
    
    n_length = n_length + ncol(profile_lengths_predictors)*pred_length
  } else {
    log_det = as.double(2*Matrix::determinant(sigma_chol)[['modulus']] - 2*sum(log(Matrix::diag(A_inv_chol))) + 
                          n_length * log(me_response))
  }

  return(-0.5 * n_length * log(2*pi)  + -0.5*log_det + -0.5 * rss)
}

compute_probabilities_gibbs <- function(scores_sample, 
                                        marginal_membership, n_profiles, n_scores, G,
                                        parameters,data_inputs){
  p_mat <- matrix(NA, nrow = n_profiles, ncol = G)
  for (x in 1:n_profiles){
    likelihoods =  sapply(1:G, function(g) {
      as.double(-0.5*sum((data_inputs[['response_profs']][[x]] -
                            data_inputs[['Phi_response_prof']][[x]] %*% 
                            parameters[['means_response']][[g]] -
                            as.double(data_inputs[['Phi_response_prof']][[x]] %*% 
                                        parameters[['Omega_response']][[g]] %*% scores_sample[x,] )
      )^2))
    }) 
    p_mat[x,] <- log_like_to_probs(likelihoods, marg = marginal_membership[x,], G = G)
  }
  return(p_mat)
}

form_joint_cov_scores <- function(variances_response, variances_predictors, Lambdas, G) {
  Gamma_list <- list()
  for (g in 1:G) {
    D_response <- Matrix::Diagonal(n= length(variances_response[[g]]),
                           variances_response[[g]])
    D_predictors <- Matrix::Diagonal(n= length(variances_predictors[[g]]),
                             variances_predictors[[g]])
    D_cc <- Matrix::tcrossprod(D_predictors, Lambdas[[g]])
    Gamma_list[[g]] <- rbind(cbind(D_response, Matrix::t(D_cc)),
                             cbind(D_cc,  D_predictors))
  }
  Gamma_list
}

create_scores_sample <- function(cluster_mat, n_profiles, n_scores, conditional_distributions){
  sample = matrix(NA, nrow=n_profiles, ncol=n_scores)
  for (x in 1:n_profiles){
    n_true <- length(conditional_distributions[[x]][[cluster_mat[x]]][[3]])
    sample[x,(n_scores - n_true + 1):n_scores] = as.numeric(conditional_distributions[[x]][[cluster_mat[x]]][[3]] + 
                                                              as.double(conditional_distributions[[x]][[cluster_mat[x]]][[2]] %*% 
                                                                          rnorm(nrow(conditional_distributions[[x]][[cluster_mat[x]]][[2]]))))
  }  
  return(sample)
}

compute_conditional_distributions <- function(response_prof, predictors_prof, phi_resp, phi_pred,
                                              phi_x_phi_resp, phi_x_phi_pred, n_samples, means_resp,
                                              means_pred, Omegas_resp, Omegas_pred, Lambdas, Sigma_eta_inv,
                                              me_resp, me_pred, vars_resp, vars_pred, basis_lengths_p,
                                              cond_probs, is_bgc, conditional_distributions){
  
  phi_dense = lapply(phi_resp, function(x) as.matrix(x))
  phi_x_phi_resp_dense = lapply(phi_x_phi_resp, function(x) as.matrix(x))
  phi_x_phi_pred = lapply(phi_x_phi_pred, function(x) as(x, 'dgCMatrix'))
  means_mat_resp = as.matrix(do.call(cbind, means_resp))
  means_mat_pred = as.matrix(do.call(cbind, means_pred))
  vars_resp = as.matrix(do.call(cbind, vars_resp))
  vars_pred = as.matrix(do.call(cbind, vars_pred))
  mode(basis_lengths_p) = 'integer'
  
  c_compute_conditional_distribution(response_prof, predictors_prof, phi_dense, phi_pred, phi_x_phi_resp_dense,
                                     phi_x_phi_pred, n_samples, means_mat_resp, means_mat_pred, Omegas_resp, Omegas_pred,
                                     Lambdas, Sigma_eta_inv, me_resp, me_pred, vars_resp, vars_pred, basis_lengths_p,
                                     cond_probs, as.integer(is_bgc), conditional_distributions)
}



log_like_to_probs <- function(log_like, G, marg = NULL, tolerance = 10^(-12)) {
  log_like[log_like == -Inf] <- min(log_like[log_like != -Inf]) - 10^8
  log_like <- log_like - mean(log_like)
  log_like_ord <- -stl_sort(log_like)
  cond_membership_vec <- rep(0, G)
  if (is.na(log_like_ord[1] - log_like_ord[2])) {
    return(rep(1/G, G))
  }
  if (((log_like_ord[1] - log_like_ord[2]) > abs(log(tolerance)))) {
    cond_membership_vec[which.max(log_like)] <- 1
    return(cond_membership_vec)
  }
  if (is.null(marg)) {
    s = sum(sapply(1:G, function(l) exp(log_like[[l]])))
    if (s == Inf) {
      cond_membership_vec[which.max(log_like)] <- 1
      return(cond_membership_vec)
    }
    if (is.na(s) | is.nan(s) | sum(is.na(log_like_ord)) > 
        0 | sum(is.nan(log_like_ord)) > 0) {
      return(rep(1/G, G))
    }
  } else {
    s = sum(sapply(1:G, function(l) exp(log_like[[l]]) * 
                     marg[l]))
    if (is.na(s) | is.nan(s) | sum(is.na(log_like_ord)) > 
        0 | sum(is.nan(log_like_ord)) > 0 | sum(is.na(marg)) > 
        0 | sum(is.nan(marg)) > 0) {
      return(rep(1/G, G))
    }
  }
  if (s == Inf) {
    cond_membership_vec[which.max(log_like)] <- 1
    return(cond_membership_vec)
  }
  for (g in 1:G) {
    if (!is.null(marg)) {
      cond_membership_vec[g] = exp(log_like[[g]]) * marg[g] / s
    } else {
      cond_membership_vec[g] = exp(log_like[[g]]) / s
    }
  }
  if (sum(is.nan(cond_membership_vec)) > 0 | (abs(sum(cond_membership_vec, 
                                                      na.rm = T) - 1) > 0.01)) {
    cond_membership_vec <- rep(1/G, G)
  }
  cond_membership_vec
}

rPotts <- function(G, neighbors, p_mat, likelihoods, theta, n_samples = 1, skip = 3){
  
  ngibbs = n_samples*skip+skip
  
  n = nrow(p_mat)
  samples = matrix(0, nrow = n, ncol = n_samples)
  potts_sample = rep(0, n)
  
  for(i in 1:n){
    potts_sample[i] = sample(1:G, size = 1, prob = p_mat[i,])
  }
  
  cluster_mat <- matrix(0, length(potts_sample), G)
  cluster_mat[cbind(1:length(potts_sample), 
                    potts_sample)] <- 1
  marginal_membership <- matrix(0, length(potts_sample), G)
  for(x in 1:length(neighbors)){
    tmp = apply(cluster_mat[neighbors[[x]],,drop=F], 2,
                function(t){exp(theta*sum(t, na.rm = T))})
    if (sum(tmp) == 0) {
      marginal_membership[x,] <- rep(1/G, G)
    } else {
      marginal_membership[x,] <- tmp/sum(tmp)
    }
  } 
  
  counter = 1
  
  for(i in 1:ngibbs){
    for(x in 1:n){
      prob = log_like_to_probs(likelihoods[x,], G, marg = marginal_membership)
      t_sample = sample(1:G, size = 1, prob = prob)
      if(t_sample != potts_sample[x]){
        cluster_mat[x,potts_sample[x]] = 0
        cluster_mat[x,t_sample] = 1
        potts_sample[x] = t_sample
        for (z in neighbors[[x]]) {
          tmp = apply(cluster_mat[neighbors[[z]],,drop=F], 2,
                      function(t){exp(theta*sum(t, na.rm = T))})
          if (sum(tmp) == 0) {
            marginal_membership[z,] <- rep(1/G, G)
          } else {
            marginal_membership[z,] <- tmp/sum(tmp)
          }
        }
      }
    }
    
    if (i %% skip == 0){
      samples[,counter] = potts_sample
      counter = counter + 1
    }
    
    if (counter == n_samples+1){
      break
    }
  }
  return(samples)
}

sample_memberships <- function(n_profiles, G, p_mat){
  new_sample = rep(NA, n_profiles)
  for (x in 1:n_profiles){
    new_sample[x] = sample(1:G, 1, prob = p_mat[x,])
  }
  return(new_sample)
}



