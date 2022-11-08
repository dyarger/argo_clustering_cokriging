compute_u_matrices <- function(Phi_x_Phi_resp, weights, cmi, cmi_BGC, n_e_step_samples, G,
                               pcs_r, Phi_x_Phi_pred = NULL, pcs_p1=NULL, pcs_p2=NULL, ind = F){
  
  u_mats = list()
  u_mats[['mean_resp']] = Matrix::bdiag(
    c_create_summed_U_matrix_sparse(Phi_x_Phi_resp, cmi_BGC, as.numeric(weights), as.integer(G),
                                    n_e_step_samples)
  )
  u_mats[['pc_resp']] = lapply(1:dim(pcs_r)[3], function(q){
    pc_weights_r = pcs_r[,,q]
    Matrix::bdiag(
      c_create_summed_U_matrix_pcs_sparse(Phi_x_Phi_resp, cmi_BGC, as.numeric(weights), pc_weights_r, as.integer(G),
                                        n_e_step_samples)
    )
  })
  if(!is.null(Phi_x_Phi_pred)){
    u_mats[['mean_pred']] = Matrix::bdiag(
      c_create_summed_U_matrix_sparse(Phi_x_Phi_pred, cmi, as.numeric(weights), as.integer(G),
                                      n_e_step_samples)
    )
    
    u_mats[['pc_pred']] = lapply(1:dim(pcs_p1)[3], function(q){
      pc_weights_p = pcs_p1[,,q]
      Matrix::bdiag(
        c_create_summed_U_matrix_pcs_sparse(Phi_x_Phi_pred, cmi, as.numeric(weights), pc_weights_p, as.integer(G),
                                            n_e_step_samples)
      )
    })
    if (!ind){
      u_mats[['lt']] = lapply(1:dim(pcs_p1)[3], function(q){
        pc_weights_p = pcs_p2[,,q]
        Matrix::bdiag(
          c_create_summed_U_matrix_pcs_sparse(Phi_x_Phi_resp, cmi_BGC, as.numeric(weights), pc_weights_p, as.integer(G),
                                              n_e_step_samples)
        )
      })
    }
  }
  return(u_mats)
}


update_means = function(U, Phi_prof, profs, cmi, pcs_array1, Omegas1, weights,
                        pen_mat, G=1, Omegas2=NULL, pcs_array2=NULL, lambda=10^-5){
  
  n_e_step_samples = as.integer(ncol(cmi))

  pcs1 = lapply(1:length(profs), function(x) matrix(pcs_array1[x,,], nrow = n_e_step_samples))
  
  pen_mat_all = Matrix::bdiag(replicate(G, pen_mat))
  
  if (is.null(Omegas2)){
    
    v_mat = c_create_summed_V_matrix_sparse_single(Phi_prof, profs, as.integer(G), Omegas1, cmi,
                                                   pcs1, n_e_step_samples, as.numeric(weights))

    V <- as.double(v_mat)
    
    res = as.double(Matrix::solve(U + lambda * pen_mat_all, V))
    
    return(split(res, rep(1:G, each=length(res)/G)))
    
  } else {
    pcs2 = lapply(1:length(profs), function(x) matrix(pcs_array2[x,,], nrow = n_e_step_samples))
    
    v_mat = c_create_summed_V_matrix_sparse(Phi_prof, profs, as.integer(G), Omegas1, Omegas2,
                                            cmi, pcs1, pcs2, n_e_step_samples, as.numeric(weights))

    V <- as.double(v_mat)
    
    res = as.double(Matrix::solve(U + lambda * pen_mat_all, V))
    
    return(split(res, rep(1:G, each=length(res)/G)))
  }
}

update_principal_components <- function(u_mats, cmi, pcs_r1, profs, Phi_prof, Omegas1,
                                        means, weights, pen_mat, G=1,  pcs_r2 = NULL,
                                        Omegas2=NULL, lambda=10^-5){
  
  # The sqeuence pcs/Omegas needs to match, i.e. if one updates the 
  # parameters$Omega_r1, then pcs_r1 = betas (Gamma %*% pcs_p)
  # Its always Omegas1 that is updated
  
  Q = dim(pcs_r1)[3]
  n_e_step_samples = as.integer(ncol(cmi))
  means_mat = as.matrix(do.call(cbind, means))
  
  pcs1 = lapply(1:length(profs), function(x) matrix(pcs_r1[x,,], nrow = n_e_step_samples))
  
  if (!is.null(Omegas2)){
    pcs2 = lapply(1:length(profs), function(x) matrix(pcs_r2[x,,], nrow = n_e_step_samples))
  }
  
  for (q in 1:Q) {
    
    pc_weights = pcs_r1[,,q]
    
    if(is.null(Omegas2)){
      v_mat = c_create_summed_V_matrix_pcs_sparse_single(Phi_prof, profs, means_mat, as.integer(G), Omegas1,
                                                         cmi, pcs1, n_e_step_samples, weights, as.integer(q-1))
    } else {
      v_mat = c_create_summed_V_matrix_pcs_sparse(Phi_prof, profs, means_mat, as.integer(G), Omegas1, Omegas2,
                                                    cmi, pcs1, pcs2, n_e_step_samples,
                                                    weights, as.integer(q-1))      
    }
    
    U <- u_mats[[q]]
    V <- as.double(v_mat)
    
    pen_mat_all = Matrix::bdiag(replicate(G, pen_mat))
    
    res = as.double(Matrix::solve(U + lambda * pen_mat_all, V))

    for (g in 1:G){
      Omegas1[[g]][,q] = split(res, rep(1:G, each = length(res)/G))[[g]]
    } 
  }
  return(Omegas1)
}


update_measurement_error <- function(cluster_mat, pcs_array1, profs, phi_prof, means, Omegas1,
                                     n_profiles, Omegas2=NULL, pcs_array2=NULL, n_vars = 1,
                                     weights = rep(1, ncol(cluster_mat)), G=1, prof_ns = NULL) {
  
  cluster_mat_i = cluster_mat - 1
  mode(cluster_mat_i) = 'integer'
  means_mat = as.matrix(do.call(cbind, means))
  
  if(!is.null(Omegas2)){

    pcs1 = lapply(1:n_profiles, function(x) {
      matrix(pcs_array1[x,,], nrow = ncol(cluster_mat))
    })
    
    pcs2 = lapply(1:n_profiles, function(x) {
      matrix(pcs_array2[x,,], nrow = ncol(cluster_mat))
    })
    
    c_update_measurement_error(cluster_mat_i, phi_prof, pcs1, pcs2, profs, means_mat, Omegas1, Omegas2, n_profiles, weights, G)
    
  } else {
    
    errors = rep(0, n_vars)
    MC <- ncol(cluster_mat)
    all_prof_lengths_equal <- all(apply(prof_ns, 1, function(x) {diff(range(x)) < .Machine$double.eps ^ 0.5}))
    for (x in 1:n_profiles){
      diff = c_compute_squared_sparse(cluster_mat_i[x,,drop= F], matrix(pcs_array1[x, ,], ncol = dim(pcs_array1)[3]),
                                      profs[[x]], as.integer(G), phi_prof[[x]], means_mat, Omegas1, weights)
      if (all_prof_lengths_equal) {
        errors <- errors + colSums(matrix(ncol = n_vars, diff))
      } else {
        vec_indexes <- as.double(
          sapply(1:length(prof_ns[x,]), function(i) rep(i, times=prof_ns[x,i]))
        )
        errors <- errors + tapply(diff, vec_indexes, sum)
      }
    }
    return(errors/(sum(weights)*colSums(prof_ns)))
  }
}

update_cluster_membership <- function(cluster_mat, n_profiles, G){
  conditional_probs = matrix(0, nrow=n_profiles, ncol=G)
  for(g in 1:G){
    conditional_probs[,g] = rowMeans(cluster_mat==g)
  }	
  cluster_membership = apply(conditional_probs, 1, which.max)
  return(list(cluster_membership, conditional_probs))
}

# Gauss Markov Random Field
# Optimize functions for MRF (pretty much taken from Liang et al)s
update_theta <- function(interval, nn_list, cluster_mat, n_profiles, G=1, theta, weights= rep(1, ncol(cluster_mat))){
  tryCatch(expr = {uniroot(gradient_mrf, interval, nn_list = nn_list, weights = weights,
                           cluster_mat = cluster_mat, tol = .Machine$double.eps^0.15, n = n_profiles, G = G)$root },
           error = function(x) {print('MRF optim failed');theta})
}

gradient_mrf <- function(theta, cluster_mat, nn_list, n_profiles, G, weights = rep(1, ncol(cluster_mat))) {
  MC <- ncol(cluster_mat)
  tmp <- matrix(0, n_profiles, G)
  for(x in 1:n_profiles){
    nn <- length(nn_list[[x]])
    for (mc in 1:MC) {
      neigh_clust <- cluster_mat[nn_list[[x]],mc]
      n_each <- tabulate(neigh_clust, nbins = G)
      denominator <-  exp(n_each * theta)
      numerator <- denominator * n_each
      if (sum(numerator) == Inf) {
        next
      }
      for(g in 1:G){
        tmp[x,g] <-tmp[x,g] + 
          weights[mc] * (cluster_mat[x,mc] == g)*(n_each[g] -sum(numerator)/sum(denominator))
      }
    }
  }
  return(sum(tmp))
}


update_variances <- function(cluster_mat, pcs_array, range_params,
                             G, n_pcs, MC, cov_fun, independent_model = F, 
                             Linv, variances, vecchia_info) {
  lapply(1:G, function(g) {
    sapply(1:n_pcs, function(q) {
      mean(sapply(1:MC, function(mc) {
        indexes <- cluster_mat[,mc] == g
        pcs_use <- as.double(pcs_array[indexes,mc,q])
        if (independent_model) {
          return(sum(pcs_use^2)/sum(indexes))
        } else {
          return(sum(as.double((Linv[[mc]][[g]][[q]] *sqrt(variances[[g]][[q]])) %*% pcs_use[vecchia_info[[mc]][['order']][[g]]])^2)/sum(indexes))
        }
      }))
    })
  })
}

# spatial range parameters
update_range_params <- function(G, n_pcs, variances, cluster_mat, pcs_array,
                                vecchia_info,
                                mc_weights = rep(1, ncol(cluster_mat)),
                                covfun_name = 'exponential_spheretime',
                                maxit = 10, start_params_previous = NULL) {
  if (covfun_name == 'exponential_spheretime_warp'){
    start_params <- log(c(.2, 25, exp(c(0,0,0,0,0))))
    lower_params <- log(c(.001, .5, exp(c(-4,-4,-4,-4,-4))))
    upper_params <- log(c(4, 800, exp(c(4,4,4,4,4))))
  } else {
    start_params <- log(c(.2, 25))
    lower_params <- log(c(.000001, .1))
    upper_params <-  log(c(4, 800))
  }
  no_cluster_variation <- all(apply(cluster_mat, 1, function(x) var(x) == 0))
  # for each cluster and principal component, use optim
  vals <- lapply(1:G, function(g) {
    sapply(1:n_pcs, function(q) {
      if (!is.null(start_params_previous)) {
        if (covfun_name == 'exponential_spheretime_warp') {
          start_params <- c(log(start_params_previous[[g]][1:2,q]),
                            start_params_previous[[g]][-c(1:2),q])
        } else {
          start_params <- log(start_params_previous[[g]][,q])
        }
      }
      # theta, var_score, cluster_mat, pcs_array, q, g2, MC_keep,
      # mc_weights, vecchia_info, covfun_name = 'exponential_spheretime',
      # no_cluster_variation
      return_vals <- optim(method = 'L-BFGS-B', range_optim_vecchia, par = start_params,
                           var_score = variances[[g]][q],
                           cluster_mat = cluster_mat, pcs_array = pcs_array, control = list(factr = 1e7, maxit = maxit),
                           mc_weights = mc_weights, 
                           q = q, g2 = g, MC_keep = dim(cluster_mat)[2], lower = lower_params, upper = upper_params,
                           covfun_name = covfun_name, 
                           vecchia_info = vecchia_info, no_cluster_variation = no_cluster_variation)$par
      return_vals[1:2] <- exp(return_vals[1:2])
      return_vals
    })
  })
  vals
}

# spatial range parameter for 1 pc and 1 cluster
range_optim_vecchia <- function(theta, var_score, cluster_mat, pcs_array, q, g2, MC_keep,
                                mc_weights, vecchia_info, covfun_name = 'exponential_spheretime',
                                no_cluster_variation) {
  scale <- exp(theta[1:2])
  if (covfun_name == 'exponential_spheretime') {
    cov_params_use <-  c(var_score, scale[1], scale[2], .00001)
  } else if (covfun_name == 'exponential_spheretime_warp') {
    cov_params_use <-  c(var_score, scale[1], scale[2], .00001, theta[3:7])
  } else if (covfun_name == 'matern_spheretime') {
    cov_params_use <-  c(var_score, scale[1], scale[2], .5, .00001)
  } else if (covfun_name == 'exponential_isotropic'){
    cov_params_use <-  c(var_score, scale[1], .00001)
  }
  
  like_comp <- tryCatch({
    
    val <- rep(0, MC_keep)
    for (mc in 1:MC_keep) {
      if (mc == 1) {
        Linv_temp <- GpGp::vecchia_Linv(cov_params_use, covfun_name = covfun_name,
                                        locs = vecchia_info[[mc]][['locs_list']][[g2]], 
                                        NNarray = vecchia_info[[mc]][['NN_list']][[g2]])
      } else if (no_cluster_variation) {
        
      } else if (!isTRUE(all.equal(cluster_mat[,mc], cluster_mat[,mc-1]))){
        Linv_temp <- GpGp::vecchia_Linv(cov_params_use, covfun_name = covfun_name,
                                        locs =vecchia_info[[mc]][['locs_list']][[g2]], 
                                        NNarray = vecchia_info[[mc]][['NN_list']][[g2]])
      }
      order_use <- vecchia_info[[mc]][['order']][[g2]]
      vec_to_sq <- GpGp::Linv_mult(Linv_temp, pcs_array[cluster_mat[,mc] == g2,mc,q][order_use],
                                   vecchia_info[[mc]][['NN_list']][[g2]])
      val[mc] <- mc_weights[mc]*(.5 *  sum(vec_to_sq^2) -
                                   sum(log(Linv_temp[,1])) + .5*nrow(Linv_temp) * log(2 * pi))
    }
    val}, 
    error = function(z) {
      message(z)
      NA
    }) 
  (as.double(0.5 * sum(like_comp)))
}

update_Lambda <- function(cluster_mat, response_pcs_array, predictor_pcs_array, G=1){
  # Transform to Matrix so that indexing is easier
  nrows = dim(response_pcs_array)[1] * dim(response_pcs_array)[2]
  response_pcs_array = matrix(response_pcs_array, nrow = nrows, ncol = dim(response_pcs_array)[3])
  nrows = dim(predictor_pcs_array)[1] * dim(predictor_pcs_array)[2]
  predictor_pcs_array = matrix(predictor_pcs_array, nrow = nrows, ncol = dim(predictor_pcs_array)[3])
  lapply(1:G, function(g){
    ind <- as.vector(cluster_mat == g)
    cross = crossprod(response_pcs_array[ind,], predictor_pcs_array[ind,])
    predictor_matrix = crossprod(predictor_pcs_array[ind,])
    tryCatch({ t(solve(predictor_matrix, t(cross)))},
             error = function(x) {
               t(solve(predictor_matrix + diag(nrow(predictor_matrix)), t(cross)))})
  })
}

orthogonalize <- function(Omegas1, variances1, inner_prod_mat1, G,
                          Omegas2=NULL, variances2=NULL, inner_prod_mat2=NULL,
                          Lambdas=NULL, cluster_mat=NULL, pcs1=NULL, pcs2=NULL){
  if(!is.null(cluster_mat)){
    ret = list()
    variances_response = list()
    variances_predictors = list()
    Sigma_eta_inv = list()
    # Transform to Matrix so that indexing is easier
    n_pcs_response = dim(pcs1)[3]
    n_pcs_predictors = dim(pcs2)[3]
    nrows = dim(pcs1)[1] * dim(pcs1)[2]
    response_pcs_array = matrix(pcs1, nrow = nrows, ncol = n_pcs_response)
    nrows = dim(pcs2)[1] * dim(pcs2)[2]
    predictors_pcs_array = matrix(pcs2, nrow = nrows, ncol = n_pcs_predictors)
    
    for (g in 1:G){
      ind <- as.vector(cluster_mat == g)
      #Should the mean be estimated as well?
      cov_estimate_response = 1/sum(ind)*crossprod(response_pcs_array[ind,])
      A = inner_prod_mat1 %*% Omegas1[[g]] %*%
        tcrossprod(cov_estimate_response, Omegas1[[g]])  %*% t(inner_prod_mat1)
      cov_estimate_predictors = 1/sum(ind)*crossprod(predictors_pcs_array[ind,])
      B = inner_prod_mat2 %*% Omegas2[[g]] %*%
        tcrossprod(cov_estimate_predictors, Omegas2[[g]])%*% 
        t(inner_prod_mat2)
      eig_A = eigen(A, symmetric = T)
      eig_B = eigen(B, symmetric = T)
      Lambdas[[g]] = crossprod(make_first_entry_positive(eig_A$vectors[,1:n_pcs_response]),
                               inner_prod_mat1 %*% Omegas1[[g]]) %*% Lambdas[[g]] %*% 
        solve(crossprod(make_first_entry_positive(eig_B$vectors[,1:n_pcs_predictors]), inner_prod_mat2 %*% Omegas2[[g]]))
      
      Omegas1[[g]] = solve(inner_prod_mat1) %*% make_first_entry_positive(eig_A$vectors[,1:n_pcs_response])
      Omegas2[[g]] = solve(inner_prod_mat2) %*%  make_first_entry_positive(eig_B$vectors[, 1:n_pcs_predictors])
      variances_response[[g]] = eig_A$values[1:n_pcs_response]
      variances_predictors[[g]] = eig_B$values[1:n_pcs_predictors]
      Sigma_eta_inv[[g]] = solve(diag(variances_response[[g]]) - Lambdas[[g]] %*% tcrossprod(diag(variances_predictors[[g]]), Lambdas[[g]]), diag(1, nrow=length(variances_response[[g]]), ncol=length(variances_response[[g]])))
    }
    ret[[1]]=Omegas1
    ret[[2]]=Omegas2
    ret[[3]]=variances_response
    ret[[4]]=variances_predictors
    ret[[5]]=Lambdas
    ret[[6]]=Sigma_eta_inv
    return(ret)
  } else if (!is.null(Omegas2)){
    lapply(1:G, function(g) {
      A = tcrossprod(inner_prod_mat2 %*% Omegas2[[g]] %*% 
                       diag(variances2[[g]]),
                     Omegas2[[g]])  %*%  t(inner_prod_mat2)
      eig_A = eigen(A, symmetric = T)
      new_variances_r = eig_A[['values']][1:ncol(Omegas2[[g]])]
      new_Omega2 = solve(inner_prod_mat2) %*% make_first_entry_positive(eig_A[['vectors']][,1:ncol(Omegas2[[g]])])
      
      A = tcrossprod(inner_prod_mat1 %*% Omegas1[[g]] %*%
                       diag(variances1[[g]]),
                     Omegas1[[g]])  %*%
        t(inner_prod_mat1)
      eig_A = eigen(A, symmetric = T)
      new_variances_p <- eig_A[['values']][1:ncol(Omegas1[[g]])]
      O_alpha <- make_first_entry_positive(eig_A[['vectors']][,1:ncol(Omegas1[[g]])])
      new_Omega1 <- solve(inner_prod_mat1) %*% O_alpha
      new_lambda <- Lambdas[[g]] %*% 
        solve(t(O_alpha) %*% inner_prod_mat1 %*% Omegas1[[g]])
      list(new_variances_p, new_lambda, new_Omega1,
           new_variances_r, new_Omega2)
    })
  } else {
    lapply(1:G, function(g) {
      A = tcrossprod(inner_prod_mat1 %*% Omegas1[[g]] %*% 
                       diag(variances1[[g]]),
                     Omegas1[[g]])  %*%  t(inner_prod_mat1)
      eig_A = eigen(A, symmetric = T)
      new_variances_r <- eig_A[['values']][1:ncol(Omegas1[[g]])]
      new_Omegas <- solve(inner_prod_mat1) %*% make_first_entry_positive(eig_A[['vectors']][,1:ncol(Omegas1[[g]])])
      list(new_variances_r, new_Omegas)
    })
  }
}

compute_marginal_probabilities <- function(cluster_membership, neighbors, theta, n_profiles, G) {
  marginal_membership <- matrix(NA, length(neighbors), G)
  cluster_mat <- matrix(0, length(cluster_membership), G)
  cluster_mat[cbind(1:length(cluster_membership), 
                    cluster_membership)] <- 1
  for(x in 1:length(neighbors)){
    tmp = apply(cluster_mat[neighbors[[x]],,drop=F], 2,
                function(t){exp(theta*sum(t, na.rm = T))})
    if (sum(tmp) == 0) {
      marginal_membership[x,] <- rep(1/G, G)
    } else {
      marginal_membership[x,] <- tmp/sum(tmp)
    }
  } 
  marginal_membership
}



