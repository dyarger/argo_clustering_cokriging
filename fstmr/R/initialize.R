initialize <- function(model, params, data_inputs){
  UseMethod('initialize')
}

initialize.one_var <- function(model, params, data_inputs){
  parameters = list()
  
  var = params[['Y']]
  id = params[['id']]
  ind_name = params[['ind_name']]
  domain = params[['domain_Y']]
  n_profiles = length(data_inputs[['response_profs']])
  G = params[['G']]
  
  parameters[['cluster_mat']] = initialize_clusters(levels = params[['levels']],
                                                    n_profiles = n_profiles,
                                                    var = var,
                                                    G = G,
                                                    init_strategy = params[['init_strategy']],
                                                    resp_profs = data_inputs[['response_profs']],
                                                    ind_vals1 = data_inputs[['ind_response']],
                                                    init_cluster = params[['init_cluster']])
  
  cmi = parameters[['cluster_mat']] -1 
  mode(cmi) = 'integer'
  
  weights = 1
  
  response_pcs_array = array(0, c(data_inputs[['n_profiles']], 1, params[['pc_response']]))
  
  u_mats = compute_u_matrices(Phi_x_Phi_resp = data_inputs[['Phi_x_Phi_response']], 
                              weights = weights,
                              cmi = NULL,
                              cmi_BGC = cmi,
                              n_e_step_samples=ncol(cmi),
                              G = G,
                              pcs_r = response_pcs_array)
  
  parameters[['means_response']] = update_means(U = u_mats[['mean_resp']],
                                                Phi_prof = data_inputs[['Phi_response_prof']],
                                                profs = data_inputs[['response_profs']],
                                                cmi = cmi,
                                                pcs_array1 = response_pcs_array,
                                                Omegas1 = lapply(1:G, function(g){
                                                  matrix(0, ncol = params[['pc_response']], nrow = data_inputs[['nbasis_response']])
                                                }),
                                                weights = weights,
                                                pen_mat = data_inputs[['penalty_mat_response']],
                                                G = G,
                                                lambda = params[['lambda_mean_response']][1])
  
  parameters[['lambda_mean_response']] = params[['lambda_mean_response']][1]

  
  p_seq <- seq(domain[1], domain[2], by = (domain[2] - domain[1])/100)
  
  if (G == 2 & params[['pc_response']] == 2) {
    pcs_initial <- list(cbind(1, p_seq), cbind(rev(p_seq)^2, rev(p_seq)))
  } else {
    pcs_initial <- replicate(params[['G']], cbind(1, p_seq, p_seq^2, p_seq^3,
                                                  p_seq^4, p_seq^5, p_seq^6,
                                                  p_seq^7, p_seq^8, p_seq^9,
                                                  p_seq^10, p_seq^11, p_seq^12, 
                                                  p_seq^13, p_seq^14, p_seq^15)[,1:params[['pc_response']]], simplify = F)
  }
  
  Phi <- fda::eval.basis(p_seq, params[['basis_response']])
  
  parameters[['Omega_response']] = lapply(1:params[['G']], function(g){
    make_first_entry_positive(Matrix::solve(data_inputs[['inner_prod_mat_response']]) %*% qr.Q(qr(Matrix::solve(Matrix::crossprod(Phi) + diag(.5,ncol(Phi) ), Matrix::crossprod(Phi, pcs_initial[[g]])))))
  })
  
  var_guess <- sqrt(mean(sapply(1:data_inputs[['n_profiles']], 
                                function(x) {
                                  mean(as.double(data_inputs[['response_profs']][[x]] - 
                                          data_inputs[['Phi_response_prof']][[x]] %*% 
                                          parameters[['means_response']][[parameters[['cluster_mat']][x]]])^2)
                                }))*2000)
  
  response_pcs_array = array(rnorm(data_inputs[['n_profiles']] * params$pc_response,
                                   sd = var_guess), 
                             dim = c(data_inputs[['n_profiles']], 1, 
                                     params$pc_response))
  
  u_mats = compute_u_matrices(Phi_x_Phi_resp = data_inputs[['Phi_x_Phi_response']], 
                              weights = 1,
                              cmi = NULL,
                              cmi_BGC = cmi,
                              n_e_step_samples=ncol(cmi),
                              G = G,
                              pcs_r = response_pcs_array)
    
  parameters[['Omega_response']] = update_principal_components(u_mats = u_mats[['pc_resp']],
                                                               cmi = cmi,
                                                               pcs_r1 = response_pcs_array,
                                                               profs = data_inputs[['response_profs']],
                                                               Phi_prof = data_inputs[['Phi_response_prof']],
                                                               Omegas1 = parameters[['Omega_response']],
                                                               means = parameters[['means_response']], 
                                                               weights = weights,
                                                               pen_mat = data_inputs[['penalty_mat_response']],
                                                               G = G,
                                                               lambda = params[['lambda_pcs_response']][1])
  
  parameters[['lambda_pcs_response']] = params[['lambda_pcs_response']][1]
  
  parameters[['Omega_response']] <- lapply(1:params[['G']], function(g) {
    make_first_entry_positive(qr.Q(qr(parameters[['Omega_response']][[g]])))
  })
  
  # Use least squares estimate for the scores in order to obtain decent start values for the measurement error variances
  pcs_array_response = array(0, c(data_inputs[['n_profiles']], 1, params[['pc_response']]))
  for (i in 1:data_inputs[['n_profiles']]){
    cluster = parameters[['cluster_mat']][[i]]
    a = data_inputs[['response_profs']][[i]]; b = data_inputs[['Phi_response_prof']][[i]] %*% parameters[['means_response']][[cluster]]
    c = data_inputs[['Phi_response_prof']][[i]] %*% parameters[['Omega_response']][[cluster]]
    pcs_array_response[i,1,] = as.double(Matrix::solve(Matrix::crossprod(c) + diag(x = .00001, ncol(c)),Matrix::crossprod(c, a - b)))
  }
  
  rp = data_inputs[['response_profs']]
  phi_prof =data_inputs[['Phi_response_prof']]
  var = 0
  for (i in 1:length(rp)){
    g = parameters[['cluster_mat']][i,1]
    centered = rp[[i]] - phi_prof[[i]] %*% (parameters$means_response[[g]] - parameters$Omega_response[[g]] %*% pcs_array_response[i,1,])
    var = var + sum(centered^2)/length(centered)
  }
  var/n_profiles
  
  parameters[['measurement_error_response']] = update_measurement_error(cluster_mat = parameters[['cluster_mat']],
                                                                        pcs_array1 = pcs_array_response,
                                                                        profs = data_inputs[['response_profs']],
                                                                        phi_prof =data_inputs[['Phi_response_prof']],
                                                                        means = parameters[['means_response']],
                                                                        Omegas1 = parameters[['Omega_response']],
                                                                        n_profiles = data_inputs[['n_profiles']],
                                                                        weights = rep(1, ncol(parameters[['cluster_mat']])),
                                                                        G=G,
                                                                        prof_ns = data_inputs[['profile_lengths_response']])
  
  parameters[['variances_response']] <- update_variances(cluster_mat = parameters[['cluster_mat']],
                                                         pcs_array = pcs_array_response,
                                                         range_params = parameters[['range_params']],
                                                         G = G,
                                                         n_pcs = dim(pcs_array_response)[3],
                                                         MC = ncol(parameters[['cluster_mat']]),
                                                         cov_fun = params[['cov_fun']],
                                                         independent_model = T,
                                                         Linv = NULL,
                                                         variances = NULL,
                                                         vecchia_info = NULL)
  # Markov random field
  parameters[['theta']] <- .2
  parameters[['marginal_probabilities']] <- compute_marginal_probabilities(parameters[['cluster_mat']], parameters[['theta']],
                                                                           neighbors = data_inputs[['nn_list']], data_inputs[['n_profiles']],
                                                                           params[['G']])
  
  if (is.null(params[['covariance_function']])) {
    n_range_params <- 2
  } else {
    n_range_params <- ifelse(params[['covariance_function']] == 'exponential_spheretime_warp',
                             7, 2)
  }
  vec_start <- c(.2, 50, 0, 0, 0, 0,0)
  parameters[['range_params']] <- lapply(1:params[['G']], function(d) {
    matrix(nrow = n_range_params, ncol = params[['pc_response']], vec_start[1:n_range_params])
  })
  
  if(!is.null(params[['lambda_mean_response']])){
    parameters[['lambda_mean_response']] = params[['lambda_mean_response']]
  } else {
    parameters[['lambda_mean_response']] = 10
  }

  return(parameters)
}

initialize.mult_var_ind <- function(model, params, data_inputs){

  if (is.null(data_inputs[['n_profiles_TS']])) {
    data_inputs[['n_profiles_TS']] = data_inputs[['n_profiles']]  
  }
  
  parameters = list()
  
  var = params[['Y']]
  preds = params[['X']]
  
  id = params[['id']]
  ind_name = params[['ind_name']]
  
  domain1 = params[['domain_Y']]
  domains2 = params[['domain_X']]
  
  n_profiles = data_inputs[['n_profiles_TS']]
  G = params[['G']]
  
  
  combined = c(params[['Y']], params[['preds']])
  
  parameters[['cluster_mat']] = initialize_clusters(n_profiles = n_profiles,
                                                    var = var, G = G,
                                                    init_strategy = params[['init_strategy']],
                                                    levels = params[['levels']],
                                                    resp_profs = data_inputs[['response_profs']],
                                                    ind_vals1 = data_inputs[['ind_response']],
                                                    preds=preds,
                                                    pred_profs=data_inputs[['predictor_profs']], 
                                                    ind_vals2=data_inputs[['ind_predictors']],
                                                    profile_lengths=data_inputs[['profile_lengths_predictors']],
                                                    equal=(data_inputs[['n_profiles_TS']]==data_inputs[['n_profiles']]),
                                                    init_cluster = params[['init_cluster']])
  
  cmi = parameters[['cluster_mat']] -1 
  mode(cmi) = 'integer'
  cmi_BGC = cmi[data_inputs[['BGC']], ,drop = F]
  mode(cmi_BGC) = 'integer'
  
  response_pcs_array = array(0, c(data_inputs[['n_profiles']], 1, params[['pc_response']]))
  predictor_pcs_array = array(0, c(data_inputs[['n_profiles_TS']], 1, params[['pc_predictors']]))
  predictor_pcs_array_BGC = predictor_pcs_array[data_inputs[['BGC']], , , drop = F]
  
  u_mats = compute_u_matrices(Phi_x_Phi_resp = data_inputs[['Phi_x_Phi_response']], 
                              weights = 1,
                              cmi = cmi,
                              cmi_BGC = cmi_BGC,
                              n_e_step_samples=ncol(cmi),
                              G = G,
                              pcs_r = response_pcs_array,
                              Phi_x_Phi_pred = data_inputs[['Phi_x_Phi_predictors']],
                              pcs_p1 = predictor_pcs_array,
                              pcs_p2 = predictor_pcs_array_BGC,
                              ind = T)
  
  
  parameters[['means_response']] = update_means(U = u_mats[['mean_resp']],
                                                Phi_prof = data_inputs[['Phi_response_prof']],
                                                profs = data_inputs[['response_profs']],
                                                cmi = cmi_BGC,
                                                pcs_array1 = response_pcs_array,
                                                Omegas1 = lapply(1:G, function(g){
                                                  matrix(0, ncol = params[['pc_response']], nrow = data_inputs[['nbasis_response']])
                                                }),
                                                weights = rep(1, ncol(parameters[['cluster_mat']])),
                                                pen_mat = data_inputs[['penalty_mat_response']],
                                                G = G,
                                                lambda = params[['lambda_mean_response']][1])
  
  parameters[['lambda_mean_response']] = params[['lambda_mean_response']][1]
  
  parameters[['means_predictors']] = update_means(U = u_mats[['mean_pred']],
                                                  Phi_prof = data_inputs[['Phi_predictors_prof']],
                                                  profs = data_inputs[['predictor_profs']],
                                                  cmi = cmi,
                                                  pcs_array1 = predictor_pcs_array,
                                                  Omegas1 = lapply(1:params[['G']], function(g) matrix(0, nrow = sum(data_inputs[['nbasis_predictors']]) , ncol = params[['pc_predictors']])),
                                                  weights = rep(1, ncol(parameters[['cluster_mat']])),
                                                  pen_mat = data_inputs[['penalty_mat_predictors']],
                                                  G = G,
                                                  lambda = params[['lambda_mean_predictors']][1])
  parameters[['lambda_mean_predictors']] = params[['lambda_mean_predictors']][1]
  
  parameters[['Omega_response']] = lapply(1:G, function(g){
    make_first_entry_positive(Matrix::solve(data_inputs[['inner_prod_mat_response']]) %*% 
                                qr.Q(qr(matrix(rnorm(data_inputs[['nbasis_response']] * params[['pc_response']]),
                                               nrow = data_inputs[['nbasis_response']], ncol=params[['pc_response']]))))
  })
  parameters[['lambda_pcs_response']] = params[['lambda_pcs_response']][1]
  
  parameters[['Omega_predictors']] = lapply(1:G, function(g){
    make_first_entry_positive(Matrix::solve(data_inputs[['inner_prod_mat_predictors']]) %*%
                                qr.Q(qr(matrix(rnorm(sum(data_inputs[['nbasis_predictors']]) * params[['pc_predictors']]),
                                               nrow = sum(data_inputs[['nbasis_predictors']]), ncol=params[['pc_predictors']]))))
  })
  parameters[['lambda_pcs_predictors']] = params[['lambda_pcs_predictors']][1]
  
  # Use least squares estimate for the scores in order to obtain decent start values for the measurement error variances
  pcs_array_response = array(0, c(data_inputs[['n_profiles']], 1, params[['pc_response']]))
  pcs_array_predictors = array(0, c(data_inputs[['n_profiles_TS']], 1, params[['pc_predictors']]))
  parameters[['cluster_mat_BGC']] <- parameters[['cluster_mat']][data_inputs[['BGC']], , drop = F]
  
  for (i in 1:data_inputs[['n_profiles']]){
    cluster = parameters[['cluster_mat_BGC']][[i]]
    a = data_inputs[['response_profs']][[i]];
    b = data_inputs[['Phi_response_prof']][[i]] %*% parameters[['means_response']][[cluster]];
    c = data_inputs[['Phi_response_prof']][[i]] %*% parameters[['Omega_response']][[cluster]]
    pcs_array_response[i,1,] = tryCatch({as.double(Matrix::solve(Matrix::crossprod(c),Matrix::crossprod(c, a - b)))},
                                        error = function(x) {as.double(Matrix::solve(Matrix::crossprod(c) + diag(nrow = ncol(c), Matrix::diag(Matrix::crossprod(c))),Matrix::crossprod(c, a - b)))})
  }
  for (i in 1:data_inputs[['n_profiles_TS']]) {
    cluster = parameters[['cluster_mat']][[i]]
    a = data_inputs[['predictor_profs']][[i]];
    b = data_inputs[['Phi_predictors_prof']][[i]] %*% parameters[['means_predictors']][[cluster]];
    c = data_inputs[['Phi_predictors_prof']][[i]] %*% parameters[['Omega_predictors']][[cluster]]
    pcs_array_predictors[i,1,] = tryCatch({as.double(Matrix::solve(Matrix::crossprod(c),Matrix::crossprod(c, a - b)))},
                                          error = function(x) {as.double(Matrix::solve(Matrix::crossprod(c) + diag(nrow = ncol(c), Matrix::diag(Matrix::crossprod(c))),Matrix::crossprod(c, a - b)))})
  }
  
  parameters[['measurement_error_response']] = update_measurement_error(cluster_mat = parameters[['cluster_mat_BGC']],
                                                                        pcs_array1 = pcs_array_response,
                                                                        profs = data_inputs[['response_profs']],
                                                                        phi_prof = data_inputs[['Phi_response_prof']],
                                                                        means = parameters[['means_response']],
                                                                        Omegas1 = parameters[['Omega_response']],
                                                                        n_profiles = data_inputs[['n_profiles']],
                                                                        G=G,
                                                                        prof_ns = data_inputs[['profile_lengths_response']])
  
  parameters[['measurement_error_predictors']] = update_measurement_error(cluster_mat = parameters$cluster_mat,
                                                                          pcs_array1 = pcs_array_predictors,
                                                                          profs = data_inputs[['predictor_profs']],
                                                                          phi_prof = data_inputs[['Phi_predictors_prof']],
                                                                          means = parameters[['means_predictors']],
                                                                          Omegas1 = parameters[['Omega_predictors']],
                                                                          n_profiles = data_inputs$n_profiles_TS,
                                                                          prof_ns = data_inputs[['profile_lengths_predictors']], 
                                                                          G=G)
  
  parameters[['Lambdas']] = update_Lambda(parameters[['cluster_mat_BGC']], pcs_array_response,
                                          pcs_array_predictors[data_inputs[['BGC']],,,drop = F], params[['G']])
  
  orthogonalize_list = orthogonalize(Omegas1 = parameters[['Omega_response']], 
                                     inner_prod_mat1 = data_inputs[['inner_prod_mat_response']],
                                     G = G,
                                     Omegas2 = parameters[['Omega_predictors']],
                                     inner_prod_mat2 = data_inputs[['inner_prod_mat_predictors']],
                                     Lambdas = parameters[['Lambdas']],
                                     cluster_mat = parameters[['cluster_mat_BGC']],
                                     pcs1 = pcs_array_response, 
                                     pcs2 = pcs_array_predictors[data_inputs[['BGC']],,,drop = F])
  
  parameters[['Omega_response']] = orthogonalize_list[[1]]
  parameters[['Omega_predictors']]  = orthogonalize_list[[2]]
  parameters[['variances_response']] = orthogonalize_list[[3]]
  parameters[['variances_predictors']] = orthogonalize_list[[4]]
  parameters[['Lambdas']] = orthogonalize_list[[5]]
  parameters[['Sigma_eta_inv']] = orthogonalize_list[[6]]
  
  n_profiles = data_inputs[['n_profiles']]
  n_profiles_TS = data_inputs[['n_profiles_TS']]
  
  # Markov random field
  parameters[['theta']] <- .2
  parameters[['marginal_probabilities']] <- compute_marginal_probabilities(parameters[['cluster_mat']], parameters[['theta']], neighbors = data_inputs[['nn_list']], data_inputs[['n_profiles']], params[['G']])
  
  # conditional distributions
  parameters[['conditional_distributions']] = lapply(1:data_inputs$n_profiles_TS, function(x){
    lapply(1:params$G, function(g){
      return(list(NULL, NULL, NULL))
    })
  })
  
  # if(is.null(params[['lambda_mean_response']])){
  #   parameters[['lambda_mean_response']] = params[['lambda_mean_response']]
  # } else {
  #   parameters[['lambda_mean_response']] = 10
  # }
  # if (is.null(params[['lambda_mean_predictors']])){
  #   parameters[['lambda_mean_predictors']] = params[['lambda_mean_predictors']]
  # } else {
  #   parameters[['lambda_mean_predictors']] = 10
  # }
  # 
  return(parameters)
  
}

initialize.mult_var_spat <- function(model, params, data_inputs){
  
  var = params[['Y']]
  preds = params[['X']]
  
  id = params[['id']]
  ind_name = params[['ind_name']]
  
  domain1 = params[['domain_Y']]
  domains2 = params[['domain_X']]
  
  n_profiles = data_inputs[['n_profiles_TS']]
  G = params[['G']]
  
  # Initialize the linear model part
  class(model) = 'mult_var_ind'
  
  params[['pc_response']] = params[['pc_response1']]
  pc_response2 = params[['pc_response2']]
  params[['pc_response1']] = NULL
  params[['pc_response2']] = NULL
  
  parameters = initialize(model, params, data_inputs)
  
  for (i in 1:params[['EM_iter_init']]){
    MCEM_res = Estep(model, params, data_inputs, parameters)
    parameters = Mstep(model, params, data_inputs, parameters, MCEM_res, i)
  }
  
  
  MCEM_res = Estep(model, params, data_inputs, parameters)

  parameters[['Omega_response1']] <- lapply(1:G, function(g) {
    parameters$Omega_response[[g]] %*% parameters$Lambdas[[g]]
  })
  
  params[['pc_response2']] = pc_response2
  params[['pc_response1']] = params[['pc_response']]
  params[['pc_response']] = NULL
  
  if (params$pc_response2 > ncol(parameters$Omega_response[[1]])) {
    difference <- params$pc_response2 - ncol(parameters$Omega_response[[1]])
    parameters$Omega_response2 <- lapply(1:params[['G']], function(g) {
      cbind(parameters$Omega_response[[g]], matrix(nrow = nrow(parameters$Omega_response[[1]]),
                                                   ncol = difference,
                                                   rnorm(difference * nrow(parameters$Omega_response[[1]]))))
    })
  } else {
    parameters$Omega_response2 <- lapply(1:params[['G']], function(g) {
      matrix(rnorm(params$pc_response2 * nrow(parameters$Omega_response[[g]]), sd = 0.01), ncol= params$pc_response2)
      #parameters$Omega_response[[g]][, 1:params$pc_response2]
    })
  }
  
  parameters$Gammas <- lapply(1:params[['G']], function(g) {
    diag(params$pc_response)
  })
  
  parameters$variances_response <- lapply(1:params[['G']], function(g) {
    eigen(Matrix::solve(parameters$Sigma_eta_inv[[g]]))$values
  })
  
  parameters$variances_response <- lapply(1:params[['G']], function(g) {
    if (params[['pc_response2']] > length(parameters$variances_response[[g]])) {
      c(parameters$variances_response[[g]], rep(1000, params[['pc_response2']] - length(parameters$variances_response[[g]]) ))
    } else {
      parameters$variances_response[[g]][1:params[['pc_response2']]]
    }
  })
  if (is.null(params[['covariance_function']])) {
    n_range_params <- 2
  } else {
    n_range_params <- ifelse(params[['covariance_function']] == 'exponential_spheretime_warp',
                             7, 2)
  }
  vec_start <- c(.1, 50, 0,0,0,0, 0)
  
  parameters$range_params_p <- lapply(1:params[['G']], function(d) {
    matrix(nrow = n_range_params, ncol = params[['pc_response1']], vec_start[1:n_range_params])
  })
  
  parameters$range_params_r <- lapply(1:params[['G']], function(d) {
    matrix(nrow = n_range_params, ncol = params[['pc_response2']], vec_start[1:n_range_params])
  })
  
  parameters$conditional_distributions <- NULL
  
  # Preparation for running the single_var model on the functional residuals
  class(model) = 'mult_var_spat'
  
  MC_max = params[['MC_max']]
  params[['MC_max']] = 1
  
  MCEM_res = Estep(model, params, data_inputs, parameters)
  
  pcs_mat = get_relevant_scores(MCEM_res$cond_mean, params$pc_predictors, params$pc_response1, data_inputs$n_profiles_TS)
  
  data_inputs[['response_profs']] = centered_response(cluster_mat = as.matrix(parameters[['cluster_mat']][data_inputs[['BGC']],]),
                                                      profs = data_inputs$response_profs,
                                                      Phi_prof = data_inputs$Phi_response_prof,
                                                      Omegas = parameters[['Omega_response1']],
                                                      pcs = pcs_mat)
  
  # Run one iteration of the single var model to initialize the spatio-temporal pcs
  class(model) = 'one_var'
  
  params[['MC_max']] = MC_max
  params[['pc_response']] = params[['pc_response2']]
  params[['pc_response1']] = NULL
  params[['pc_response2']] = NULL
  params[['init_cluster']] = parameters[['cluster_mat_BGC']][,1, drop = F]
  
  if (!is.null(data_inputs[['locs']])) {
    data_inputs[['locs']] <- data_inputs[['locs']][data_inputs[['BGC']],]
  }
  if (!is.null(data_inputs[['time']])) {
    data_inputs[['time']] <- data_inputs[['time']][data_inputs[['BGC']]]
  }
  if (!is.null(data_inputs[['day']])) {
    data_inputs[['day']] <- data_inputs[['day']][data_inputs[['BGC']]]
  }
  data_inputs[['nn_list']] <- data_inputs[['nn_list']][data_inputs[['BGC']]]
  data_inputs[['nn_list']] <- 
    lapply(data_inputs[['nn_list']], function(x) {
      cumsum(data_inputs[['BGC']])[x[x %in% which(data_inputs[['BGC']])]]
    })

  parameters_temp = initialize(model, params, data_inputs)
  
  MCEM_res = Estep(model, params, data_inputs, parameters_temp)
  parameters_temp = Mstep(model, params, data_inputs, parameters_temp, MCEM_res, 1)
  
  parameters[['Omega_response2']] = parameters_temp[['Omega_response']]
  
  if(is.null(params[['lambda_lt']])){
    parameters[['lambda_lt']] = 1
  } else {
    parameters[['lambda_lt']] = params[['lambda_lt']][1]
  }
  if (!is.null(params[['lambda_pcs_predictors']])){
    parameters[['lambda_pcs_predictors']] = params[['lambda_pcs_predictors']][1]
  } else {
    parameters[['lambda_pcs_predictors']] = 10
  }
  # if(!is.null(params[['lambda_mean_response']])){
  #   parameters[['lambda_mean_response']] = params[['lambda_mean_response']]
  # } else {
  #   parameters[['lambda_mean_response']] = 10
  # }
  # if (!is.null(params[['lambda_mean_predictors']])){
  #   parameters[['lambda_mean_predictors']] = params[['lambda_mean_predictors']]
  # } else {
  #   parameters[['lambda_mean_predictors']] = 10
  # }
  
  return(parameters)
}


