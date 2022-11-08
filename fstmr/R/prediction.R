
##### Prediction set up functions #####
init_pred_params <- function(params, verbose) {
  if(is.null(params[['m_pred_predictors']]) | is.na(params[['m_pred_predictors']])) {
    if (verbose) {
      print('You have not provided params$m_pred_predictors, defaulting to 10')
    }
    params[['m_pred_predictors']] <- 10
  }
  if(is.null(params[['m_pred_response']]) | is.na(params[['m_pred_response']])) {
    if (verbose) {
      print('You have not provided params$m_pred_response, defaulting to 10')
    }
    params[['m_pred_response']] <- 10
  }
  if(is.null(params[['compute_variances']])) {
    if (verbose) {
      print('You have not specified whether to compute variances, defaulting yes')
    }
    params[['compute_variances']] <- T
  }
  params
}

# takes in df_list, sister function to prepare_EM_data
# df_list[[1]] contains the data with ind where you want to predict, if applicable
# df_list[[2]] contains any predictor data you may want to use, or additional locations without data
# include params and data_inputs to find nearest neighbors, etc.
prepare_pred_data <- function(model, df_list, params) {
  if (class(model) != 'mult_var_spat') {
    return(NULL)
  }
  if (class(df_list) == 'data.frame') {
    df_list <- list(df_list, df_list)
  }
  if (nrow(df_list[[1]]) == 0) {
    print('not planning prediction')
    return(NULL)
  }
  data_inputs <- model[['data_inputs']]
  var = params[['Y']]
  preds = params[['X']]
  id = params[['id']]
  ind_name = params[['ind_name']]
  loc_names = params[['loc_names']]
  time = params[['time']]
  doy = params[['dayofyear']]
  
  gdf_all = split(df_list[[2]], df_list[[2]][[id]])
  pred_inputs = list()
  pred_inputs[['has_data']] <- sapply(gdf_all, function(x) sum(!is.na(x[[preds[1]]])) > 0)
  gdf_new_pred <- gdf_all[!pred_inputs[['has_data']]]
  gdf <- gdf_all[pred_inputs[['has_data']]]
  pred_inputs[['n_predictions']] <- length(pred_inputs[['has_data']])
  
  if (!is.null(loc_names)) {
    pred_inputs[['locs']] <- t(sapply(gdf_all, function(x) as.double(x[,loc_names][1,])))
    colnames(pred_inputs[['locs']]) <- loc_names
  }
  if (!is.null(doy)){
    pred_inputs[['day']] <- sapply(gdf_all, function(x) as.double(x[c('dayofyear')][1,]))
  }
  if (!is.null(time)){
    pred_inputs[['time']] <- sapply(gdf_all, function(x) julian(as.Date(x[[time]][1],format = '%Y-%m-%d'), 
                                                                origin = as.Date('2000-01-01')))
  }
  
  pred_inputs[['neighbors']] <- get_nearest_neighbors(loc_names = loc_names,
                                                      locs1 = data_inputs[['locs']], 
                                                      locs2 = pred_inputs[['locs']],
                                                      time1 =  data_inputs[['time']],
                                                      time2 = pred_inputs[['time']],
                                                      day1 = data_inputs[['day']],
                                                      day2 = pred_inputs[['day']], 
                                                      nn = params[['nn']],
                                                      nn_range = params[['nn_range']],
                                                      nn_range_lat = params[['nn_range_lat']],
                                                      nn_range_lon = params[['nn_range_lon']],
                                                      nn_range_time = params[['nn_range_time']],
                                                      nn_type = params[['nn_type']],
                                                      nn_strategy = params[['nn_strategy']],
                                                      nn_lat_param = params[['nn_lat_param']],
                                                      nn_time_param = params[['nn_time_param']],
                                                      remove_between_land = params[['remove_between_land']])
  
  
  pred_inputs[['predictor_profs']] = lapply(gdf, function(x){
    if (!(preds[1] %in% colnames(x))) {
      return(NULL)
    }
    ret = c()
    for (i in 1:length(preds)){
      ret = c(ret, as.double(unlist(stats::na.omit(x[,preds[i]]))))
    }
    return(ret)
  })
  data_inputs[['predictor_profs']] = list() 
  data_inputs[['ind_predictors']] = list()
  counter = 1
  
  profile_lengths = matrix(0, nrow =length(gdf), ncol = length(preds));
  for(x in gdf){
    ret1 = c()
    ret2 = c()
    for (i in 1:length(preds)){
      t_df = x[!is.na(x[,preds[i]]),]
      ret1 = c(ret1, as.double(unlist(t_df[,preds[i]])))
      ret2 = c(ret2, as.double(unlist(t_df[,ind_name])))
      profile_lengths[counter, i] = nrow(t_df)
    }
    pred_inputs[['predictor_profs']][[counter]] = ret1
    pred_inputs[['ind_predictors']][[counter]] = ret2
    counter = counter + 1
  }
  pred_inputs[['profile_lengths_predictors']] = profile_lengths
  if (!is.null(pred_inputs[['ind_predictors']])){
    pred_inputs[['Phi_predictors_prof']] = get_phi_prof(indexes = pred_inputs[['ind_predictors']], 'predictors', 
                                                        params[['basis_predictors']],
                                                        prof_ns = pred_inputs[['profile_lengths_predictors']])
  } else {
    pred_inputs[['Phi_predictors_prof']] = list()
  }
  pred_inputs[['Phi_x_Phi_predictors']] <- lapply(pred_inputs[['Phi_predictors_prof']], function(x) as(Matrix::crossprod(x), 'dgCMatrix'))
  
  pred_inputs[['n_profiles_data']] = length(pred_inputs[['predictor_profs']])
  
  gdf_resp = split(df_list[[1]], df_list[[1]][[id]])
  pred_inputs[['has_data_resp']] <- sapply(gdf_resp, function(x) sum(!is.na(x[[var[1]]])) > 0)
  gdf_resp <- gdf_resp[pred_inputs[['has_data_resp']]]
  
  pred_inputs[['response_profs']] = lapply(gdf_resp, function(x) as.double(unlist(x[!is.na(x[,var]), var])))
  pred_inputs[['ind_response']] = lapply(gdf_resp, function(x) as.double(unlist(x[!is.na(x[,var]), ind_name])))
  
  pred_inputs[['Phi_response_prof']] =  get_phi_prof(indexes = pred_inputs[['ind_response']], 
                                                     kind = 'response', 
                                                     basis_list = params[['basis_response']])
  
  pred_inputs
}

set_up_prediction <- function(estimated_model, params,
                              core_data = NULL, grid = NULL) {
  
  params <- init_pred_params(params, verbose = T)
  G <- params[['G']]
  MC_max <- params[['MC_max']]
  data_inputs <- estimated_model[["data_inputs"]]
  parameters <- estimated_model[["parameters"]]
  
  training_in_pred <- !data_inputs[['BGC']]
  
  # compute marginal probabilities of clusters
  # if we are predicting using core data and a grid
  if (!is.null(core_data) & !is.null(grid)) {
    
    # deal with core data
    pred_inputs_core <- prepare_pred_data(estimated_model, list(core_data, core_data), params)
    
    # predict core memberships
    marginal_probs_core <- compute_marginal_probabilities(parameters[['cluster_membership']], 
                                                          pred_inputs_core[['neighbors']],
                                                          parameters[['theta']],
                                                          pred_inputs_core[['n_predictions']],
                                                          G = params[['G']])
    memberships_core <- predict_group_membership_from_preds(estimated_model,
                                                            pred_inputs_core, 
                                                            marginal_probs_core)[[1]]
    
    # predict grid memberships for initialization
    neighbors_grid <- get_nearest_neighbors(loc_names = c('longitude', 'latitude'), 
                                            locs1 = rbind(data_inputs[["locs"]], pred_inputs_core[['locs']]), locs2 = grid[,c('longitude', 'latitude')], 
                                            time1 = c(data_inputs[["time"]], pred_inputs_core[['time']]), time2 = grid[['dayofyear']], 
                                            day1 = c(data_inputs[["day"]],pred_inputs_core[['day']]), day2 = grid[['dayofyear']], 
                                            nn = params[["nn"]], nn_range = params[["nn_range"]], 
                                            nn_range_lat = params[["nn_range_lat"]], nn_range_lon = params[["nn_range_lon"]], 
                                            nn_range_time = params[["nn_range_time"]], nn_type = params[["nn_type"]], 
                                            nn_strategy = params[["nn_strategy"]], nn_lat_param = params[["nn_lat_param"]], 
                                            nn_time_param = params[["nn_time_param"]], remove_between_land = params[["remove_between_land"]])
    
    marginal_probs_grid <- compute_marginal_probabilities(c(parameters[['cluster_membership']], memberships_core), 
                                                          neighbors_grid,
                                                          parameters[['theta']],
                                                          length(neighbors_grid),
                                                          G = params[['G']])
    memberships <- apply(marginal_probs_grid, 1, which.max)
    
    locations_all <- rbind(data_inputs[['locs']], pred_inputs_core[['locs']], grid[,c('longitude', 'latitude')])
    day_all <- c(data_inputs[['day']], pred_inputs_core[['day']], grid[['day']])
    
    probs <- return_probabilities(vars_r = parameters[['variances_response']],
                                  response_profs = c(data_inputs[['response_profs']]),
                                  is_bgc = c(rep(T, length(data_inputs[['response_profs']])),
                                             rep(F, length(pred_inputs_core[['predictor_profs']]))),
                                  phi_resp = data_inputs[['Phi_response_prof']],
                                  means_resp = parameters[['means_response']],
                                  me_resp = parameters[['measurement_error_response']],
                                  Omegas_r1 = parameters[['Omega_response1']],
                                  Omegas_r2 = parameters[['Omega_response2']],
                                  means_pred = parameters[['means_predictors']],
                                  me_pred = parameters[['measurement_error_predictors']],
                                  vars_p = parameters[['variances_predictors']],
                                  predictor_profs = c(data_inputs[['predictor_profs']],
                                                      pred_inputs_core[['predictor_profs']]),
                                  phi_pred = c(data_inputs[['Phi_predictors_prof']],
                                               pred_inputs_core[['Phi_predictors_prof']]),
                                  Omegas_pred = parameters[['Omega_predictors']],
                                  profile_lengths_p = rbind(data_inputs[['profile_lengths_predictors']],
                                                            pred_inputs_core[['profile_lengths_predictors']]),
                                  G = G,
                                  n_profiles = data_inputs[['n_profiles_TS']] + pred_inputs_core[['n_predictions']], 
                                  n_preds = 2)
    has_data <- c(rep(T, data_inputs[['n_profiles_TS']] + pred_inputs_core[['n_predictions']]), 
                  rep(F, nrow(grid)))
    
    prediction_model <- estimated_model
    prediction_model[['data_inputs']][['BGC']] <- c(data_inputs[['BGC']],
                                                    rep(F, length(pred_inputs_core[['predictor_profs']])))
    prediction_model[['data_inputs']][['Phi_response_prof']] <- data_inputs[['Phi_response_prof']]
    prediction_model[['data_inputs']][['Phi_predictors_prof']] <- c(data_inputs[['Phi_predictors_prof']],
                                                                    pred_inputs_core[['Phi_predictors_prof']])
    prediction_model[['data_inputs']][['response_profs']] <- data_inputs[['response_profs']]
    prediction_model[['data_inputs']][['predictor_profs']] <- c(data_inputs[['predictor_profs']], 
                                                                pred_inputs_core[['predictor_profs']])
    prediction_model[['data_inputs']][['Phi_x_Phi_response']] <- data_inputs[['Phi_x_Phi_response']]
    prediction_model[['data_inputs']][['Phi_x_Phi_predictors']] <- c(data_inputs[['Phi_x_Phi_predictors']],
                                                                     pred_inputs_core[['Phi_x_Phi_predictors']])
    prediction_model[['data_inputs']][['nbasis_predictors']] <- data_inputs[['nbasis_predictors']]
    prediction_model[['data_inputs']][['n_profiles']] <- data_inputs[['n_profiles']]
    prediction_model[['data_inputs']][['n_profiles_TS']] <- data_inputs[['n_profiles_TS']] + 
      pred_inputs_core[['n_predictions']]
    prediction_model[['data_inputs']][['profile_lengths_predictors']] <- rbind(data_inputs[['profile_lengths_predictors']],
                                                                               pred_inputs_core[['profile_lengths_predictors']])
    prediction_model[['data_inputs']][['locs']] <- rbind(data_inputs[['locs']],
                                                         pred_inputs_core[['locs']],
                                                         grid[,c('longitude', 'latitude')])
    
    prediction_model[['data_inputs']][['time']] <- c(data_inputs[['time']], pred_inputs_core[['time']], 
                                                     grid[['time']])
    time_used <- c(data_inputs[['time']], pred_inputs_core[['time']], 
                   grid[['time']])
    
    is_prediction <- which(c(training_in_pred,
                             rep(F, length(pred_inputs_core[['time']])),
                             rep(T, length(grid[['time']]))))
  } else { # in a leave out setting with no additional data
    memberships <- NULL
    locations_all <- data_inputs[['locs']]
    day_all <- data_inputs[['day']]
    
    probs <- return_probabilities(vars_r = parameters[['variances_response']],
                                  response_profs = data_inputs[['response_profs']],
                                  is_bgc = data_inputs[['BGC']],
                                  phi_resp = c(data_inputs[['Phi_response_prof']]),
                                  means_resp = parameters[['means_response']],
                                  me_resp = parameters[['measurement_error_response']],
                                  Omegas_r1 = parameters[['Omega_response1']],
                                  Omegas_r2 = parameters[['Omega_response2']],
                                  means_pred = parameters[['means_predictors']],
                                  me_pred = parameters[['measurement_error_predictors']],
                                  vars_p = parameters[['variances_predictors']],
                                  predictor_profs = data_inputs[['predictor_profs']],
                                  phi_pred = data_inputs[['Phi_predictors_prof']],
                                  Omegas_pred = parameters[['Omega_predictors']],
                                  profile_lengths_p = data_inputs[['profile_lengths_predictors']],
                                  G = G,
                                  n_profiles = data_inputs[['n_profiles_TS']], 
                                  n_preds = 2)
    is_prediction <- which(training_in_pred)
    has_data <- c(rep(T, data_inputs[['n_profiles_TS']]))
    time_used <- data_inputs[['time']]
    prediction_model <- estimated_model
  }
  
  neighbors_all <- get_nearest_neighbors(loc_names = c('longitude', 'latitude'), 
                                         locs1 = locations_all, 
                                         day1 = day_all,
                                         nn = params[["nn"]], nn_range = params[["nn_range"]], 
                                         nn_range_lat = params[["nn_range_lat"]], nn_range_lon = params[["nn_range_lon"]], 
                                         nn_range_time = params[["nn_range_time"]], nn_type = params[["nn_type"]], 
                                         nn_strategy = params[["nn_strategy"]], nn_lat_param = params[["nn_lat_param"]], 
                                         nn_time_param = params[["nn_time_param"]], remove_between_land = params[["remove_between_land"]])
  
  p_mat <- probs[[1]]; like_mat <- probs[[2]]
  
  list('p_mat' = p_mat, 'like_mat' = like_mat,
       'neighbors_all' = neighbors_all, 
       'memberships' = memberships, 
       'has_data' = has_data, 
       'prediction_model' = prediction_model,
       'is_prediction' = is_prediction,
       'time' = time_used)
}




# based on predictor information in pred_inputs
# compute conditional probabilities and memberships
predict_group_membership_from_preds <- function(model, pred_inputs, marginal_probs){
  parameters <- model[['parameters']]
  n_profiles = length(pred_inputs[['predictor_profs']])
  n_predictions = pred_inputs[['n_predictions']]
  n_groups = length(parameters[['means_predictors']])
  memberships = rep(NA, n_predictions)
  likelihoods = rep(NA, n_groups)
  x_data <- 0
  cond_prob <- matrix(nrow = n_predictions, ncol = n_groups)
  for(x in 1:n_predictions){
    if (pred_inputs[['has_data']][x] == F) {
      memberships[x] = which.max(marginal_probs[x,])
      cond_prob[x,] <- marginal_probs[x,]
      next
    }
    x_data <- x_data + 1
    measurement_errors_vector = rep(parameters[['measurement_error_predictors']],  pred_inputs[['profile_lengths_predictors']][x_data,])
    for(g in 1:n_groups){
      phi_lambda = pred_inputs[['Phi_predictors_prof']][[x_data]] %*% parameters[['Omega_predictors']][[g]]
      mean = pred_inputs[['Phi_predictors_prof']][[x_data]] %*% parameters[['means_predictors']][[g]]
      likelihoods[g] = 
        compute_likelihood(pred_inputs[['predictor_profs']][[x_data]], mean, measurement_errors_vector, phi_lambda, parameters[['variances_predictors']][[g]])
    }
    cond_prob[x,] <- log_like_to_probs(likelihoods, marg = marginal_probs[x,], G = n_groups)
    memberships[x] = which.max(cond_prob[x,])
  }
  return(list(memberships, cond_prob))
}

compute_likelihood <- function(x, mean, sigma, U, W){
  ret = c_lik_eigen_sherman_pred(as.numeric(x), as.numeric(mean), as.numeric(sigma), as.matrix(U), as.numeric(W))
  return(ret)
}



rPotts_prediction <- function(G, neighbors, has_data, p_mat, likelihoods,
                              theta, n_samples = 1, skip = 3, init_pred_cluster = NULL){
  ngibbs = n_samples*skip+skip
  nobs = nrow(p_mat)
  ntotal = length(neighbors)
  
  samples = matrix(0, nrow = ntotal, ncol = n_samples)
  potts_sample = rep(0, ntotal)
  
  if (is.null(init_pred_cluster)) {
    potts_sample[!has_data] = sample(1:G, size = ntotal-nobs, replace = T)
  } else {
    potts_sample[!has_data] = init_pred_cluster
  }
  data_indexes <- which(has_data)
  for (i in 1:nobs) {
    potts_sample[data_indexes[i]] = sample(1:G, size = 1, prob = p_mat[i,], replace = T)
  }
  cluster_mat <- matrix(0, length(potts_sample), G)
  cluster_mat[cbind(1:length(potts_sample), potts_sample)] <- 1
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
  
  uniform = -log(rep(G, G))
  counter = 1
  
  for(i in 1:ngibbs){
    dc = 1
    for(x in 1:ntotal){
      if (has_data[x]){
        prob = log_like_to_probs(likelihoods[dc,], G, marg = marginal_membership[x,])
        dc = dc+1
      } else {
        prob = log_like_to_probs(uniform, G, marg = marginal_membership[x,])
      }
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
  samples
}

pred_E_step <- function(pred_set_up, params, cluster_mat, prediction_only = F) {
  n_total <- length(pred_set_up[['neighbors_all']])
  
  pcs_array_predictors <- array(0, dim = c(n_total, params[['MC_max']], params[['pc_predictors']]))
  pcs_array_response <- array(0, dim = c(n_total, params[['MC_max']], params[['pc_response2']]))
  
  data_inputs <- pred_set_up[['prediction_model']][['data_inputs']]
  parameters <- pred_set_up[['prediction_model']][['parameters']]
  
  V_val_list <- list()
  # save vecchia information so that we don't have to redo it in update_variances/range_params
  E_step_list <- list()
  like_spatial <- rep(0, params[['MC_max']])
  like_ind <- rep(0, params[['MC_max']])
  start <- (n_total)*params[['pc_predictors']]
  
  for (mc in 1:params[['MC_max']]){
    print(mc)
    print(Sys.time())
    cmi = cluster_mat[,mc] -1
    mode(cmi) = 'integer'
    
    if (mc == 1 || !all(cluster_mat[,mc] == cluster_mat[,mc - 1])){
      if (mc > 1 & length(data_inputs[['Phi_x_Phi_predictors']]) > 50000) {
        rm(E_step_list)
        E_step_list <- list()
        gc()
      }
      
      E_step_list[['UTU']] = c_compute_UTU(phi_x_phi_r = data_inputs[['Phi_x_Phi_response']],
                                           phi_x_phi_p = data_inputs[['Phi_x_Phi_predictors']],
                                           Omegas_r1 = parameters[['Omega_response1']],
                                           Omegas_r2 = parameters[['Omega_response2']],
                                           Omegas_p = parameters[['Omega_predictors']],
                                           me_r = as.numeric(parameters[['measurement_error_response']]),
                                           me_p = as.numeric(parameters[['measurement_error_predictors']]),
                                           n_basis_p = as.integer(data_inputs[['nbasis_predictors']]), 
                                           cmi[pred_set_up[['has_data']]], 
                                           data_inputs[['n_profiles']],
                                           data_inputs[['n_profiles_TS']],
                                           params[['G']],
                                           as.integer(data_inputs[['BGC']]))
      
      
      means_mat_r = as.matrix(do.call(cbind, parameters[['means_response']]))
      means_mat_p = as.matrix(do.call(cbind, parameters[['means_predictors']]))
      UTX <- c_compute_UTX(data_inputs[['Phi_response_prof']], data_inputs[['Phi_predictors_prof']],
                           data_inputs[['response_profs']], data_inputs[['predictor_profs']],
                           parameters[['Omega_response1']], parameters[['Omega_response2']], parameters[['Omega_predictors']],
                           means_mat_r, means_mat_p,
                           cmi[pred_set_up[['has_data']]], as.numeric(parameters[['measurement_error_response']]),
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
                                         cmi[pred_set_up[['has_data']]], 
                                         as.integer(data_inputs[['nbasis_predictors']]), 
                                         as.integer(data_inputs[['BGC']]))
      
      # Create the spatial covariance matrices
      cluster_mat_BGC <- cluster_mat[pred_set_up[['has_data']],][data_inputs[['BGC']],mc]
      E_step_list[['UTU']] = UTU = Matrix::forceSymmetric(E_step_list[['UTU']], uplo = 'L')
      
      if (is.null(data_inputs[['time']])){
        locs = as.matrix(data_inputs[['locs']])
        locs_data = locs[pred_set_up[['has_data']],]
      } else {
        locs = as.matrix(cbind(data_inputs[['locs']], time = pred_set_up[['time']]))
        locs_data = locs[pred_set_up[['has_data']],]
      }
      
      pred_scores_chol <- spatial_inv_covariance_scores(z = cluster_mat[,mc], n_pcs = params[['pc_predictors']], 
                                                        G = params[['G']], variances = parameters[['variances_predictors']],
                                                        range_params = parameters[['range_params_p']],
                                                        model=ifelse(!is.null(params[['covariance_function']]), 'vecchia', 'independence'),                                                        covfun_name = params[['covariance_function']], 
                                                        locs = locs,
                                                        m = params[['m_pred_predictors']])
      
      resp_scores_chol <- spatial_inv_covariance_scores(cluster_mat[, mc], n_pcs = params[['pc_response2']],
                                                        G = params[['G']], variances = parameters[['variances_response']],
                                                        range_params = parameters[['range_params_r']],
                                                        model=ifelse(!is.null(params[['covariance_function']]), 'vecchia', 'independence'),                                                        covfun_name = params[['covariance_function']],
                                                        locs = locs,
                                                        m = params[['m_pred_response']])
      
      # reorder things back to the way it was before
      pred_scores_full <- Matrix::crossprod(pred_scores_chol[['Linv_all']])
      if (!is.null(params[['covariance_function']])) {
        reorder_pred <- order(pred_scores_chol[['vecchia_info']][['order_overall']])
      } else {
        reorder_pred <- 1:nrow(pred_scores_full)
      }
      pred_scores_full_reordered <- pred_scores_full[reorder_pred, reorder_pred]
      
      resp_scores_full <- Matrix::crossprod(resp_scores_chol[['Linv_all']])
      if (!is.null(params[['covariance_function']])) {
        reorder_resp <- order(resp_scores_chol[['vecchia_info']][['order_overall']])
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
      
      n_total_pred <- nrow(pred_scores_full_reordered)
      n_scores_data_pred <- data_inputs[['n_profiles_TS']] * params[['pc_predictors']]
      
      fill_vec <- n_total_pred + 
        which(rep(data_inputs[['BGC']], each = params[['pc_response2']]))
      full_fill <- c(1:n_scores_data_pred, fill_vec)
      
      UTU_filled <- create_expanded_UTU(UTU, full_fill, nrow(E_step_list[['A_inv']]))
      UTX_filled <- rep(0, nrow(E_step_list[['A_inv']]))
      UTX_filled[full_fill] <- UTX
      
      # Form joint score matrix of training and predictions
      E_step_list[['Sigma_inv']] <- E_step_list[['A_inv']] + UTU_filled
      E_step_list[['Sigma_inv_chol']]  <- 
        Matrix::Cholesky(E_step_list[['Sigma_inv']], perm = T, super = NA, LDL = F)
      
      cond_mean = Matrix::solve(E_step_list[['Sigma_inv_chol']], UTX_filled, system = 'A')
      
      # compute likelihoods
      if (mc == 1 || mean(cluster_mat[pred_set_up[['has_data']],mc] != 
                          cluster_mat[pred_set_up[['has_data']],mc-1]) > 0.001) {
        ######## DATA ONLY
        pred_chol_data <- spatial_inv_covariance_scores(z = cluster_mat[pred_set_up[['has_data']],mc], 
                                                        n_pcs = params[['pc_predictors']], 
                                                        G = params[['G']], variances = parameters[['variances_predictors']],
                                                        range_params = parameters[['range_params_p']],
                                                        model=ifelse(!is.null(params[['covariance_function']]), 'vecchia', 'independence'),
                                                        covfun_name = params[['covariance_function']], 
                                                        locs = locs_data,
                                                        m = params[['m_pred_predictors']])
        
        resp_chol_data <- spatial_inv_covariance_scores(cluster_mat[pred_set_up[['has_data']], mc][data_inputs[['BGC']]], n_pcs = params[['pc_response2']],
                                                        G = params[['G']], variances = parameters[['variances_response']],
                                                        range_params = parameters[['range_params_r']],
                                                        model=ifelse(!is.null(params[['covariance_function']]), 'vecchia', 'independence'),
                                                        covfun_name = params[['covariance_function']],
                                                        locs = locs_data[data_inputs[['BGC']],],
                                                        m = params[['m_pred_response']])
        
        # reorder things back to the way it was before
        pred_scores_data <- Matrix::crossprod(pred_chol_data[['Linv_all']])
        if (!is.null(params[['covariance_function']])) {
          reorder_pred <- order(pred_chol_data[['vecchia_info']][['order_overall']])
        } else {
          reorder_pred <- 1:nrow(pred_scores_data)
        }
        pred_scores_data <- pred_scores_data[reorder_pred, reorder_pred]
        
        resp_scores_data <- Matrix::crossprod(resp_chol_data[['Linv_all']])
        if (!is.null(params[['covariance_function']])) {
          reorder_resp <- order(resp_chol_data[['vecchia_info']][['order_overall']])
        } else {
          reorder_resp <- 1:nrow(resp_scores_data)
        }
        resp_scores_data <- resp_scores_data[reorder_resp, reorder_resp]
        
        
        # create this matrix, we will need its determinant, but it is not reordered
        E_step_list[['A_inv_chol_data_wrong_order']] <- Matrix::bdiag(pred_chol_data[['Linv_all']], 
                                                                      resp_chol_data[['Linv_all']])
        
        E_step_list[['A_inv_data']] <- Matrix::bdiag(pred_scores_data,
                                                     resp_scores_data)
        E_step_list[['Sigma_inv_data']] <- E_step_list[['A_inv_data']] + UTU
        E_step_list[['Sigma_inv_data_chol']]  <-
          Matrix::Cholesky(E_step_list[['Sigma_inv_data']], perm = T, super = NA, LDL = F)
        cond_mean_data = Matrix::solve(E_step_list[['Sigma_inv_data_chol']], UTX, system = 'A')
        
        data_ll <- compute_likelihood_for_IS(sigma_chol = E_step_list[['Sigma_inv_data_chol']],
                                             A_inv_chol = E_step_list[['A_inv_chol_data_wrong_order']], 
                                             me_response = parameters[['measurement_error_response']],
                                             response_profs = data_inputs[['response_profs']],
                                             cond_mean = cond_mean_data,
                                             UTX = UTX,
                                             sse_mean = sse_mean,
                                             predictor_profs = data_inputs[['predictor_profs']],
                                             profile_lengths_predictors = data_inputs[['profile_lengths_predictors']],
                                             me_predictors = parameters[['measurement_error_predictors']])
        normalizing_sums <- apply(pred_set_up[['like_mat']], 1, function(x) max(x) + log(sum(exp(x - max(x)))))
        ind_ll <-sum(pred_set_up[['like_mat']][cbind(1:nrow(cluster_mat[pred_set_up[['has_data']],]), 
                                                     cluster_mat[pred_set_up[['has_data']],mc])])- 
          sum(normalizing_sums)
      }
      # compute variances
      V_val_list[[mc]] <- list()
      Sigma_sparse_band <- invert_selected_entries(pred_set_up[['is_prediction']],
                                                   params, cmi, E_step_list[['Sigma_inv']])
      # only take blocks for prediction locations to save memory space
      for (i in 1:length(pred_set_up[['is_prediction']])) { 
        predictor_score_index <- ( (i-1) * (params[['pc_predictors']]) + 1):( i * (params[['pc_predictors']]))
        response_score_index <-length(pred_set_up[['is_prediction']]) * params[['pc_predictors']] +
          ( (i-1) * (params[['pc_response2']]) + 1):(i * (params[['pc_response2']]))
        indexes <- c(predictor_score_index, response_score_index)
        V_val_list[[mc]][[i]] <- Sigma_sparse_band[indexes,indexes]
      }
      
    } else {
      # compute variances
      V_val_list[[mc]] <- V_val_list[[mc-1]]
    }
    like_ind[mc] <- ind_ll
    like_spatial[mc] <- data_ll
    
    pcs_array_predictors[,mc,] = matrix(cond_mean[1:start], nrow = n_total, byrow = T)
    pcs_array_response[,mc,] = matrix(cond_mean[-c(1:start)], nrow = n_total, byrow = T)
  }
  weights <- like_spatial - like_ind
  weights_normal <- vector('numeric', length(weights))
  for (i in 1:length(weights)) {
    weights_normal[i] <- 1/(sum(exp(weights - weights[i])))
  }
  wt <- weights_normal * params[['MC_max']]
  
  if (prediction_only) {
    pcs_array_predictors <- pcs_array_predictors[pred_set_up[['is_prediction']],,]
    pcs_array_response <- pcs_array_response[pred_set_up[['is_prediction']],,]
  }
  
  list('wt' = wt, 'pcs_array_predictors' = pcs_array_predictors,
       'pcs_array_response' = pcs_array_response, 'like_ind' = like_ind,
       'like_spatial' = like_spatial,
       'variances' = V_val_list)
}

pressure_level_prediction <- function(held_out_profs, params, pred_set_up, pred_score_results,
                                      cluster_mat) {
  gdf_resp = split(held_out_profs, held_out_profs[[params[['id']]]])
  ind_response_pred <- lapply(gdf_resp, function(x) 
    as.double(unlist(x[!is.na(x[,params[['var']]]), params[['ind_name']]])))
  response_pred <- lapply(gdf_resp, function(x) 
    as.double(unlist(x[!is.na(x[,params[['var']]]), params[['var']]])))
  Phi_response_prof <- get_phi_prof(indexes = ind_response_pred, 
                                    kind = 'response', 
                                    basis_list = params[['basis_response']])
  prof_names <- names(ind_response_pred)
  prof_lengths <- sapply(ind_response_pred,length)
  pred_list <- list()
  for(mc in 1:params[['MC_max']]) {
    pred_full <- predict_at_pressure_levels(pred_set_up[['prediction_model']], 
                                            membership = cluster_mat[pred_set_up[['is_prediction']],mc],
                                            E_alpha_mat = pred_score_results[['pcs_array_predictors']][pred_set_up[['is_prediction']],mc,],
                                            E_eta_mat =  pred_score_results[['pcs_array_response']][pred_set_up[['is_prediction']],mc,],
                                            Phi_response_prof = Phi_response_prof,
                                            prediction_part = "full")
    pred_TS <- predict_at_pressure_levels(pred_set_up[['prediction_model']], 
                                          membership = cluster_mat[pred_set_up[['is_prediction']],mc],
                                          E_alpha_mat = pred_score_results[['pcs_array_predictors']][pred_set_up[['is_prediction']],mc,],
                                          E_eta_mat = pred_score_results[['pcs_array_response']][pred_set_up[['is_prediction']],mc,],
                                          Phi_response_prof = Phi_response_prof,
                                          prediction_part = "TS")
    pred_mean <- predict_at_pressure_levels(pred_set_up[['prediction_model']], 
                                            membership = cluster_mat[pred_set_up[['is_prediction']],mc],
                                            E_alpha_mat = pred_score_results[['pcs_array_predictors']][pred_set_up[['is_prediction']],mc,],
                                            E_eta_mat = pred_score_results[['pcs_array_response']][pred_set_up[['is_prediction']],mc,],
                                            Phi_response_prof = Phi_response_prof,
                                            prediction_part = "mean")
    var_at_pressure <- var_at_pressure_levels(pred_set_up[['prediction_model']], var_spatial_scores = pred_score_results[['variances']][[mc]],
                                              membership = cluster_mat[pred_set_up[['is_prediction']],mc], 
                                              Phi_response_prof = Phi_response_prof)
    pred_list[[mc]] <- data.frame(pred_full, pred_TS, pred_mean, 
                                  var_at_pressure= unlist(var_at_pressure),
                                  mc= mc,
                                  profile_unique = rep(prof_names, prof_lengths),
                                  index = 1:length(pred_full),
                                  weight = pred_score_results[['wt']][mc])
  }
  data_info <- data.frame('Y' = unlist(response_pred),
                          'ind' = unlist(ind_response_pred), 
                          'index' = 1:length(pred_full),
                          'profile_unique' = rep(prof_names, prof_lengths))
  pred_df<- dplyr::bind_rows(pred_list)
  mean_prediction <- pred_df %>%
    group_by(index) %>%
    summarise(avg_full = mean(pred_full*weight),
              avg_TS = mean(pred_TS*weight),
              avg_mean = mean(pred_mean*weight),
              avg_var = mean(var_at_pressure*weight),
              var_avg = mean((pred_full - mean(pred_full*weight))^2 * weight),
              profile_unique = profile_unique[1])
  list('pred_df' = pred_df, 'mean_prediction' = mean_prediction,
       'data_info' = data_info)
}


# utility function for prediction
create_expanded_UTU <- function(UTU, full_fill, nrow_A) {
  # create a bigger UTU that includes predictions
  UTU_dgT <- as(UTU, 'dgTMatrix')
  UTU_filled <- as(Matrix::Matrix(0, nrow_A, nrow_A, doDiag = F), 'dgTMatrix')
  UTU_filled@i <- as.integer(full_fill[UTU_dgT@i+1]-1)
  UTU_filled@j <- as.integer(full_fill[UTU_dgT@j+1]-1)
  UTU_filled@x <- UTU_dgT@x
  Matrix::forceSymmetric(UTU_filled, uplo = 'L')
}

# utility function for prediction variances
invert_selected_entries <- function(is_prediction, params, cmi_all, Sigma_inv) {
  predictor_score_index <- as.vector(sapply(is_prediction, function(x) ( (x-1) * (params[['pc_predictors']]) + 1):( x * (params[['pc_predictors']]))))
  response_score_index <- as.vector(sapply(is_prediction, function(x) {length(cmi_all) * params[['pc_predictors']] +
      ( (x-1) * (params[['pc_response2']]) + 1):( x * (params[['pc_response2']]))}))
  indexes <- c(predictor_score_index, response_score_index)
  Sigma_inv_subset <- Sigma_inv[indexes, indexes]
  Sigma_inv_subset_chol <- Matrix::Cholesky(Sigma_inv_subset, perm = T, super = NA, LDL = F)
  Sigma_sparseinv <- sparseinv::Takahashi_Davis(Q = Sigma_inv_subset,
                                                cholQp = as(Sigma_inv_subset_chol, 'Matrix'), 
                                                P = Matrix::t(as(Sigma_inv_subset_chol, 'pMatrix')))
  suppressWarnings({mat_return <- Matrix::band(Sigma_sparseinv, 
                                               k1 =-(params[['pc_predictors']] + params[['pc_response2']]),
                                               k2 =params[['pc_predictors']] + params[['pc_response2']] )})
  mat_return
}

# takes predict scores, predicts at pressures defined by pred_inputs[['Phi_response_prof']]
predict_at_pressure_levels <- function(model, membership,Phi_response_prof,
                                       prediction_part = 'full', 
                                       E_alpha_mat, E_eta_mat){
  if (class(model) != 'mult_var_spat') {
    return(NULL)
  }
  parameters <- model[['parameters']]
  
  G <- length(parameters[['means_response']])
  predictions = lapply(1:length(Phi_response_prof), function(x){
    mean_vals <- as.double(Phi_response_prof[[x]] %*% parameters[['means_response']][[membership[x]]])
    TS_vals <- as.double(Phi_response_prof[[x]] %*% (parameters[['Omega_response1']][[membership[x]]] 
                                                     %*% E_alpha_mat[x,]))
    oxy_vals <- as.double(Phi_response_prof[[x]] %*% (parameters[['Omega_response2']][[membership[x]]] 
                                                      %*% E_eta_mat[x,]))
    if (prediction_part == 'full') {
      mean_vals + TS_vals + oxy_vals
    } else if (prediction_part == 'mean') {
      mean_vals
    } else if (prediction_part == 'TS') {
      mean_vals + TS_vals
    }
  })
  as.double(unlist(sapply(predictions, as.double)))
}

# takes predicted scores/variances, predicts variance 
# at pressures defined by pred_inputs[['Phi_response_prof']]
var_at_pressure_levels <- function(model, 
                                   var_spatial_scores,Phi_response_prof,
                                   membership){
  if (class(model) != 'mult_var_spat') {
    return(NULL)
  }
  parameters <- model[['parameters']]
  if (is.null(var_spatial_scores) | length(var_spatial_scores) < 1) {
    print('You have not supplied variance information of the scores')
    return(NA)
  }
  G <- length(parameters[['means_response']])
  
  predictions = lapply(1:length(Phi_response_prof), function(x){
    Phi_all <- cbind(Phi_response_prof[[x]] %*% parameters[['Omega_response1']][[membership[x]]],
                     Phi_response_prof[[x]] %*% parameters[['Omega_response2']][[membership[x]]])
    
    Phi_var <- Phi_all %*% var_spatial_scores[[x]]
    Matrix::rowSums(Phi_var * Phi_all) + parameters[['measurement_error_response']]
  })
}

