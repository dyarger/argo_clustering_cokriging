
# for one set of predicted scores, compute prediction at all locations at p
# type can denote temperature, salinity, or parts of the oxygen predictions
predict_at_pressure <- function(p, parameters, pred, params, type = 'full') {
  Phi <- as(eval.basis(basisobj = params[['basis_response']], p), 'sparseMatrix')
  Phi_pred <- bdiag(lapply(1:length(params[['basis_predictors']]), function(x) {
    as(eval.basis(basisobj = params[['basis_predictors']][[x]], p), 'sparseMatrix')
  }))
  mean_predictor_g <- lapply(1:params[['G']], function(g) {
    as.double(Phi_pred %*% parameters[['means_predictors']][[g]])
  })
  pc_predictor_g <- lapply(1:params[['G']], function(g) {
    as.matrix(Phi_pred %*% parameters[['Omega_predictors']][[g]])
  })
  mean_response_g <- lapply(1:params[['G']], function(g) {
    as.double(Phi %*% parameters[['means_response']][[g]])
  })
  pc_response1_g <- lapply(1:params[['G']], function(g) {
    as.matrix(Phi %*% parameters[['Omega_response1']][[g]])
  })
  pc_response2_g <- lapply(1:params[['G']], function(g) {
    as.matrix(Phi %*% parameters[['Omega_response2']][[g]])
  })
  pred_at_p <- rep(NA, nrow(pred))
  E_alpha_mat <- as.matrix(pred[,paste0('E_alpha.', 1:ncol(parameters[['Omega_response1']][[1]]))])
  E_eta_mat <- as.matrix(pred[,paste0('E_eta.', 1:ncol(parameters[['Omega_response2']][[1]]))])
  memberships <- pred$membership
  for (i in 1:nrow(pred)) {
    pred_scores <- E_alpha_mat[i,]
    resp_scores <- E_eta_mat[i,]
    if (type == 'temperature') {
      pred_at_p[i] <- mean_predictor_g[[memberships[i]]][1] + 
        (pc_predictor_g[[memberships[i]]] %*% pred_scores)[1,1]
    } else if(type == 'salinity') {
      pred_at_p[i] <- mean_predictor_g[[memberships[i]]][2] + 
        (pc_predictor_g[[memberships[i]]] %*% pred_scores)[2,1]
    } else  if (type == 'full') {
      pred_at_p[i] <- (mean_response_g[[memberships[i]]] + 
                         pc_response1_g[[memberships[i]]] %*%
                         pred_scores + 
                         pc_response2_g[[memberships[i]]] %*%
                         resp_scores)[1]
    } else if (type == 'mean') {
      coef <- parameters[['means_response']][[memberships[i]]]
      pred_at_p[i] <- mean_response_g[[memberships[i]]]
    }
  }
  pred_at_p
}

# same as predict_at_pressure, but predicts for each potential cluster membership
predict_at_pressure_all_groups <- function(p, parameters, pred, params, type = 'full') {
  sapply(1:params[['G']], function(g){ 
    prediction_data_subset <- pred
    cols <- colnames(prediction_data_subset)
    last_numbers <- sapply(strsplit(cols, '\\.'), function(x) if (length(x) <3) {F} else {x[[3]] == g})
    new_cols <- which(substr(cols, 1, 1) == 'E' &last_numbers)
    colnames(prediction_data_subset)[new_cols] <- apply(sapply(strsplit(cols[new_cols], split = '\\.'), function(r) r[1:2]),
                                                        2, paste0, collapse = '.')
    prediction_data_subset$membership = ifelse(prediction_data_subset$is_prediction,
                                               g,  prediction_data_subset$membership)
    predict_at_pressure(p, parameters, prediction_data_subset, params, type)
  })
}

# predict_at_pressure but for variances
# type - will specify temperature, salinity, or type of oxygen
# cond_type - if conditional or unconditional variances should be computed

variance_at_pressure <- function(pressure, parameters, var_scores, memberships, params,
                                        pred_df, type = 'full', cond_type = 'conditional'){
  Phi <- eval.basis(basisobj = params[['basis_response']], pressure)
  Phi_g <- lapply(1:params[['G']], function(g) {
    as.matrix(cbind(Phi %*% parameters[['Omega_response1']][[g]], Phi %*% parameters[['Omega_response2']][[g]] ))
  })
  Phi_g1 <- lapply(1:params[['G']], function(g) {
    as.matrix(Phi %*% parameters[['Omega_response1']][[g]])
  })
  Phi_g2 <- lapply(1:params[['G']], function(g) {
    as.matrix(Phi %*% parameters[['Omega_response2']][[g]])
  })
  Phi_predictors_g <- lapply(1:params[['G']], function(g) {
    as.matrix(bdiag(Phi, Phi) %*% parameters[['Omega_predictors']][[g]])
  })
  #check_diag <- sapply(var_scores, isDiagonal)
  check_diag <- rep(T, length(var_scores))
  
  var_scores_g <- lapply(1:params[['G']], function(g) {
    bdiag(Diagonal(n = params[['pc_response1']], parameters[['variances_predictors']][[g]]),
          Diagonal(n = params[['pc_response2']], parameters[['variances_response']][[g]]))
  })
  var_vec <- rep(0, length(memberships))
  for (x in 1:length(memberships)) {
    if (type == 'full') {
      Phi_all <- Phi_g[[memberships[x]]]
      if (cond_type == 'conditional') {
        var_scores_use <- var_scores[[x]]
      } else if (cond_type == 'unconditional') {
        var_scores_use <- var_scores_g[[memberships[x]]]
      }
    } else if (type == 'TS') {
      Phi_all <- Phi_g1[[memberships[x]]]
      if (cond_type == 'conditional') {
        var_scores_use <- var_scores[[x]][1:params[['pc_response1']], 1:params[['pc_response1']]]
      } else if (cond_type == 'unconditional') {
        var_scores_use <- var_scores_g[[memberships[x]]][1:params[['pc_response1']], 1:params[['pc_response1']]]
      }
    } else if (type %in% c('temperature', 'salinity')) {
      Phi_all <- Phi_predictors_g[[memberships[x]]]
      if (cond_type == 'conditional') {
        var_scores_use <- var_scores[[x]][1:params[['pc_response1']], 1:params[['pc_response1']]]
      } else if (cond_type == 'unconditional') {
        var_scores_use <- var_scores_g[[memberships[x]]][1:params[['pc_response1']], 1:params[['pc_response1']]]
      }
    } else if (type == 'spatial') {
      Phi_all <- Phi_g2[[memberships[x]]]
      if (cond_type == 'conditional') {
        var_scores_use <- var_scores[[x]][-(1:params[['pc_response1']]), -(1:params[['pc_response1']])]
      } else if (cond_type == 'unconditional') {
        var_scores_use <- var_scores_g[[memberships[x]]][-(1:params[['pc_response1']]), -(1:params[['pc_response1']])]
      }
    }
    
    if (type == 'temperature') {
      Phi_var <- as.matrix(Phi_all %*% var_scores_use)
      var_vec[x] <- sum(Phi_var[1,] * Phi_all[1,]) + parameters[['measurement_error_predictors']][1]
    } else if (type == 'salinity') {
      Phi_var <- as.matrix(Phi_all %*% var_scores_use)
      var_vec[x] <- sum(Phi_var[2,] * Phi_all[2,]) + parameters[['measurement_error_predictors']][2]
    } else {
      if (check_diag[x]) {
        var_vec[x] <- sum(var_scores_use@x * Phi_all^2) + parameters[['measurement_error_response']]
      } else {
        Phi_var <- as.matrix(Phi_all %*% var_scores_use)
        var_vec[x] <- sum(Phi_var * Phi_all) + parameters[['measurement_error_predictors']][1]
      }
    }
  }
  var_vec
}

# predict_at_pressure_all_groups but for variances
# type - will specify temperature, salinity, or type of oxygen
# cond_type - if conditional or unconditional variances should be computed
variance_at_pressure_all_groups <- function(p, parameters, pred, params, pred_var,
                                            type = 'full', cond_type = 'conditional') {
  sapply(1:params[['G']], function(g){ 
    variance_at_pressure(p, parameters, pred_var[[g]], 
                        memberships = rep(g, nrow(pred)),
                        params, pred_df, type = type, 
                        cond_type = cond_type)
  })
} 

# give both parts of variances
# returns a list with entry:
# 1: expected (conditional variance of variable on cluster)
# 2: variance (conditional expectation of variable on cluster)
# requires cond_probs to also be specified
compute_both_parts_of_variance <- function(p, parameters, pred_use, params, resp_scores_var_use, cond_probs,
                                           type = 'full',
                                           cond_type = 'conditional') {
  all_predictions <- variance_at_pressure_all_groups(p, parameters, pred_use, params, resp_scores_var_use, type = type,
                                                     cond_type = cond_type)
  all_predictions_expectation <- predict_at_pressure_all_groups(p, parameters, pred_use, params, 
                                                                type = ifelse(type == 'unconditional', 'mean', type))
  second_part_covariance <- sapply(1:nrow(pred_use), function(x) {
    mat_use <- -tcrossprod(cond_probs[x, ])
    diag(mat_use) <- cond_probs[x, ] * (1-cond_probs[x, ])
    as.double(t(all_predictions_expectation[x,]) %*% mat_use %*% all_predictions_expectation[x,])
  })
  list(rowSums(all_predictions * cond_probs), second_part_covariance)
}

# get grid info for the grid used, to specify height and width in plots

get_grid_info <- function(grid_size) {
  library(ncdf4)
  grid_file <- nc_open('paper/data/grid.nc', write = F)
  grid_long <- ncvar_get(grid_file, 'XC')[,1]
  grid_lat <- ncvar_get(grid_file, 'YC')[1,]
  grid_depth <- ncvar_get(grid_file, 'Depth')[seq(1, length(grid_long), by = round(6*grid_size)),
                                              seq(1, length(grid_lat), by = round(6*grid_size))]
  grid_long <- grid_long[seq(1, length(grid_long), by = round(6*grid_size))]
  grid_lat <- grid_lat[seq(1, length(grid_lat), by = round(6*grid_size))]
  
  height <- grid_lat[2:length(grid_lat)] - grid_lat[1:(length(grid_lat) - 1)]
  height <- data.frame("latitude" = grid_lat, height = c(height[1], height))
  
  grid_df <- cbind(expand.grid('longitude' = grid_long, 'latitude' = grid_lat),
                   depth = as.double(grid_depth))
  grid_df[['width']] <- grid_size
  grid_df <- merge(grid_df, height, by = "latitude")
}


