prepare_EM_data <- function(data, params){
  
  data_inputs = list()

  var = params[['Y']]
  id = params[['id']]
  ind_name = params[['ind_name']]
  loc_names = params[['loc_names']]
  time = params[['time']]
  doy = params[['dayofyear']]
  
  if (is.data.frame(data)){
    gdf = split(data, data[,id])
    gdf_additional <- gdf
    data_inputs[['n_profiles']] = length(unique(data[,id]))
    data_inputs[['n_profiles_TS']] = data_inputs[['n_profiles']]
    data_inputs[['BGC']] = names(gdf) %in% names(gdf)
  } else {
    gdf = split(data[[1]], data[[1]][,id])
    gdf_additional = split(data[[2]], data[[2]][,id])
    data_inputs[['BGC']] = names(gdf_additional) %in% names(gdf)
    data_inputs[['n_profiles']] = length(gdf)
    data_inputs[['n_profiles_TS']] = length(unique(data[[2]][[id]]))
  }
  
  if (!is.null(time)){
    data_inputs[['time']] <- sapply(gdf_additional, function(x) julian(as.Date(x[[time]][1], format = '%Y-%m-%d'),
                                                            origin = as.Date('2000-01-01')))
  }
  
  if(!is.null(loc_names)){
    data_inputs[['locs']] <- t(sapply(gdf_additional, function(x) as.double(x[,loc_names][1,])))
    colnames(data_inputs[['locs']]) <- loc_names
  }
  
  if (!is.null(doy)){
    data_inputs[['day']] <- sapply(gdf_additional, function(x) as.double(x[,doy][1,]))
  }
    
  data_inputs[['nn_list']] = get_nearest_neighbors(loc_names = loc_names,
                                                   locs1 = data_inputs[['locs']],
                                                   nn = params[['nn']],
                                                   nn_range = params[['nn_range']],
                                                   nn_range_lat = params[['nn_range_lat']],
                                                   nn_range_lon = params[['nn_range_lon']],
                                                   nn_range_time = params[['nn_range_time']],
                                                   nn_type = params[['nn_type']],
                                                   nn_strategy = params[['nn_strategy']],
                                                   nn_lat_param = params[['nn_lat_param']],
                                                   nn_time_param = params[['nn_time_param']],
                                                   time1 = data_inputs[['time']],
                                                   day1 = data_inputs[['day']],
                                                   remove_between_land = params[['remove_between_land']])
  
  data_inputs[['response_profs']] = lapply(gdf, function(x) as.double(unlist(x[!is.na(x[,var]), var])))
  data_inputs[['ind_response']] = lapply(gdf, function(x) as.double(unlist(x[!is.na(x[,var]), ind_name])))
    
  data_inputs[['Phi_response_prof']] = get_phi_prof(indexes = data_inputs[['ind_response']], 
                                                    kind = 'response', 
                                                    basis_list = params[['basis_response']])

  data_inputs[['Phi_x_Phi_response']] <- lapply(data_inputs[['Phi_response_prof']], function(x) as(Matrix::crossprod(x), 'dgCMatrix'))
  
  data_inputs[['penalty_mat_response']] <- get_pen_mat('response', params[['basis_response']])
  data_inputs[['inner_prod_mat_response']] <- as.matrix(Matrix::chol(fda::inprod(params[['basis_response']], params[['basis_response']])))
  data_inputs[['nbasis_response']] = ncol(data_inputs[['Phi_response_prof']][[1]])
  
  resp_lengths = matrix(0, nrow = data_inputs[['n_profiles']], ncol = 1)
  counter = 1
  for (x in data_inputs[['response_profs']]){
    resp_lengths[counter,] = length(x)
    counter = counter + 1
  }
  data_inputs[['profile_lengths_response']] = resp_lengths
  
  
  if (!is.null(params[['X']])){
    
    preds = params[['X']]
    
    data_inputs[['predictor_profs']] = list() 
    data_inputs[['ind_predictors']] = list()
    
    profile_lengths = matrix(0, nrow = data_inputs[['n_profiles_TS']], ncol = length(preds));
    counter = 1
    if (is.data.frame(data)){
      for(x in gdf){
        ret1 = c()
        ret2 = c()
        for (i in 1:length(preds)){
          t_df = x[!is.na(x[,preds[i]]),]
          ret1 = c(ret1, as.double(unlist(t_df[,preds[i]])))
          ret2 = c(ret2, as.double(unlist(t_df[,ind_name])))
          profile_lengths[counter, i] = nrow(t_df)
        }
        data_inputs[['predictor_profs']][[counter]] = ret1
        data_inputs[['ind_predictors']][[counter]] = ret2
        counter = counter + 1
      }
    } else {
      for(x in gdf_additional){
        ret1 = c()
        ret2 = c()
        for (i in 1:length(preds)){
          t_df = x[!is.na(x[,preds[i]]),]
          ret1 = c(ret1, as.double(unlist(t_df[,preds[i]])))
          ret2 = c(ret2, as.double(unlist(t_df[,ind_name])))
          profile_lengths[counter, i] = nrow(t_df)
        }
        data_inputs[['predictor_profs']][[counter]] = ret1
        data_inputs[['ind_predictors']][[counter]] = ret2
        counter = counter + 1
      }
    }
    
    data_inputs[['profile_lengths_predictors']] = profile_lengths
    data_inputs[['Phi_predictors_prof']] = get_phi_prof(indexes = data_inputs[['ind_predictors']], 
                                                        kind = 'predictors', 
                                                        basis_list = params[['basis_predictors']],
                                                        prof_ns = data_inputs[['profile_lengths_predictors']])
    data_inputs[['Phi_x_Phi_predictors']] <- lapply(data_inputs[['Phi_predictors_prof']], function(x) as(Matrix::crossprod(x), 'dgCMatrix'))
    
   
    data_inputs[['penalty_mat_predictors']] <- get_pen_mat('predictors', params[['basis_predictors']])
    data_inputs[['nbasis_predictors']] = sapply(params[['basis_predictors']], function(i){
      i$nbasis
    })
    data_inputs[['inner_prod_mat_predictors']] <-  as.matrix(
        Matrix::bdiag(lapply(1:length(params[['basis_predictors']]), function(x) {
          Matrix::chol(fda::inprod(params[['basis_predictors']][[x]], params[['basis_predictors']][[x]]))
        })
    ))
  }
  return(data_inputs)
}


get_nearest_neighbors <- function(locs1, loc_names, nn, nn_range, nn_range_lat, nn_range_lon, nn_range_time,
                                  nn_type, nn_strategy, nn_lat_param, nn_time_param, locs2 = NULL,
                                  time1 = NULL, time2 = NULL, day1 = NULL, day2 = NULL,
                                  remove_between_land = F){
  if (is.null(remove_between_land)) {
    remove_between_land <- F
  }
  one_set_of_locations <- is.null(locs2)
  
  if (!one_set_of_locations) {
    one_set_of_locations <- nrow(locs2) == nrow(locs1)
  }
  
  if (one_set_of_locations) {
    locs2 <- locs1
    time2 <- time1
    day2 <- day1
  }
  
  neighbors <- list()
  locs1_long_lat <- as.matrix(locs1[, loc_names])
  locs2_long_lat <- as.matrix(locs2[, loc_names])
  locs1_long_lat[,1] <- ifelse(locs1_long_lat[,1] > 180, locs1_long_lat[,1] - 360, locs1_long_lat[,1]) 
  locs2_long_lat[,1] <- ifelse(locs2_long_lat[,1] > 180, locs2_long_lat[,1] - 360, locs2_long_lat[,1]) 
  
  mean_vals <- colMeans(locs1_long_lat)
  
  for (i in 1:nrow(locs2)) {
    if (nn_type == 'space') {
      dist <- fields::rdist.earth.vec(locs2_long_lat[i,, drop = F],
                                      locs1_long_lat, miles = F)
    } else if (nn_type == 'spacewarp') {
      lon_dist <- fields::rdist.earth.vec(cbind(locs2_long_lat[i,1], mean_vals[2]),
                                          cbind(locs1_long_lat[,1], mean_vals[2]), miles = F)
      lat_dist <- fields::rdist.earth.vec(cbind( mean_vals[1], locs2_long_lat[i,2]),
                                          cbind( mean_vals[1], locs1_long_lat[,2]), miles = F)
      dist <- lon_dist + nn_lat_param*lat_dist
    } else if (nn_type == 'spacewarp_time') {
      lon_dist <- fields::rdist.earth.vec(cbind(locs2_long_lat[i,1], mean_vals[2]),
                                          cbind(locs1_long_lat[,1], mean_vals[2]), miles = F)
      lat_dist <- fields::rdist.earth.vec(cbind( mean_vals[1], locs2_long_lat[i,2]),
                                          cbind( mean_vals[1], locs1_long_lat[,2]), miles = F)
      time_dist <- abs(time1 - time2[i])
      dist <- lon_dist + nn_lat_param*lat_dist + nn_time_param*time_dist
    } else if (nn_type == 'space_time') {
      dist <- fields::rdist.earth.vec(locs2_long_lat[i,, drop = F],
                                      locs1_long_lat, miles = F)
      time_dist <- abs(time1 - time2[i])
    } else if (nn_type == 'space_day') {
      dist <- fields::rdist.earth.vec(locs2_long_lat[i,, drop = F],
                                      locs1_long_lat, miles = F)
      time_dist1 <- abs(day1 - day2[i])
      time_dist2 <- abs(day1 - day2[i] + 365.25)
      time_dist3 <- abs(day1 - day2[i] - 365.25)
      time_dist <- time_dist3
      b1 <- time_dist1 < time_dist2 & time_dist1 < time_dist3
      b2 <- time_dist2 < time_dist1 & time_dist2 < time_dist3
      time_dist[b1] <- time_dist1[b1]
      time_dist[b2] <- time_dist2[b2]
      dist <- dist + nn_time_param*time_dist
    } else if (nn_type == 'spacewarp_day') {
      lon_dist <- fields::rdist.earth.vec(cbind(locs2_long_lat[i,1], mean_vals[2]),
                                          cbind(locs1_long_lat[,1], mean_vals[2]), miles = F)
      lat_dist <- fields::rdist.earth.vec(cbind( mean_vals[1], locs2_long_lat[i,2]),
                                          cbind( mean_vals[1], locs1_long_lat[,2]), miles = F)
      time_dist1 <- abs(day1 - day2[i])
      time_dist2 <- abs(day1 - day2[i] + 365.25)
      time_dist3 <- abs(day1 - day2[i] - 365.25)
      time_dist <- time_dist3
      b1 <- time_dist1 < time_dist2 & time_dist1 < time_dist3
      b2 <- time_dist2 < time_dist1 & time_dist2 < time_dist3
      time_dist[b1] <- time_dist1[b1]
      time_dist[b2] <- time_dist2[b2]
      dist <- lon_dist + nn_lat_param*lat_dist + nn_time_param*time_dist
    }
    
    if (one_set_of_locations) {dist[i] <- Inf}
    if (remove_between_land) {
      dist[!remove_land_locs(locs1_long_lat, locs2_long_lat[i,])] <- Inf
    }

    if (nn_strategy == 'nn') {
      value <- Rfast::nth(dist, k = nn)
      neighbors[[i]] <- which(dist <= value)[1:nn]
    } else if (nn_strategy == 'range') {
      if (nn_type == 'space') {
        neighbors[[i]] <- which(dist <= nn_range)
      } else if (nn_type == 'space_warp') {
        neighbors[[i]] <- which(lon_dist <= nn_range_lon & 
                                lat_dist <= nn_range_lat & 
                                dist != Inf)
      } else if (nn_type == 'space_time') {
        neighbors[[i]] <- which(lon_dist <= nn_range_lon & 
                                lat_dist <= nn_range_lat & 
                                time_dist <= nn_range_time & 
                                dist != Inf)
      }
    } else if (nn_strategy == 'nn_range') {
      dist_use_indexes <- which(time_dist <= nn_range_time)
      dist_use <- dist[dist_use_indexes]
      value <- Rfast::nth(dist_use, k = nn)
      neighbors[[i]] <- dist_use_indexes[which(dist_use <= value)]
    }
  }
  
  if (one_set_of_locations) {
    neighbors = make_graph_undirected(neighbors)
  }
  
  neighbors
}

make_graph_undirected <- function(neighbors, union = T){
  for (i in 1:length(neighbors)){
    for (j in 1:length(neighbors[[i]])){
      nbs = neighbors[[i]]
      if (!(i %in% neighbors[[nbs[j]]])){
        neighbors[[nbs[j]]] = c(neighbors[[nbs[j]]], i)
      }
    }
  }
  return(neighbors)
}

get_phi_prof <- function(indexes, kind, basis_list, knots = NULL, prof_ns = NULL, sparse=T){
  if (kind == 'response'){
    return_value = lapply(indexes, function(x){
      if (sparse){
        if (!is.null(knots)) {
          as(splineMatrixC(order = 4, knots = knots, x = x), 'dgCMatrix')
        } else {
          as(fda::eval.basis(x, basis_list), 'dgCMatrix')
        }
      } else {
        as.matrix(fda::eval.basis(x, basis_list))
      }
    })
  } else {
    return_value = lapply(1:length(indexes), function(x){
      start = 1
      ret = list()
      for (i in 1:length(basis_list)){
        ind = start:(start+prof_ns[x,i]-1)
        if (sparse){
          if (!is.null(knots)) {
            ret[[i]] = as(splineMatrixC(order = 4, knots = knots, x = x), 'dgCMatrix')
          } else {
            ret[[i]] = as(fda::eval.basis(indexes[[x]][ind], basis_list[[i]]), 'dgCMatrix')
          }
        } else {
          ret[[i]] = as.matrix(fda::eval.basis(indexes[[x]][ind], basis_list[[i]]))
        }
        start = prof_ns[x,i] + 1
      }
      Matrix::bdiag(ret)
    })
  }
}

get_pen_mat<- function(kind, basis_list){
  if (kind == 'response'){
    return_value = as(fda::inprod(basis_list, basis_list, 2, 2), 'sparseMatrix')
  } else {
    return_value = Matrix::bdiag(lapply(1:length(basis_list), function(x) {
      as(fda::inprod(basis_list[[x]], basis_list[[x]], 2, 2), 'sparseMatrix')
    }) )
  }
}
remove_land_locs <- function(locs1, locs2_vec) {
  longitude <- locs1[,1]
  latitude <- locs1[,2]
  long_0 <- locs2_vec[1]
  lat_0 <- locs2_vec[2]
  # it is bad if (-100 < long1 < -60), and (-60 < lon2 < -20)
  # it is bad if (long1 < -20), and (long2 > -100)
  one <- !((longitude < -60 & longitude > -100 & latitude < -62.5) & 
             (long_0> -60 & long_0 < -20))
  two <- !((longitude > -60  & longitude < -20 & latitude < -62.5) & 
             (long_0 < -60 & long_0 > -100))
  three <- !((longitude < -73 & longitude > -115 & latitude > -55) &
               (long_0 > -73 & long_0 < -35))
  four <- !((longitude > -73 & longitude < -35 & latitude > -55) & 
              (long_0 < -73 & long_0 < -115))
  one & two & three & four 
}


