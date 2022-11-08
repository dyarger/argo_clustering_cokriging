initialize_clusters <- function(levels, n_profiles, var, G, init_strategy, ind_vals1,
                                ind_vals2 = NULL, resp_profs = NULL, preds=NULL, pred_profs=NULL,
                                profile_lengths=NULL, equal=T, init_cluster=NULL){
  
  if (!is.null(init_cluster)){
    return(init_cluster)
  } else if (init_strategy == 'random'){
    matrix(sample(1:G, n_profiles, replace = T), ncol = 1)
  } else if (init_strategy == 'kmeans'){
    if (is.null(preds)){
      init_mat = interpolate_onto_p_levels(levels = levels, ind_vals1 = ind_vals1,
                                          var = var, resp_profs = resp_profs)
    } else if (equal){
      init_mat = interpolate_onto_p_levels(levels = levels, ind_vals1 = ind_vals1,
                                          var = var, resp_profs = resp_profs,
                                          preds = preds, pred_profs = pred_profs, ind_vals2 = ind_vals2, 
                                          profile_lengths = profile_lengths)
    } else {
      init_mat = interpolate_onto_p_levels(levels=levels, ind_vals2=ind_vals2,
                                          preds=preds, pred_profs=pred_profs, profile_lengths=profile_lengths)
    }
    matrix(stats::kmeans(init_mat, centers = G)$cluster, ncol = 1)
  }
}

interpolate_onto_p_levels <- function(levels, ind_vals1, var=NULL, resp_profs=NULL, preds=NULL, 
                                      pred_profs=NULL, ind_vals2, profile_lengths=NULL) {
  
  ncol = length(do.call('c', levels))
  n = ifelse(is.null(resp_profs), length(pred_profs), length(resp_profs))
  init_mat <- matrix(0, ncol = ncol, nrow = n)
  
  start1=rep(1, n)
  
  if(!is.null(resp_profs)){
    for (i in 1:n){
      ind1 = start1[i]:(start1[i]+length(levels[[var]])-1)
      init_mat[i,ind1] = approx(ind_vals1[[i]], resp_profs[[i]], xout=levels[[var]], rule=2,
                                ties = list('ordered', mean))$y
    }
  }
  
  if(!is.null(pred_profs)){
    for (i in 1:n){
      start2 = 1
      for (j in 1:length(preds)){
        ind1 = start1[i]:(start1[i]+length(levels[[preds[j]]])-1)
        
        p_length = profile_lengths[i,j]
        ind2 = start2:(start2+p_length-1)
        
        init_mat[i,ind1]=approx(ind_vals2[[i]][ind2], pred_profs[[i]][ind2], xout=levels[[preds[j]]], rule=2,
                                ties = list('ordered', mean))$y
        
        start2 = start2 + p_length + 1
        start1[i] = start1[i] + length(levels[[preds[j]]]) + 1
      }
    }
  }
  init_mat
}

centered_response <- function(cluster_mat, profs, Phi_prof, Omegas, pcs) {
  lapply(1:length(profs), function(x){
    g_ind = cluster_mat[x,]
    as.numeric(profs[[x]] - Phi_prof[[x]] %*%  (Omegas[[g_ind]] %*% pcs[x,]))
  })
}

get_relevant_scores <- function(cond_mean, n_pc1, n_pc2, n){
  ret = matrix(0, nrow = n, ncol = n_pc2)
  for (i in 1:n){
    ret[i,] = as.numeric(cond_mean[((i-1)*n_pc1+1):((i-1)*n_pc1 + n_pc2)])
  }
  return(ret)
}

make_first_entry_positive <- function(m){
  for (i in 1:ncol(m)){
    if (m[1,i] < 0){m[,i] = -m[,i]}
  }
  return(m)
}


