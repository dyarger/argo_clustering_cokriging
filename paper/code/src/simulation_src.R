library(PottsUtils)
library(fields)
library(FNN)
library(ggplot2)
library(fda)
library(mvtnorm)

simulate_data <- function(params=list(), show_locs=F, show_pcs=F, show_profs=F){
  
  params = init_simu_params(params)
  
  # Frequently used
  index_df = params$locs
  G = params$G
  n_pcs = params$n_pcs
  
  dists = compute_dists(index_df, params$deptype, params$anisotropic, params$use_gcd)
  
  if(params$G > 1){
    if (params$deptype!='spatial' | !is.null(params$dist_weights)){
      agg_dists = aggregate_dists(dists, params$dist_weights)
      Z = simu_mrf(G, agg_dists, params$theta, params$knn, params$range)
    } else {
      Z = simu_mrf(G, dists, params$theta, params$knn, params$range)
    } 
  } else {
    Z = rep(1, params$n)
  }
  
  index_df$cluster = as.factor(Z)
  
  for (q in 1:n_pcs[1]){
    for (g in 1:G){
      name = paste0('pc_resp', as.character(q))
      ind = which(index_df$cluster==g)
      index_df[ind, name] = simu_single_scores(index_df[ind,],
                                               params$deptype,
                                               variance=params$var_resp[g,q],
                                               st_range=params$st_range_resp[[g]][[q]],
                                               params$anisotropic)
    } 
  }
  
  if(params$var > 1){
    names1 = c()
    names2 = c()
    # for(q in 1:n_pcs[3]){
    #   names2 = c(names2, paste0('pc_pred2', as.character(q)))
    # }
    for (g in 1:G){
      ind = which(index_df$cluster==g)
      names1 = c()
      for (q in 1:n_pcs[2]){
        name = paste0('pc_pred', as.character(q))
        names1 = c(names1, name)
        index_df[ind, name] = simu_single_scores(index_df[ind,], params$deptype, variance=params$var_pred[[1]][g,q],
                                                 st_range=params$st_range_pred[[g]][[q]], params$anisotropic)
      }
      # index_df[ind,names2] = 0
      # index_df[ind, names2] = simu_paired_scores(as.matrix(index_df[ind,names1]),
      #                                            params$var_pred[[1]][g,],
      #                                            params$var_pred[[2]][g,], params$C[[g]])
    }
  }
  
  data = create_data(index_df, params$G, params$var, params$n_pcs, params$means,
                     params$pc_funs, params$lt, params$domain,
                     params$n_obs, params$random_n_obs,
                     params$random_locs, params$me_var, params$var_pred[[1]],
                     params$var_pred[[2]], params$C)
  #params$proj_var = data$proj_var
  
  data$params = params
  
  visualize(data, show_locs, show_pcs, show_profs)
  
  return(data)
}

simu_locs <- function(n, sp_range, time_range, deptype){
  lat = runif(n, min = sp_range[[1]][1], max = sp_range[[1]][2])
  lon = runif(n, min = sp_range[[2]][1], max = sp_range[[2]][2])
  if (deptype != 'spatial'){
    time = runif(n, min = time_range[1], max = time_range[2]) 
    return(data.frame(time=time, lat=lat, lon=lon))
  }
  return(data.frame(lat=lat, lon=lon))
}

compute_dists <-function(locs, deptype, anisotropic=FALSE, use_gcd=F){
  if (deptype == 'spatial'){
    if (!anisotropic){
      if (!use_gcd){
        return(space = as.matrix(fields::rdist(as.matrix(locs))))
      } else {
        return(space = as.matrix(fields::rdist.earth(as.matrix(locs), miles = F)))
      } 
    } else {
      return(list(
        lat = as.matrix(fields::rdist(as.matrix(locs$lat))),
        lon = as.matrix(fields::rdist(as.matrix(locs$lon)))
      ))
    }
  } else {
    if (!anisotropic){
      if(!use_gcd){
        return(list(
          space = as.matrix(fields::rdist(as.matrix(locs[,c('lat', 'lon')]))),
          time = as.matrix(fields::rdist(as.matrix(locs[,'time'])))
        ))
      } else {
        return(list(
          space = as.matrix(fields::rdist.earth(as.matrix(locs[,c('lat', 'lon')]), miles = F)),
          time = as.matrix(fields::rdist(as.matrix(locs$time)))
        ))
      }
    } else {
      return(list(
        lat = as.matrix(fields::rdist(as.matrix(locs$lat))),
        lon = as.matrix(fields::rdist(as.matrix(locs$lon))),
        time = as.matrix(fields::rdist(as.matrix(locs$time)))
      ))
    }
  }
}

simu_mrf <- function(G, dists, theta, k=NULL, st_range=NULL) {
  nns = list()
  if(is.null(st_range)){
    nn = FNN::knn.index(dists, k)
    for (i in 1:nrow(nn)){
      nns[[i]] = nn[i,]
    }
    n = nrow(nn)
  } else {
    if (class(dists[1]) == "list"){
      n = nrow(dists[[1]])
    } else {
      n = nrow(dists)
    }
    ind = matrix(1, nrow = n, ncol = n)
    for (name in names(dists)){
      ind = ind * (dists[[name]] < st_range[[name]])
    }
    for (i in 1:n){
      nns[[i]] = which(ind[i,] == 1)
    }
  }
  
  edges = NULL
  for(i in 1:n){
    ind = which(nns[[i]] > i)
    for(j in ind){
      edges = rbind(edges,c(i,nns[[i]][j]))
    }
  }
  
  mrf = SW2(1, nvertex=n, ncolor=G, edges, beta=theta)
  iter = 1
  issue_sims = which(apply(mrf, 2, var) == 0)
  while(length(issue_sims)>0){
    mrf[,issue_sims] = SW2(n=sum(issue_sims), nvertex=n, ncolor=G, edges, beta=theta)
    issue_sims = which(apply(mrf, 2, var) == 0)
    iter = iter + 1
    if(iter > 50) stop("Couldn't produce enough meaningful MRF's in time.")
  }
  mrf
}

simu_single_scores <- function(index_df, deptype = 'spatial', variance=1, st_range=1, anisotropic=F, use_gcd=F){
  n = nrow(index_df)
  if (anisotropic | deptype != "spatial"){
    dist = matrix(0, nrow = n, ncol = n)
    for (name in c('lat', 'long', 'time')){
      dist = dist + as.matrix(dist(index_df[[name]]))^2 / st_range[[name]]
    }
    dist = sqrt(dist)
  } else {
    if (!use_gcd){
      dist = as.matrix(fields::rdist(as.matrix(index_df[,c('lat', 'lon')]))) / st_range
    } else {
      dist = as.matrix(fields::rdist.earth(as.matrix(index_df[,c('lat', 'lon')]), miles = F)) / st_range
    } 
  }
  cov = matrix(Exponential(as.numeric(dist), phi = variance), ncol = n)
  L = chol(cov)
  return(crossprod(L, rnorm(n)))
}

simu_paired_scores <- function(scores, vars1, vars2, C){
  return(tcrossprod(scores, C) + mvtnorm::rmvnorm(nrow(scores), sigma=diag(vars2) - C %*% diag(vars1) %*% t(C)))
}

create_pc_funs <- function(G, n_pcs, range){
  pcs = list()
  n_knots = 2*n_pcs + 2
  knots = seq(range[1], range[2], length.out = n_knots)
  basis = fda::create.bspline.basis(knots)
  inprod_mat = inprod(basis, basis)
  inv_inprod_chol = solve(chol(inprod_mat))
  coefs_mat = matrix(c(1,5,2,5,1,2), ncol = 2)
  for (g in 1:G){
    pcs[[g]] = list()
    for (q in 1:n_pcs){
      #which_pcs = sample(1:supp, replace = F, size=supp-1)
      #which_pcs = seq(3+(q-1)*supp, 3+(supp-1)+(q-1)*supp, by = 1)
      which_pcs = c(2 + (q-1)*2, 3 + (q-1)*2)
      #coefs = runif(length(which_pcs), min = 4, max = 5)
      #coefs = as.numeric(rmvnorm(1, mean = c(4,4), sigma = matrix(c(4,-3.5, -3.5, 4), ncol=2)))
      coefs = coefs_mat[g,]
      coefs = coefs / norm(coefs)
      #coefs_temp = c(tail(coefs, -shifts[q,g]), head(coefs, shifts[q,g]))
      pc = list(basis=basis, ortho=inv_inprod_chol, which_pcs=which_pcs, coefs=coefs)
      #if(evaluate_pcs(range[1], pc) < 0){
      #  coefs = -coefs
      #}
      pcs[[g]][[q]] = list(basis=basis, ortho=inv_inprod_chol, which_pcs=which_pcs, coefs=coefs)
    }
  }
  pcs
}

evaluate_pcs <- function(p, pc){
  if (class(pc) == "function"){
    return(pc(p))
  } else {
    if (length(pc$coefs) == 1){
      return(((fda::eval.basis(p, pc$basis) %*% pc$ortho)[,pc$which_pcs] * pc$coefs))
    }
    return(((fda::eval.basis(p, pc$basis) %*% pc$ortho)[,pc$which_pcs] %*% pc$coefs))
  }
}

create_mean_functions <- function(var, G){
  
  means = list()
  means[[1]] = function(p) {2 * exp(2-2*p) * cos(5*(p-0.1))/2}
  means[[2]] = function(p) {3*cos(1.5*p)}
  means[[3]] = function(p) {4*exp(p) * cos(p)/2}
  means[[4]] = function(p) {5*cos(5*p)} 
  if (var != 1){
    means_ret = list()
    means_ret[['resp']] = means[1:G]
    means_ret[['pred1']] = means[1:G]
    means_ret[['pred2']] = means[4:(4-G)]
    return(means_ret)
  } else {
    return(means[1:G])
  }
  
  # means = list()
  # means[[1]] = function(p) {2 * exp(2-2*p) * cos(5*(p-0.1))/2}
  # means[[2]] = function(p) {5*cos(5*p)-2} 
  # means[[3]] = function(p) {5*sin(5*p)} 
  # means[[4]] = function(p) {3*cos(1.5*p)+3}
  # means[[5]] = function(p) {4*exp(p)*cos(p)/2}
  # means[[6]] = function(p) {cos(5*p)+3} 
  # means[[7]] = function(p) {3*sin(1.5*p)+3}
  # means[[8]] = function(p) {4*exp(p)*sin(p)/2}
  # means[[9]] = function(p) {sin(5*p)+3} 
  # if (var != 1){
  #   means_ret = list()
  #   means_ret[['resp']] = means[1:G]
  #   means_ret[['pred1']] = means[(G+1):(2*G)]
  #   means_ret[['pred2']] = means[(2*G+1):(3*G)]
  #   return(means_ret)
  # } else {
  #   return(means[1:G])
  # }
}

create_lt_functions <- function(G, n_pcs){
  lts = list()
  lts[[1]] = function(p){return(sin(6*p))}
  lts[[2]] = function(p){return(cos(6*p))}
  lts[[3]] = function(p){return(exp(p)*cos(5*p))}
  lts[[4]] = function(p){return(2*sin(4*p)*p*(1/(p+1) + 0.5))}
  lts[[5]] = function(p){return(sqrt(3)*(0.8-(p-0.1)^2) + cos(4*p))}
  lts[[6]] = function(p){return(4*(0.8-(p-0.5)^2) - 2)}
  lts[[7]] = function(p){return(sqrt(6)*p*cos(5*p))}
  lts[[8]] = function(p){return(0.9*exp(p)*sin(5*p))}
  lts[[9]] = function(p){return(-(0.75-p)^2*(p-0.6)*5)}
  lts[[10]] = function(p){return(1/(p-1.3) + 1)}
  lts[[11]] = function(p){return(8*(p-0.5)^2 * (p-0.75))}
  lts[[12]] = function(p){return(((1.2*p-0.5)^4) * 4)}
  lts_ret = list()
  
  lts_ret = list()
  for (g in 1:G){
    lts_ret[[g]] = list()
    for (q in 1:n_pcs){
      lts_ret[[g]][[q]] = lts[[(g-1)*G + q]]
    }
  }
  return(lts_ret)
}

create_data <- function(index_df, G, var, n_pcs, mean_functions, pc_functions, lt, range,
                        n_obs, random_n_obs, random_locs, me_var, vars1, vars2, C){
  data_df = data.frame()
  if (var == 1){
    proj_var <- list()
    for(i in 1:nrow(index_df)){
      g_temp = index_df$cluster[i]
      n_pcs = length(pc_functions[[g_temp]])
      if(random_n_obs){
        temp_n = n_obs + sample(1, c(-1,1)) * sample(1, 1:round(0.1*n_obs))
      } else {
        temp_n = n_obs
      }
      
      if(random_locs){
        ps = sort(runif(temp_n, min = range[1], max = range[2]))
      } else {
        ps = seq(range[1], range[2], length.out = temp_n)
      }
      
      obs = mean_functions[[g_temp]](ps) 
      for (q in 1:n_pcs){
        name = paste0('pc_resp', as.character(q))
        obs = obs + as.numeric(index_df[i, name]) * evaluate_pcs(ps, pc_functions[[g_temp]][[q]])
      }
      obs = obs + rnorm(temp_n, sd = sqrt(me_var[1]))
      
      temp_df = data.frame(pressure = ps,
                           var1 = obs,
                           profile_unique = i,
                           latitude = index_df$lat[i],
                           longitude = index_df$lon[i],
                           cluster = index_df$cluster[i])
      if(!is.null(index_df$time[i])){
        temp_df$time = index_df$time[i]
      }
      data_df = rbind(data_df, temp_df)
    }
  } else {
    
    # proj = list()
    # proj_var = list()
    # for (g in 1:G){
    #   cov_mat = rbind(cbind(diag(vars1[g,]), C[[g]] %*% diag(vars2[g,])), cbind(tcrossprod(diag(vars2[g,]), C[[g]]), diag(vars2[g,])))
    #   eig = eigen(cov_mat)
    #   proj[[g]] = eig$vectors
    #   if (min(eig$values) <= 0){
    #     stop("Non positive definiteness problem")
    #   } 
    #   proj_var[[g]] = eig$values
    # }
    
    
    # names = c()
    # for(j in 1:(var-1)){
    #   for (q in 1:n_pcs[j+1]){
    #     name = paste0('pc_pred', as.character(j), as.character(q))
    #     names = c(names, name)
    #   } 
    # }
    
    for(i in 1:nrow(index_df)){
      g_temp = index_df$cluster[i]
      
      if(random_n_obs){
        temp_n = n_obs + sample(1, c(-1,1)) * sample(1, 1:round(0.1*n_obs))
      } else {
        temp_n = n_obs
      }
      
      if(random_locs){
        ps = sort(runif(temp_n, min = range[['resp']][1], max = range[['resp']][2]))
      } else {
        ps = seq(range[['resp']][1], range[['resp']][2], length.out = temp_n)
      }
      
      obs = mean_functions[['resp']][[g_temp]](ps) 
      
      for (q in 1:n_pcs[2]){
        coef = index_df[i,paste0('pc_pred', as.character(q))]
        obs = obs + as.numeric(coef) * lt[[g_temp]][[q]](ps)
      }
      for (q in 1:n_pcs[1]){
        name = paste0('pc_resp', as.character(q))
        obs = obs + as.numeric(index_df[i, name]) * evaluate_pcs(ps, pc_functions[['resp']][[g_temp]][[q]])
      }
      obs = obs + rnorm(temp_n, sd = sqrt(me_var[1]))
      
      temp_df = data.frame(pressure = ps,
                           var1 = obs,
                           profile_unique = i,
                           latitude = index_df$lat[i],
                           longitude = index_df$lon[i],
                           cluster = index_df$cluster[i],
                           pred1 = NA, pred2 = NA)
      if(!is.null(index_df$time[i])){
        temp_df$time = index_df$time[i]
      }
      
      data_df = rbind(data_df, temp_df)
      
      for(j in 1:2){
        
        sname = paste0('pred', as.character(j))
        
        if(random_n_obs){
          temp_n = n_obs + sample(1, c(-1,1)) * sample(1, 1:round(0.1*n_obs))
        } else {
          temp_n = n_obs
        }
        
        if(random_locs){
          ps = sort(runif(temp_n, min = range[[sname]][1], max = range[[sname]][2]))
        } else {
          ps = seq(range[[sname]][1], range[[sname]][2], length.out = temp_n)
        }
        
        obs = mean_functions[[sname]][[g_temp]](ps) 
        
        for (q in 1:n_pcs[2]){
          name = paste0('pc_pred', as.character(q))
          obs = obs + as.numeric(index_df[i, name]) * evaluate_pcs(ps, pc_functions[[sname]][[g_temp]][[q]])
        }
        
        obs = obs + rnorm(temp_n, sd = sqrt(me_var[j+1]))
        
        temp_df = data.frame(pressure = ps,
                             profile_unique = i,
                             latitude = index_df$lat[i],
                             longitude = index_df$lon[i],
                             cluster = index_df$cluster[i],
                             var1 = NA, pred1 = NA, pred2 = NA)
        temp_df[,sname] = obs
        
        if(!is.null(index_df$time[i])){
          temp_df$time = index_df$time[i]
        }
        data_df = rbind(data_df, temp_df)
      }
    }
  }
  return(list(data_df=data_df, index_df=index_df))
}

init_simu_params <- function(params){
  if(is.null(params[['n']])){
    params$n = 200
  }
  if(is.null(params$G)){
    params$G = 3
  } 
  if(is.null(params$deptype)){
    params$deptype = 'spatial'
  }
  if(is.null(params[['var']])){
    params$var = 1
  }
  if(is.null(params$n_pcs)){
    if(params$var > 1){
      params$n_pcs = c(2,2,2,4)
    } else {
      params$n_pcs = 3
    }
  }
  if(is.null(params$sp_range)){
    params$sp_range = list(c(0,1), c(0,1))
  }
  if(is.null(params$time_range) & params$deptype != 'spatial'){
    params$time_range = c(0,1)
  }
  if(is.null(params$domain)){
    if (params$var > 1){
      params$domain[['resp']] = c(0,1)
      params$domain[['pred1']] = c(0,1)
      params$domain[['pred2']] = c(0,1)
    } else {
      params$domain = c(0,1)
    }
  }
  if(is.null(params$n_obs)){
    params$n_obs = 50
  }
  if(is.null(params$random_n_obs)){
    params$random_n_obs = F
  }  
  if(is.null(params$random_locs)){
    params$random_locs = F
  }
  if (is.null(params$locs)){
    params$locs = simu_locs(params$n, params$sp_range, params$time_range, params$deptype)
  }
  if(is.null(params$var_resp)){
    params$var_resp = matrix(rep(params$n_pcs[1]:1, params$G), nrow = params$G, byrow = T) / params$n_pcs[1]
  }
  if(is.null(params$st_range_resp)){
    for (g in 1:params$G){
      params$st_range_resp[[g]] = list()
      for (q in 1:params$n_pcs[1]){
        params$st_range_resp[[g]][[q]] = runif(1)/100
      }
    }
  }
  if(params$var > 1){
    if(is.null(params$var_pred)){
      params$var_pred[[1]] = matrix(rep(params$n_pcs[2]:1, params$G), nrow = params$G, byrow = T)
    }  
    if(is.null(params$st_range_pred)){
      for (g in 1:params$G){
        params$st_range_pred[[g]] = list()
        for (q in 1:params$n_pcs[2]){
          params$st_range_pred[[g]][[q]] = runif(1)/100
        }
      }
    }
    # if(is.null(params$C)){
    #   params$C = list()
    #   for(g in 1:params$G){
    #     params$C[[g]] = 1/10*(matrix(rnorm(params$n_pcs[2] * params$n_pcs[3]), ncol = params$n_pcs[2], nrow = params$n_pcs[3]))
    #   } 
    # } 
  }
  if (is.null(params$pc_funs)){
    if (params$var == 1){
      params$pc_funs = create_pc_funs(params$G, params$n_pcs, params$domain)
    } else {
      params$pc_funs = list()
      params$pc_funs[['resp']] = create_pc_funs(params$G, params$n_pcs[1], params$domain[['resp']])
      params$pc_funs[['pred1']] = create_pc_funs(params$G, params$n_pcs[2], params$domain[['pred1']])
      params$pc_funs[['pred2']] = create_pc_funs(params$G, params$n_pcs[2], params$domain[['pred2']])
    }
  } 
  if(params$var > 1){
    if(is.null(params$lt)){
      params$lt = create_lt_functions(params$G, params$n_pcs[2])
    } 
  }
  if (is.null(params$means)){
    params$means = create_mean_functions(params$var, params$G)
  }
  if(is.null(params$theta)){
    params$theta = 0.6
  }
  if(is.null(params$anisotropic)){
    params$anisotropic=F
  }
  if(is.null(params$use_gcd)){
    params$use_gcd=F
  }
  if(is.null(params$knn)){
    params$knn = 5
  }
  if(is.null(params$me_var)){
    if (params$var > 1){
      params$me_var = rep(0.2, 3)
    } else {
      params$me_var = 0.2
    }
  }
  return(params)
}

visualize <- function(data, show_locs, show_pcs, show_profs){
  params = data$params
  data_df = data$data_df
  index_df = data$index_df
  if(show_locs){
    print(ggplot(index_df[,c('lat', 'lon', 'cluster')], aes(x = lat, y = lon, color = cluster)) +
            geom_point() + theme_bw())
  }
  
  if(show_pcs){
    if (params$var == 1){
      ps = runif(200, min=params$domain[1], max=params$domain[2])
      t_df = data.frame()
      for(g in 1:params$G){
        for(q in 1:params$n_pcs[1]){
          tt_df = data.frame(p = ps)
          tt_df[['Value']] = evaluate_pcs(ps, params$pc_funs[[g]][[q]])
          tt_df[['Cluster']] = g
          tt_df[['PC']] = q
          tt_df[['WhichPC']] = paste0(as.character(g), as.character(q))
          t_df = rbind(t_df, tt_df)
        } 
      }
      t_df$Cluster = as.factor(t_df$Cluster)
      t_df$PC = as.factor(t_df$PC)
      print(ggplot(t_df, aes(x = p, y = Value, color = Cluster, group = WhichPC, linetype = PC)) + 
              geom_line() + theme_bw())
    } else {
      ps = runif(200, min=params$domain[['resp']][1], max=params$domain[['resp']][2])
      t_df = data.frame()
      for(g in 1:params$G){
        for(q in 1:params$n_pcs[1]){
          tt_df = data.frame(p = ps)
          tt_df[['Value']] = evaluate_pcs(ps, params$pc_funs[['resp']][[g]][[q]])
          tt_df[['Cluster']] = g
          tt_df[['PC']] = q
          tt_df[['WhichPC']] = paste0(as.character(g), as.character(q))
          t_df = rbind(t_df, tt_df)
        } 
      }
      t_df$Cluster = as.factor(t_df$Cluster)
      t_df$PC = as.factor(t_df$PC)
      print(summary(t_df$Cluster))
      print(ggplot(t_df, aes(x = p, y = Value, color = Cluster, group = WhichPC, linetype = PC)) + 
              geom_line() + theme_classic())
    }
  }
  
  if(show_profs){
    ind = sample(1:nrow(index_df), size = 30)
    t_df = data_df[(data$data_df)$profile_unique %in% ind,]
    print(ggplot(t_df[!is.na(t_df[,'var1']),], aes(x = pressure, y = var1, color = cluster, group = profile_unique)) +
            geom_line() + theme_bw())
    if(params$var > 1){
      print(ggplot(t_df[!is.na(t_df[,'pred1']),], aes(x = pressure, y = pred1, color = cluster, group = profile_unique)) +
              geom_line() + theme_bw())
      print(ggplot(t_df[!is.na(t_df[,'pred2']),], aes(x = pressure, y = pred2, color = cluster, group = profile_unique)) +
              geom_line() + theme_bw())
    }
  }
}


SW2 <- function(n, nvertex, ncolor, edges, beta, weights=  1) {
  if (ncol(edges) != 2) 
    stop("'edges' with two columns have to be provided.")
  nedge <- nrow(edges)
  if (nedge < length(weights)) 
    stop("The number of 'edges' is less than the number of 'weights'.")
  if (length(weights) < nedge) {
    weights <- rep(weights, length = nrow(edges))
  }
  bondProbs <- 1 - exp(weights * (-beta))
  oneIteration <- sample(x = 1:ncolor, nvertex, replace = TRUE) - 
    1
  oneIteration <- structure(as.integer(oneIteration), dim = dim(oneIteration))
  edges <- edges - 1
  edges <- structure(as.integer(edges), dim = dim(edges))
  colors <- .Call("sw", PACKAGE = "PottsUtils", bondProbs, oneIteration, edges, nedge, 
                  as.integer(n), as.integer(nvertex), as.integer(ncolor))
  colors + 1
}

check_permutations <- function(true, test, n = 3){
  clusters = seq(1,n, by = 1)
  true = as.factor(true)
  perms = expand.grid(clusters,clusters,clusters)
  perms = perms[apply(perms, 1, anyDuplicated) == 0, ]
  results = rep(0, nrow(perms))
  for(i in 1:nrow(perms)){
    results[i] =  sum(true == factor(test, labels = perms[i,]))/length(true)
  }
  return(max(results))
}

create_cv_splits <- function(obs_df, id, cv_num){
  ids = unique(obs_df[,id])
  start = 1
  skip = ceiling(length(ids)/cv_num)
  ids_list = lapply(1:cv_num, function(c){
    start = (c-1)*skip + 1
    ids[start:min(start+skip-1, length(ids))]
  })
  df_list = lapply(1:cv_num, function(c){
    obs_df[obs_df[,id] %in% ids_list[[c]],]
  })
  return(df_list)
}

create_train_and_test_cv <- function(df_list, id, c){
  ret = list()
  counter = 1
  for (i in 1:length(df_list)){
    if (i == c){
      next
    }
    ret[[counter]] = df_list[[i]]
    counter = counter + 1
  }
  return(list(train = do.call(rbind, ret), test = df_list[[c]]))
}