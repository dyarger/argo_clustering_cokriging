library(Rcpp)
library(FNN)
library(Matrix)
library(fda)
library(Rfast)
library(dplyr)
library(foreach)
library(doParallel)
library(scales)
library(fstmr)

# EM Setup (cf. simu_mult_var.R)
# 'Formula'
EM_params = list()
EM_params[['Y']] = 'var1'
EM_params[['X']] = c('pred1', 'pred2')
EM_params[['compute_likelihood']] = F
EM_params[['id']] = 'profile_unique'
EM_params[['ind_name']] = 'pressure'

# Markov Random Field 
EM_params[['nn_strategy']] = 'nn'
EM_params[['nn']] = 5
EM_params[['nn_type']] = 'space'
EM_params[['loc_names']] = c('latitude', 'longitude')

# Bases
EM_params[['knots']] = seq(0, 1, length.out = 148)
EM_params[['basis_response']] = fda::create.bspline.basis(EM_params[['knots']])
EM_params[['basis_predictors']] = list(
  fda::create.bspline.basis(EM_params[['knots']]),
  fda::create.bspline.basis(EM_params[['knots']]))

# Initialization 
EM_params[['domain_Y']] = c(0,1)
EM_params[['init_strategy']] = 'kmeans'
EM_params[['levels']] = list('var1' = seq(0,1,length.out = 30),
                             'pred1'= seq(0,1,length.out = 30),
                             'pred2'= seq(0,1,length.out = 30))

# Covariance estimation
EM_params[['covariance_function']] <- 'exponential_isotropic'
EM_params[['m_train_response']] <- 25
EM_params[['m_train_predictors']] <- 25

# Number of clusters and pcs
EM_params[['G']] = 1
EM_params[['pc_response1']] <- 4
EM_params[['pc_response2']] <- 4
EM_params[['pc_predictors']] <- 4

# Monte Carlo estimation
EM_params[['MC_max']] = 100
EM_params[['EM_iter']] = 60
EM_params[['EM_iter_init']] = 10
EM_params[['maxit']] = 35

# Smoothing
EM_params[['lambda_mean_response']] <- c(1, 1e-5, 1e-1, 10, 100)
EM_params[['lambda_mean_predictors']] <- c(1, 1e-5, 1e-1, 10, 100)
EM_params[['lambda_lt']] <- c(1, 1e-5, 1e-1, 10, 100)
EM_params[['lambda_pcs_response']] <- c(1, 1e-5, 1e-1, 10, 100)
EM_params[['lambda_pcs_predictors']] <- c(1, 1e-5, 1e-1, 10, 100)
EM_params[['init_strategy']] <- 'kmeans'
EM_params[['lonlat']] = F

EM_params[["m_AIC_response"]] = 25
EM_params[["m_AIC_predictors"]] = 25
EM_params[["MC_max_AIC"]] = 10
EM_params[['cv_skip']] = 10

EM_params = fstmr:::init_params(EM_params, T)

AIC_list = list()
for(i in 1:100){
  AIC_list[[i]] = matrix(1e20, 4, 4)
}

for(m in 3:6){
  load(paste0('~/argo_functional_regression/paper/results/mult_var_simu_pcs', as.character(m),'.RData'))
  r_copy = r
  for(j in 3:6){
    for(i in 1:length(r)){
      print(i)
      EM_params[['pc_response1']] <- m
      EM_params[['pc_response2']] <- j
      EM_params[['pc_predictors']] <- m
      r[[i]]$parameters$Omega_response2[[1]] = r_copy[[i]]$parameters$Omega_response2[[1]][,1:j]
      temp_m = fstmr:::compute_information_criteria(r[[i]], EM_params, r[[i]]$data_inputs, r[[i]]$parameters)
      ind = i
      AIC_list[[ind]][j-2, m-2] = temp_m$AIC
    }
  }  
}

# Check which is the correct number of PCs
rows = rep(0, length(AIC_list))
columns = rep(0, length(AIC_list))
for (i in 1:length(AIC_list)){
  a = AIC_list[[i]]
  columns[i] = ceiling(which.min(a) / nrow(a))
  rows[i] = which.min(a) - (columns[i]-1)*nrow(a) 
}
