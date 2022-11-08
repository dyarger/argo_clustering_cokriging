library(Rcpp)
library(FNN)
library(Matrix)
library(fda)
library(Rfast)
library(dplyr)
library(foreach)
library(doParallel)
library(mvtnorm)
library(fstmr)

source('paper/code/src/simulation_src.R')

################################################################################
#### Simulation setup ##########################################################
################################################################################
params = list()

params[['n']] = 500
params[['var']] = 3
params[['G']] = 1
params[['n_pcs']] = c(4,4)
params[['knn']] = 10
params[['n_obs']] = 150

n_pcs = params[['n_pcs']]
range = c(0,1)
pcs = list()
n_knots = n_pcs[1]-2
knots = seq(range[1], range[2], length.out = n_knots)
basis = fda::create.bspline.basis(knots)
inprod_mat = fda::inprod(basis, basis)
inv_inprod_chol = Matrix::solve(chol(inprod_mat))
which_pc = c(1,2,3,4)
for (g in 1:params$G){
  pcs[[g]] = list()
  for (q in 1:(n_pcs[1])){
    which_pcs = which_pc[q]
    coefs = 1
    pc = list(basis=basis, ortho=inv_inprod_chol,
              which_pcs=which_pcs, coefs=coefs)
    pcs[[g]][[q]] = list(basis=basis, ortho=inv_inprod_chol,
                         which_pcs=which_pcs, coefs=coefs)
  }
}
params[['pc_funs']][['resp']] = pcs

pcs = list()
n_knots = n_pcs[2] - 2
knots = seq(range[1], range[2], length.out = n_knots)
basis = fda::create.bspline.basis(knots)
inprod_mat = fda::inprod(basis, basis)
inv_inprod_chol = Matrix::solve(chol(inprod_mat))
for (g in 1:params$G){
  pcs[[g]] = list()
  for (q in 1:(n_pcs[2])){
    which_pcs = which_pc[q]
    coefs = 1/sqrt(2)
    pc = list(basis=basis, ortho=inv_inprod_chol,
              which_pcs=which_pcs, coefs=coefs)
    pcs[[g]][[q]] = list(basis=basis, ortho=inv_inprod_chol,
                         which_pcs=which_pcs, coefs=coefs)
  }
}
params[['pc_funs']][['pred1']] = pcs

pcs = list()
for (g in 1:params$G){
  pcs[[g]] = list()
  for (q in 1:(n_pcs[2])){
    which_pcs = which_pc[n_pcs[2]-q+1]
    coefs = 1/sqrt(2)
    pc = list(basis=basis, ortho=inv_inprod_chol,
              which_pcs=which_pcs, coefs=coefs)
    pcs[[g]][[q]] = list(basis=basis, ortho=inv_inprod_chol,
                         which_pcs=which_pcs, coefs=coefs)
  }
}
params[['pc_funs']][['pred2']] = pcs

params[['var_resp']] = matrix(rep(c(1.6, 1.1, .6, .1), params$G),
                              nrow = params$G, byrow = T)
params[['st_range_resp']][[1]] = list(.1, .1, .11, .09)

params[['var_pred']][[1]] = matrix(rep(c(1.3, .5, .2, .1), params$G),
                                   nrow = params$G, byrow = T)
params[['st_range_pred']][[1]] = list(.05, .06, .05, .06)

params[['me_var']] = c(0.5, 0.3, 0.4)

params[['lt']] = list()
params[['lt']][[1]] = list()
params[['lt']][[1]][[1]] = function(p){return(-1.5*sin(5*(p+0.2)) + p)}
params[['lt']][[1]][[2]] = function(p){return(-2*cos(7*p))}
params[['lt']][[1]][[3]] = function(p){return(3*sin(-5.5*(1-p))*(1.1-p)^2)}
params[['lt']][[1]][[4]] = function(p){return(2*exp((1-p))*log(1-p+0.1) + 1)}

################################################################################
#### MCEM setup ################################################################
################################################################################

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

# Setup parallel processing
cores = min(20, detectCores())
registerDoParallel(cores)

r = list()

r= c(r, foreach(bs = 1:100, .verbose = T, .errorhandling = 'pass') %dopar% {
  
  set.seed(1000*bs)
  data = simulate_data(params)
  
  obs_df = data$data_df
  index_df = data$index_df
  params = data$params
  
  model = fstmr(obs_df, EM_params, compute_diagnostics = T, verbose = 5)
  return(model)
})

save(r, file = '~/argo_functional_regression/paper/results/mult_var_simu_pcs.RData')

