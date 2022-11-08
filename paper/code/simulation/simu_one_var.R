array_id = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID')) 
library(fstmr)
library(dplyr)
library(mclust)
source('paper/code/src/simulation_src.R')
library(Matrix)

output_folder <- 'paper/results/simulation_one_var/'
system(paste0('mkdir ', output_folder))
# simulation settings that we might want to vary
spatial_strength <- c('weak')
n_obs <- c(20, 100)
n_profiles <- c(200)
mean_type <- c('same')
pc_type <- c('near_0')
pc_type2 <- c('different')
initialization <- c(#'10 iterations of independent model',
                    '5 iterations of independent model'#,
                    #'Random'
)
n_simu <- 100 # number of simulations per setting
n_simu_total <-  n_simu * length(spatial_strength) * length(n_obs) * length(mean_type) *
  length(pc_type) * length(pc_type2)* length(n_profiles) * length(initialization)
set.seed(55)
seeds <- sample(1:100000, n_simu_total)
# make sure each methodology can get the same data
grid_simu <- cbind('simu_unique' = 1:n_simu_total, 'seed' = seeds,
                   expand.grid('simu' = 1:n_simu, 'n_obs' = n_obs,
                               'n_profiles' = n_profiles, 
                               'spat' = spatial_strength,
                               'mean_type' = mean_type,
                               'pc_type'  = pc_type, 
                               'pc_type2'  = pc_type2, 
                               'initialization' = initialization))

# what to do for each simulation
models_to_use <- data.frame('cluster_sampling' = c('importance_sampling','gibbs',
                                                   'P_I', 'independence', 
                                                   'independence'),
                            'use_MRF' = c(rep(T, 4), F),
                            'model_unique' = 1:5)
grid_run <- expand.grid('simu_unique' = 1:n_simu_total, 'model_unique' = 1:nrow(models_to_use)) %>%
  left_join(models_to_use, by = 'model_unique')

# combine simulaations and methods
simu_setup_df <- dplyr::left_join(grid_simu, grid_run, by = 'simu_unique') 

vec_params <- simu_setup_df[array_id, ]
set.seed(vec_params[['seed']])
#### Specify params
params = list()
params[['n_obs']] <- vec_params[['n_obs']]
params[['var']] <- 1
params[['n']] <- vec_params[['n_profiles']]
params[['G']] <- 2
params[['n_pcs']] <- 3
params[['theta']] <- .5
params[['knn']] <- 5
params[['me_var']] <- 1
params[['var_resp']] <- matrix(rep(c(2, 1, 0.5)/6, params[['G']]), nrow = params[['G']], byrow = T)

# spatial range parameters
if (vec_params[['spat']] == 'weak') {
  st_range_mat <- matrix(c(.2, .14, .1)/2, nrow = params[['G']], ncol = params[['n_pcs']])
} else if (vec_params[['spat']] == 'strong') {
  st_range_mat <- matrix(c(.2, .14, .1)*2, nrow = params[['G']], ncol = params[['n_pcs']])
}
params[['st_range_resp']] <- lapply(1:params[['G']], function(g) {
  lapply(1:params[['n_pcs']], function(q) {
    st_range_mat[g,q]
  })
})

# means
params[['means']] = list()
params[['means']][[1]] = function(p) {exp(2-2*p) * cos(5*(p-0.1))}
if (vec_params[['mean_type']] == 'different') {
  params[['means']][[2]] = function(p) {3*cos(1.5*p) }
} else {
  params[['means']][[2]] = params[['means']][[1]]
}

# pcs
params[['pc_funs']] = list()
n_knots = 9
range = c(0,1)
knots = seq(range[1], range[2], length.out = n_knots)
basis = fda::create.bspline.basis(knots)
inprod_mat = inprod(basis, basis)
inv_inprod_chol = solve(chol(inprod_mat))
for (g in 1:params$G){
  params[['pc_funs']][[g]] = list()
  for (q in 1:params[['n_pcs']]){
    if (vec_params[['pc_type2']] == 'same') {
      if (vec_params[['pc_type']] == 'near_0') {
        which_pcs = 1 + q*2
      } else {
        which_pcs = 1 + q*2
      }
    } else {
      if (vec_params[['pc_type']] == 'near_0') {
        which_pcs = g + q*2 
      } else {
        which_pcs = g + q*2
      }
    }
    params[['pc_funs']][[g]][[q]] = list(basis=basis, ortho=inv_inprod_chol, which_pcs=which_pcs, coefs=1)
  }
}

# plot PC functions
if (array_id == 1) {
  p <- seq(0, 1, length.out = 200)
  setup_df <- data.frame(p, 'mean.1' = params[['means']][[1]](p),
                         'mean.2' = params[['means']][[2]](p),
                         'pc.1_1' = evaluate_pcs(p, params[['pc_funs']][[1]][[1]]),
                         'pc.1_2' = evaluate_pcs(p, params[['pc_funs']][[1]][[2]]),
                         'pc.1_3' = evaluate_pcs(p, params[['pc_funs']][[1]][[3]]),
                         'pc.2_1' = evaluate_pcs(p, params[['pc_funs']][[2]][[1]]),
                         'pc.2_2' = evaluate_pcs(p, params[['pc_funs']][[2]][[2]]),
                         'pc.2_3' = evaluate_pcs(p, params[['pc_funs']][[2]][[3]]))
  
  
  setup_long <- tidyr::pivot_longer(setup_df, cols = starts_with('pc'), #names_to = 'Cluster',
                                    names_to = c('Cluster', 'PC'),
                                    names_sep = '_',
                                    names_prefix = 'pc.', values_to = 'value')
  
  ggplot(data = setup_long, aes(x = p, y = value, group = paste0(Cluster, PC), color = Cluster,
                                linetype = Cluster))+
    geom_line() +
    labs(x  = 't', y = 'Principal Component Functions') + 
    facet_wrap(~PC, nrow = 1, labeller = label_both) + 
    theme_bw() + 
    theme(legend.position = 'bottom')
  
  ggsave(filename = 'paper/images/one_var_sim_setup.png', height = 3.5/1.1, width = 11/1.1)
}
data = simulate_data(params)

obs_df = data[['data_df']]
index_df = data[['index_df']]
params = data[['params']]

# now, we set up our estimation strategy
EM_params = list()
EM_params[['Y']] = 'var1'
EM_params[['id']] = 'profile_unique'
EM_params[['ind_name']] = 'pressure'
EM_params[['nn']] = 5
EM_params[['nn_type']] = 'space'
EM_params[['nn_strategy']] = 'nn'
EM_params[['loc_names']] = c('longitude', 'latitude')
EM_params[['domain_Y']] = c(0,1)
EM_params[['levels']] = list('var1' = seq(0,1,length.out = 15))
EM_params[['pc_response']] = 3
EM_params[['lambda_pcs_response']] = 10^c(-4:0)
EM_params[['lambda_mean_response']] <- 10^c(-4:0)
# EM_params[['lambda_pcs_response']] = 10^-2
# EM_params[['lambda_mean_response']] <- 10^-2
EM_params[['G']] = 2

EM_params[['MC_max']] = 30
EM_params[['EM_iter']] = 60
EM_params[['EM_iter_init']] = 5
EM_params[['maxit']] = 20
EM_params[['init_strategy']] = 'random' # works better when clusters have same means
knots = seq(0,1, length.out = 20)
EM_params[['basis_response']] = create.bspline.basis(knots)
EM_params[['m_AIC_response']] = NA
EM_params[['m_AIC_predictors']] = NA
EM_params[['MC_max_AIC']] = NA

# run independence model for a few iterations
percent_correct <- function(x, y) {
  max(c(mean(x == y),
        mean(x != y)))
}
EM_params[['cv_skip']] <- 4
EM_params[['cluster_sampling']] = 'independence'
EM_params[['use_MRF']] = vec_params[['use_MRF']]
data_inputs = fstmr:::prepare_EM_data(obs_df, EM_params)
model = fstmr:::determine_model(params)
parameters = fstmr:::initialize.one_var(model, EM_params, data_inputs)
parameters[['lambda_mean_response']] <- parameters[['lambda_mean_response']][1]
parameters[['lambda_pcs_response']] <- parameters[['lambda_pcs_response']][1]
ARI_ind <- perc_correct_ind <- theta_vals_ind <-  vector(length = EM_params[['EM_iter_init']] + 1, 'numeric')
perc_correct_ind[1] <-  percent_correct(parameters[['cluster_mat']][,1], index_df[['cluster']])
ARI_ind[1] <-  mclust::adjustedRandIndex(parameters[['cluster_mat']][,1], index_df[['cluster']])
theta_vals_ind[1] <- parameters[['theta']]
EM_params[['MC_max_AIC']] <- 15
EM_params[['m_AIC_response']] <- 15
for (i in 1:EM_params[['EM_iter_init']]) {
  MCEM_res = fstmr:::Estep.one_var(model, EM_params, data_inputs, parameters)
  parameters = fstmr:::Mstep.one_var(model, EM_params, data_inputs, parameters, 
                                     MCEM_res, iter = i)
  perc_correct_ind[i+1] <- percent_correct(parameters[['cluster_membership']], index_df[['cluster']])
  ARI_ind[i+1] <-  mclust::adjustedRandIndex(parameters[['cluster_membership']], index_df[['cluster']])
  theta_vals_ind[i+1] <-  parameters[['theta']]
}
EM_params[['cluster_sampling']] = vec_params[['cluster_sampling']]
perc_correct_ind
# 4 spatial covariance parameters
if (EM_params[['cluster_sampling']] != 'independence') {
  EM_params[['covariance_function']] <- 'exponential_isotropic'
  EM_params[['m_train_response']] = 15
}
# p <- seq(0, 1, by = .01)
# Phi <- eval.basis(seq(0, 1, by = .01), EM_params$basis_response)
# plot(p, Phi %*% parameters$means_response[[1]], type = 'l')
# plot(p, Phi %*% parameters$Omega_response[[2]][,3], type = 'l')

perc_correct <- ARI <- theta_vals <-  vector(length = EM_params[['EM_iter']], 'numeric')
ess <- vector(length = EM_params[['EM_iter']], 'numeric')
for (i in 1:EM_params[["EM_iter"]]) {
  MCEM_res = fstmr:::Estep.one_var(model, EM_params, data_inputs, parameters)
  parameters = fstmr:::Mstep.one_var(model, EM_params, data_inputs, parameters, 
                                     MCEM_res, iter = i)
  perc_correct[i] <- percent_correct(parameters[['cluster_membership']], index_df[['cluster']])
  ess[i] <- sum(MCEM_res[['weights']])^2/sum(MCEM_res[['weights']]^2)
  ARI[i] <-  mclust::adjustedRandIndex(parameters[['cluster_membership']], index_df[['cluster']])
  theta_vals[i] <- parameters[['theta']]
  print(perc_correct)
}

IC <- fstmr:::compute_information_criteria(model, EM_params, data_inputs, parameters)

simu_results <- cbind('array_id' = array_id, simu_setup_df[array_id, ], data.frame('iter' = 0:(EM_params[['EM_iter_init']] + EM_params[['EM_iter']]),
                                                                                   correct = c(perc_correct_ind, perc_correct),
                                                                                   ess = c(rep(NA, EM_params$EM_iter_init + 1),
                                                                                           ess),
                                                                                   ARI = c(ARI_ind, ARI),
                                                                                   theta = c(theta_vals_ind, theta_vals)))
save(simu_results, IC,
     file = paste0(output_folder, 'simu_', array_id, '.RData'))



