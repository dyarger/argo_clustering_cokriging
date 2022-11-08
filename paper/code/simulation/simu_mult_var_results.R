library(Rcpp)
library(FNN)
library(Matrix)
library(fda)
library(Rfast)
library(dplyr)
library(foreach)
library(doParallel)
library(mvtnorm)
library(scales)
library(fstmr)
library(ggplot2)
library(dplyr)
library(tikzDevice)

setwd('/home/kortest/')
source('argo_functional_regression/paper/code/src/simulation_src.R')

###################################################
# Simulation parameters for reference
###################################################
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
params = init_simu_params(params)

#############################################
# Generate plots
#############################################

load('~/argo_functional_regression/paper/results/mult_var_simu.RData')
#load('~/r_new6.RData')
r = lapply(r, function(x)x$parameters)

#ps = seq(0, 1, length.out = 100)
ps = seq(0,1,length.out = 148)
basis_response = fda::create.bspline.basis(ps)
mat = fda::eval.basis(ps, basis_response)

mean_df = data.frame(type = 'True', ps = ps, val = params$means$resp[[1]](ps), iter=-1)
pc_df = data.frame(type = 'True', ps = ps, val = evaluate_pcs(ps, params$pc_funs$resp[[1]][[1]]), pc = 1, iter = -1)
lt_df = data.frame(type = 'True', ps = ps, val = params$lt[[1]][[1]](ps), pc = 1, iter = -1)

mean_df$iter = as.factor(mean_df$iter)
pc_df$iter = as.factor(pc_df$iter)
pc_df$pc = as.factor(pc_df$pc)

for (q in 2:params$n_pcs[1]){
  pc_df = rbind(pc_df, data.frame(type = 'True',
                                  ps = ps,
                                  val = evaluate_pcs(ps, params$pc_funs$resp[[1]][[q]]),
                                  pc = as.factor(q),
                                  iter = as.factor(-1)))
}
for (q in 2:(params$n_pcs[2])){
  lt_df = rbind(lt_df, data.frame(type = 'True',
                                  ps = ps,
                                  val = params$lt[[1]][[q]](ps),
                                  pc = as.factor(q),
                                  iter = as.factor(-1)))
}
c = 1
counter = 0
for(p in r){
  print(is.null(p$measurement_error_response))
  if (is.null(p$measurement_error_response)){
    counter = counter + 1
    next
  }
  mean_df = rbind(mean_df, data.frame(type = 'Estimates',
                                      ps = ps,
                                      val = mat %*% p$means_response[[1]],
                                      iter = as.factor(c)))
  for (q in 1:params$n_pcs[1]){
    pc_df = rbind(pc_df, data.frame(type = 'Estimates',
                                    ps = ps,
                                    val = -mat %*% p$Omega_response2[[1]][,q],
                                    pc = as.factor(q),
                                    iter = as.factor(c)))
  }
  for (q in 1:(params$n_pcs[2])){
    lt_df = rbind(lt_df, data.frame(type = 'Estimates',
                                    ps = ps,
                                    val = mat %*% p$Omega_response1[[1]][,q],
                                    pc = as.factor(q),
                                    iter = as.factor(c)))
  }
  c = c+1
}

gdf = split(pc_df, pc_df$iter)
truth = gdf[[1]]
ggdf = split(truth, truth$pc)
truths = list()
for (q in 1:params$n_pcs[1]){
  truths[[q]] = ggdf[[q]]$val
}
pc_df2 = gdf[[1]]
for (x in gdf){
  ggdf = split(x, x$pc)
  for (q in 1:params$n_pcs[1]){
    if (mean((ggdf[[q]]$val - truths[[q]])^2) > mean((-(ggdf[[q]]$val) - truths[[q]])^2)){
      ggdf[[q]]$val = - (ggdf[[q]]$val)
    }
    pc_df2 = rbind(pc_df2, ggdf[[q]])
  }
}

gdf = split(lt_df, lt_df$iter)
truth = gdf[[1]]
ggdf = split(truth, truth$pc)
truths = list()
for (q in 1:(params$n_pcs[2])){
  truths[[q]] = ggdf[[q]]$val
}
lt_df2 = gdf[[1]]
for (x in gdf){
  ggdf = split(x, x$pc)
  for (q in 1:(params$n_pcs[2])){
    if (mean((ggdf[[q]]$val - truths[[q]])^2) > mean((-(ggdf[[q]]$val) - truths[[q]])^2)){
      ggdf[[q]]$val = - (ggdf[[q]]$val)
    }
    lt_df2 = rbind(lt_df2, ggdf[[q]])
  }
}

pc_df21 = pc_df2 %>%
  filter(pc == 1)

pc_df22 = pc_df2 %>%
  filter(pc == 2)

pc_df23 = pc_df2 %>%
  filter(pc == 3)

pc_df24 = pc_df2 %>%
  filter(pc == 4)



lt_df21 = lt_df2 %>%
  filter(pc == 1)

lt_df22 = lt_df2 %>%
  filter(pc == 2)

lt_df23 = lt_df2 %>%
  filter(pc == 3)

lt_df24 = lt_df2 %>%
  filter(pc == 4)


for (j in 1:4){
  tr = "True"
  data_df = paste0('pc_df2', as.character(j), ' %>%
  filter(type != tr)')
  plot_df1 = eval(parse(text=data_df))
  gcdf = split(plot_df1, plot_df1$iter)
  
  res = matrix(0, ncol=nrow(gcdf[[plot_df1$iter[1]]]), nrow = length(gcdf))
  counter = 1
  for(i in unique(plot_df1$iter)){
    res[counter,] = gcdf[[i]]$val 
    counter = counter +1 
  }
  data_df = paste0('pc_df2', as.character(j), ' %>%
  filter(type == tr)')
  truthc = eval(parse(text=data_df))
  truthc = truthc$val[1:148]
  plot_df = data.frame(p=seq(0,1, length.out = 148), val = truthc, type = 'True')
  plot_df = rbind(plot_df, data.frame(p=seq(0,1, length.out = 148), val = colMeans(res), type = 'Simulation Average'))
  plot_df = rbind(plot_df, data.frame(p=seq(0,1, length.out = 148), val = colMeans(res) + 2*sqrt(colVars(res)),type = 'Upper'))
  plot_df = rbind(plot_df, data.frame(p=seq(0,1, length.out = 148), val = colMeans(res) - 2*sqrt(colVars(res)),type = 'Lower'))
  
  filename = paste0("argo_functional_regression/paper/images/pc",as.character(j), ".tex")
  guidepose = c(0.7, 0.1)
  if (j == 1){
    guidepose = c(0.7, 0.9)
  }
  tikz(file = filename)
  lbsize = 20
  titlesize = 20
  legsize = 20
  abc = ggplot(data=plot_df) +
    geom_line(data = plot_df %>% filter(type %in% c("True", "Simulation Average")), aes(x = p, y = val, color = type)) +
    #geom_line(data = plot_df %>% filter(type == "Estimate"), aes(x = p, y = val, color = 'blue')) +
    geom_ribbon(data = plot_df %>% filter(type == "Simulation Average"),
                aes(x = p, ymin=(plot_df %>% filter(type == "Lower"))$val,
                    ymax=(plot_df %>% filter(type == "Upper"))$val),
                fill="grey70", alpha=0.3) +
    xlab("FindMeXlab") + ylab("FindMeYlab") + 
    theme_bw() + 
    scale_y_continuous(breaks= pretty_breaks()) +
    theme(legend.position = guidepose, legend.text=element_text(size=legsize), legend.key.size = unit(1.8, 'lines')) + 
    theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(axis.text.x = element_text(size = lbsize), axis.text.y = element_text(size = lbsize)) + 
    theme(axis.title.y = element_text(size = titlesize), axis.title.x = element_blank()) 
  print(abc)
  dev.off()  
}


for (j in 1:4){
  tr = "True"
  data_df = paste0('lt_df2', as.character(j), ' %>%
  filter(type != tr)')
  plot_df1 = eval(parse(text=data_df))
  gcdf = split(plot_df1, plot_df1$iter)
  
  res = matrix(0, ncol=nrow(gcdf[[plot_df1$iter[1]]]), nrow = length(gcdf))
  counter = 1
  for(i in unique(plot_df1$iter)){
    res[counter,] = gcdf[[i]]$val 
    counter = counter +1 
  }
  data_df = paste0('lt_df2', as.character(j), ' %>%
  filter(type == tr)')
  truthc = eval(parse(text=data_df))
  truthc = truthc$val[1:148]
  plot_df = data.frame(p=seq(0,1, length.out = 148), val = truthc, type = 'True')
  plot_df = rbind(plot_df, data.frame(p=seq(0,1, length.out = 148), val = colMeans(res), type = 'Simulation Average'))
  plot_df = rbind(plot_df, data.frame(p=seq(0,1, length.out = 148), val = colMeans(res) + 2*sqrt(colVars(res)),type = 'Upper'))
  plot_df = rbind(plot_df, data.frame(p=seq(0,1, length.out = 148), val = colMeans(res) - 2*sqrt(colVars(res)),type = 'Lower'))
  
  filename = paste0("argo_functional_regression/paper/images/lt",as.character(j), ".tex")
  tikz(file = filename)
  guidepose = c(0.7, 0.1)
  if (j > 2){
    guidepose = c(0.7, 0.9)
  }
  tikz(file = filename)
  lbsize = 20
  titlesize = 20
  legsize = 20
  abc = ggplot(data=plot_df) +
    geom_line(data = plot_df %>% filter(type %in% c("True", "Simulation Average")), aes(x = p, y = val, color = type)) +
    #geom_line(data = plot_df %>% filter(type == "Estimate"), aes(x = p, y = val, color = 'blue')) +
    geom_ribbon(data = plot_df %>% filter(type == "Simulation Average"),
                aes(x = p, ymin=(plot_df %>% filter(type == "Lower"))$val,
                    ymax=(plot_df %>% filter(type == "Upper"))$val),
                fill="grey70", alpha=0.3) +
    xlab("FindMeXlab") + ylab("FindMeYlab") + 
    theme_bw() + 
    scale_y_continuous(breaks= pretty_breaks()) +
    theme(legend.position = guidepose, legend.text=element_text(size=legsize), legend.key.size = unit(1.8, 'lines')) + 
    theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(axis.text.x = element_text(size = lbsize), axis.text.y = element_text(size = lbsize)) + 
    theme(axis.title.y = element_text(size = titlesize), axis.title.x = element_blank()) 
  print(abc)
  dev.off()  
}

tr = "True"
plot_df1 = mean_df %>%
  filter(type != tr)
gcdf = split(plot_df1, plot_df1$iter)

res = matrix(0, ncol=nrow(gcdf[[plot_df1$iter[1]]]), nrow = length(gcdf))
counter = 1
for(i in unique(plot_df1$iter)){
  res[counter,] = gcdf[[i]]$val 
  counter = counter +1 
}
truthc = mean_df %>%
  filter(type == tr)
truthc = truthc$val[1:148]
plot_df = data.frame(p=seq(0,1, length.out = 148), val = truthc, type = 'True')
plot_df = rbind(plot_df, data.frame(p=seq(0,1, length.out = 148), val = colMeans(res), type = 'Simulation Average'))
plot_df = rbind(plot_df, data.frame(p=seq(0,1, length.out = 148), val = colMeans(res) + 2*sqrt(colVars(res)),type = 'Upper'))
plot_df = rbind(plot_df, data.frame(p=seq(0,1, length.out = 148), val = colMeans(res) - 2*sqrt(colVars(res)),type = 'Lower'))

filename = paste0("argo_functional_regression/paper/images/mean.tex")
guidepose = c(0.7, 0.9)
tikz(file = filename)
lbsize = 20
titlesize = 20
legsize = 20
abc = ggplot(data=plot_df) +
  geom_line(data = plot_df %>% filter(type %in% c("True", "Simulation Average")), aes(x = p, y = val, color = type)) +
  #geom_line(data = plot_df %>% filter(type == "Estimate"), aes(x = p, y = val, color = 'blue')) +
  geom_ribbon(data = plot_df %>% filter(type == "Simulation Average"),
              aes(x = p, ymin=(plot_df %>% filter(type == "Lower"))$val,
                  ymax=(plot_df %>% filter(type == "Upper"))$val),
              fill="grey70", alpha=0.3) +
  xlab("FindMeXlab") + ylab("FindMeYlab") + 
  theme_bw() + 
  scale_y_continuous(breaks= pretty_breaks()) +
  theme(legend.position = guidepose, legend.text=element_text(size=legsize), legend.key.size = unit(1.8, 'lines')) + 
  theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.text.x = element_text(size = lbsize), axis.text.y = element_text(size = lbsize)) + 
  theme(axis.title.y = element_text(size = titlesize), axis.title.x = element_blank()) 
print(abc)
dev.off()  



## Generate table 1

load('~/argo_functional_regression/paper/results/mult_var_simu.RData')
r = lapply(r, function(x)x$parameters)

res = matrix(0, 1, 4)
res_var = matrix(0, 1, 4)
res_r = matrix(0, 1, 4)

res_p_var = matrix(0, 1, 4)
res_p_r = matrix(0, 1, 4)
res_p = matrix(0, 1, 4)
c = 1
for (p in r){
  if (max(sort(1 / p$range_params_r[[1]][1,] * p$variances_response[[1]])) > 100){
    next
  }
  res = rbind(res, 1 / p$range_params_r[[1]][1,] * p$variances_response[[1]])
  res_var = rbind(res_var, p$variances_response[[1]])
  res_r = rbind(res_r,p$range_params_r[[1]][1,])
  c = c + 1
  res_p = rbind(res_p, 1 / p$range_params_p[[1]][1,] * p$variances_predictors[[1]])
  res_p_var = rbind(res_p_var, p$variances_predictors[[1]])
  res_p_r = rbind(res_p_r, p$range_params_p[[1]][1,])
}

res = res[2:nrow(res),]
res_p = res_p[2:nrow(res_p),]
res_var = res_var[2:nrow(res_var),]
res_r = res_r[2:nrow(res_r),]
res_p_var = res_p_var[2:nrow(res_p_var),]
res_p_r = res_p_r[2:nrow(res_p_r),]

var_est = round(colMeans(res_var), digits = 2)
standd = round(sqrt(colVars(res_var)), digits = 2)
for(i in 1:length(var_est)){
  var_est[i] = paste0(as.character(var_est[i]), ' (', as.character(standd[i]),')')
}

range_est = round(colMeans(res_r), digits = 2)
standd = round(sqrt(colVars(res_r)), digits = 2)
for(i in 1:length(range_est)){
  range_est[i] = paste0(as.character(range_est[i]), ' (', as.character(standd[i]),')')
}

ratio_est = round(colMeans(res), digits = 2)
standd = round(sqrt(colVars(res)), digits = 2)
for(i in 1:length(range_est)){
  ratio_est[i] = paste0(as.character(ratio_est[i]), ' (', as.character(standd[i]),')')
}

results_table_resp = data.frame(var = round(as.numeric(params$var_resp), digits = 2),
                                var_est = var_est,
                                range = round(as.numeric(params$st_range_resp[[1]]), digits = 2),
                                range_est = range_est,
                                ratio = round(as.numeric(params$var_resp/as.numeric(params$st_range_resp[[1]])), digits =2),
                                ratio_est = ratio_est)
stargazer(results_table_resp, summary = FALSE, out = 'argo_functional_regression/paper/tables/table1.tex', digits=2)


var_est = round(colMeans(res_p_var), digits = 2)
standd = round(sqrt(colVars(res_p_var)), digits = 2)
for(i in 1:length(var_est)){
  var_est[i] = paste0(as.character(var_est[i]), ' (', as.character(standd[i]),')')
}

range_est = round(colMeans(res_p_r), digits = 2)
standd = round(sqrt(colVars(res_p_r)), digits = 2)
for(i in 1:length(range_est)){
  range_est[i] = paste0(as.character(range_est[i]), ' (', as.character(standd[i]),')')
}

ratio_est = round(colMeans(res_p), digits = 2)
standd = round(sqrt(colVars(res_p)), digits = 2)
for(i in 1:length(range_est)){
  ratio_est[i] = paste0(as.character(ratio_est[i]), ' (', as.character(standd[i]),')')
}

results_table_pred = data.frame(var = as.numeric(params$var_pred[[1]]),
                                var_est = var_est,
                                range = as.numeric(params$st_range_pred[[1]]),
                                range_est = range_est,
                                ratio = as.numeric(params$var_pred[[1]]/as.numeric(params$st_range_pred[[1]])),
                                ratio_est = ratio_est)

stargazer(results_table_pred, summary = FALSE, out = 'argo_functional_regression/paper/tables/table2.tex', digits=2)

