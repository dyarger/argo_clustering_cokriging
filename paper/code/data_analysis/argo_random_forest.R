Sys.setenv(R_LIBS_SITE = '/sw/arcts/centos7/stacks/gcc/8.2.0/Rtidyverse/4.1.0/2021-06-28')
library(fstmr)
library(dplyr)
library(fda)
library(randomForest)
# 1 variables and basic setup
set.seed(56)
params = list()
params[['Y']] = Sys.getenv('var')
params[['results_folder']] <- "/scratch/stats_dept_root/stats_dept1/shared_data/argo/"
# 2 type of leave out data
params[['leave_out_type']] = Sys.getenv('leave_out_type')
params[['n_folds']] = as.integer(Sys.getenv('n_folds'))
params[['folds']] <-  as.integer(Sys.getenv('fold'))
params[['leave_out_prediction_core']] = F

##### Load and subset data #####
load(paste0('paper/data/soccom_processed_', params[['Y']], '_05_05_21.RData'))
source('paper/code/src/prepare_argo_data_src.R')
data_leave_out <- leave_out_data(df_list, params[['leave_out_type']], 
                                 params[['folds']], params[['n_folds']], params[['leave_out_prediction_core']],
                                 0, params[['core_sampling_strategy']],
                                 reduce_data = F
)
df_list <- data_leave_out[[1]]; df_list_pred <- data_leave_out[[2]]
held_out_profs <- df_list_pred[[1]]; prediction_profs <- df_list_pred[[2]]
profiles_for_prediction <- data_leave_out[[3]]
rm(data_leave_out)
df_unique <- df_list[[2]][!duplicated(df_list[[2]][['profile_unique']]),]

# random forest
p = 6 # T, S, lat, long, year, month
m = floor(p/3)
B = 500

df_train_150 <- df_list[[1]] %>%
  filter(pressure < 155, pressure > 145) %>%
  mutate(year = as.double(substr(day, 1, 4)),
         month = as.double(substr(day, 6, 7)))
df_test_150 <- df_list_pred[[2]] %>%
  filter(pressure < 155, pressure > 145) %>%
  mutate(year = as.double(substr(day, 1, 4)),
         month = as.double(substr(day, 6, 7)))
df_test_150 <- df_test_150[!is.na(df_test_150[[params[['Y']]]]),]

params[['X']] = c('temp', 'psal', 'longitude', 'latitude', 'year', 'month')

X <- df_train_150[, params[['X']]]
X_out <- df_test_150[, params[['X']]]
rf_limit <- randomForest::randomForest(y = df_train_150[[params[['Y']]]], x = X, ntree = B,
                                       mtry = m)
pred_limit <- predict(rf_limit, X_out)

print(median(abs(pred_limit - df_test_150$oxy)))
print(sqrt(mean((pred_limit - df_test_150$oxy)^2)))

rf_results <- data.frame(fold = params[['folds']], 
                         profile_unique = df_test_150$profile_unique,
                         pressure = df_test_150$pressure, 
                         oxy = df_test_150$oxy,
                         residual = df_test_150$oxy - pred_limit)
system('mkdir paper/results/randomForest/')
file_name <- paste0('paper/results/randomForest/', params[['leave_out_type']], '_', params[['Y']], '_',
                    params[['n_folds']], '_', params[['folds']], '_rf.RData')
save(rf_results, file = file_name)
