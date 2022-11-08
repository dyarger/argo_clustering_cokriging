
# df_list: data input
# leave_out_type: profiles, floats, none
# n_folds: number of folds
# folds: which fold
# leave_out_prediction_core: in leave out experiments, should temperature and salinity be left out?
# n_core_profiles: if you want, include core profiles
# core_sampling_strategy: only relevant if n_core_profiles > 0
# reduce_data: do everything with only 1000 profiles for testing
leave_out_data <- function(df_list, leave_out_type, folds, n_folds, leave_out_prediction_core,
                           n_core_profiles = 0, core_sampling_strategy, reduce_data = F) {
  # for testing look at 1000 probiles
  if (reduce_data) {
    profile_vec <- df_list[[1]][['profile_unique']];BGC_profiles <- unique(profile_vec);folds <- 3
    df_list[[1]] <- df_list[[1]][df_list[[1]][['profile_unique']] %in% BGC_profiles[1:1000],]
    df_list[[2]] <- df_list[[2]][df_list[[2]][['profile_unique']] %in% BGC_profiles[1:1000],]
  }
  if (leave_out_type == 'profiles') {
    LO_vec <- df_list[[1]][['profile_unique']]
    LO_vec_pred <- df_list[[2]][['profile_unique']]
  } else if (leave_out_type == 'floats') {
    LO_vec <- sapply(strsplit(df_list[[1]][['profile_unique']], '_'), function(x) x[1])
    LO_vec_pred <- sapply(strsplit(df_list[[2]][['profile_unique']], '_'), function(x) x[1])
  } else if (leave_out_type == 'none') {
    folds <- 2
    LO_vec <- 1
    LO_vec_pred <- 1
    n_folds <- 1
  }
  LO_unique <- unique(LO_vec)
  set.seed(50)
  held_out_groups <- sample(1:n_folds, length(LO_unique), replace = T)
  held_out_vals <- LO_unique[held_out_groups == folds]
  held_out_profs <- df_list[[1]][LO_vec %in% held_out_vals,]
  df_list[[1]] <- df_list[[1]][!(LO_vec %in% held_out_vals),]
  
  prediction_profs <- df_list[[2]][LO_vec_pred %in% held_out_vals,]

  if (leave_out_prediction_core) {
    df_list[[2]] <- df_list[[2]][!(LO_vec_pred %in% held_out_vals),]
  }
  
  # add core profiles
  if (n_core_profiles > 0) {
    load('paper/data/core_processed_06_21.RData')
    core_profiles <- unique(core_data[['profile_unique']])
    if (n_core_profiles == 'all') {
      core_profiles_use <- core_profiles
      if (core_sampling_strategy == 'south') {
        # take all core profiles south of -55, samples from other latitudes
        core_data_not_duplicated <- core_data[!duplicated(core_data[['profile_unique']]),]
        core_data_south <- core_data_not_duplicated[core_data_not_duplicated[['latitude']] < -50,]
        core_data_mid <- core_data_not_duplicated[core_data_not_duplicated[['latitude']] >= -50 & 
                                                    core_data_not_duplicated[['latitude']] < -45,]
        core_data_north <- core_data_not_duplicated[core_data_not_duplicated[['latitude']] >= -45,]
        core_profiles_use <- 
          c(core_data_south[['profile_unique']], 
            sample(core_data_mid[['profile_unique']], size = 11000),
            sample(core_data_north[['profile_unique']], size = 2000)
          )
      }
    } else {
      core_profiles_use <- sample(core_profiles, n_core_profiles, replace = F)
    }
    core_indexes <- core_data[['profile_unique']] %in% core_profiles_use
    core_data[['day']] <- core_data[['date']]
    df_list[[2]] <- rbind(cbind(df_list[[2]][,c('profile_unique', 'latitude', 'longitude', 'day', 'dayofyear', 'pressure', 'temp', 'psal')],
                                type = 'BGC'), 
                          cbind(core_data[core_indexes,c('profile_unique', 'latitude', 'longitude', 'day', 'dayofyear', 'pressure', 'temp', 'psal')],
                                type = 'core'))
  }
  df_list_pred <- list(held_out_profs %>% arrange(profile_unique, pressure),
                       prediction_profs  %>% arrange(profile_unique, pressure))
  list(df_list, df_list_pred, names(split(df_list_pred[[2]], df_list_pred[[2]]$profile_unique)))
}
