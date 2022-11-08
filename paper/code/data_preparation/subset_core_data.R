# now we look at core data
library(tidyverse)
soccom_date <- Sys.getenv('soccom_date')
argo_date <- Sys.getenv('argo_date')
load(paste0('paper/data/soccom_all_variables_', soccom_date, '.RData'))
locs <- soccom %>%
  mutate(profile_unique = paste0(float, '_', profile)) %>%
  filter(!duplicated(profile_unique)) %>%
  dplyr::select(longitude, latitude) %>%
  as.matrix()
BGC_floats <- unique(soccom[['float']])
core_files <- list.files('paper/data/core_profiles/')
core_data <- list()
for (i in 1:length(core_files)) {
  load(paste0('paper/data/core_profiles/', core_files[i]))
  head(single_prof)
  pressure <- single_prof[['pressure']]
  psal <- single_prof[['salinity']]
  
  # do not include profiles if certain conditions were not met
  if (nrow(single_prof) < 15 | max(pressure) < 1450 |
      max(abs(lead(pressure) - pressure), na.rm = T) > 200 | 
      min(pressure) > 100 | 
      max(psal) > 37.2 | min(psal) < 33 | 
      single_prof[['float']][1] %in% BGC_floats | 
      sum(single_prof[['salinity_qc']] %in% c(3,4) > 0) | 
      sum(single_prof[['temperature_qc']] %in% c(3,4) > 0) | 
      sum(single_prof[['pressure_qc']] %in% c(3,4) > 0) | 
      max(single_prof[['pres_adj_error']]) > 16 | 
      is.na(single_prof[['lat']][1]) ) {
    next
  }
  distances_core <- min(fields::rdist.earth.vec(locs,
                                                single_prof[1, c('long', 'lat')],
                                                miles = F), na.rm = T)
  if (distances_core > 500 & single_prof$lat[1] > -55){
    next
  }
  core_data[[i]] <- single_prof %>%
    filter(pressure < 2000, pressure > 0) %>%
    rename(psal = salinity, temp = temperature,
           longitude = long, latitude = lat) %>%
    mutate(profile_unique = paste(float, cycle, sep = '_')) %>%
    dplyr::select(pressure, temp, psal, profile_unique, longitude, latitude, 
                  day, date, date_time)
  if (i %% 1000 == 0) {
    print(i)
    gc()
  }
}
core_data <- dplyr::bind_rows(core_data) %>%
  mutate(longitude = ifelse(longitude > 180, longitude - 360, longitude)) %>%
  mutate(dayofyear = julian(date, origin = as.Date('2000-01-01')) %% 365.25)
save(core_data, file = 'paper/data/core_processed_', argo_date, '.RData')
