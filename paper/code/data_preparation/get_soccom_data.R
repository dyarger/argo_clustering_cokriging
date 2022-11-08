library(ncdf4)
file_folder <- Sys.getenv('soccom_location')
files <- list.files(file_folder, pattern = '.nc')
data_list <- list()
for (j in 1:length(files)) {
  file <- nc_open(paste0(file_folder, files[j]))
  
  # core info
  lat <- ncvar_get(file, 'Lat')
  lon <- ncvar_get(file, 'Lon')
  lon_qc <- ncvar_get(file, 'Lat_QFA')
  temp <- ncvar_get(file, 'Temperature')
  psal <- ncvar_get(file, 'Salinity')
  pressure <- ncvar_get(file, 'Pressure')
  pressure_QC <- ncvar_get(file, 'Pressure_QFA')
  temp_QC <- ncvar_get(file, 'Temperature_QFA')
  psal_QC <- ncvar_get(file, 'Salinity_QFA')
  day <- ncvar_get(file, 'mon_day_yr')
  time <- ncvar_get(file, 'hh_mm')
  
  # BGC info
  oxy <- ncvar_get(file, 'Oxygen')
  oxy_sat <- ncvar_get(file, 'OxygenSat') # Calculation assumes atmospheric pressure= 1013.25 mbar
  oxy_QC <- ncvar_get(file, 'Oxygen_QFA')
  oxy_sat_QC <- ncvar_get(file, 'OxygenSat_QFA')
  pdens <- ncvar_get(file, 'Sigma_theta')
  pdens_QC <- ncvar_get(file, 'Sigma_theta_QFA')
  depth <- ncvar_get(file, 'Depth')
  nitrate <- ncvar_get(file, 'Nitrate')
  nitrate_QC <- ncvar_get(file, 'Nitrate_QFA')
  pH <- ncvar_get(file, 'pHinsitu')
  pH_QC <- ncvar_get(file, 'pHinsitu_QFA')
  chl <- ncvar_get(file, 'Chl_a')
  chl_QC <- ncvar_get(file, 'Chl_a_QFA')
  poc <- ncvar_get(file, 'POC')
  poc_QC <- ncvar_get(file, 'POC_QFA')
  
  nprof <- ncol(pressure)
  npres <- nrow(pressure)
  
  if (length(dim(pressure)) == 1) {
    repeated_vars <- data.frame('float' =as.numeric(ncvar_get(file, 'Cruise')),
                                'profile' =  as.numeric(ncvar_get(file, 'Station'))-1,
                                'latitude' =                    rep(lat, each = npres),
                                'longitude'=                rep(lon, each = npres),
                                'longitude_QC' =                 rep(lon_qc, each = npres),
                                'day'=                rep(day, each = npres),
                                'time' = rep(time, each = npres))
  } else {
    repeated_vars <- data.frame('float' = rep(as.numeric(ncvar_get(file, 'Cruise')), times = npres * nprof),
                                'profile' =  rep(as.numeric(ncvar_get(file, 'Station')-1), each = npres),
                                'latitude' =                    rep(lat, each = npres),
                                'longitude'=                rep(lon, each = npres),
                                'longitude_QC' =                 rep(lon_qc, each = npres),
                                'day'=                rep(day, each = npres),
                                'time' = rep(time, each = npres))
  }
  
  mat_test <- data.frame(
    repeated_vars,
    # core
    'pressure' = as.double(pressure),
    'temp' = as.double(temp),
    'psal' = as.double(psal),
    'pdens' = as.double(pdens), 
    #BGC
    'oxy' = as.double(oxy),
    'oxy_sat' = as.double(oxy_sat),
    'nitrate' = as.double(nitrate),
    'pH' = as.double(pH),
    'chl' = as.double(chl),
    'poc' = as.double(poc),
    # Quality flags
    'pressure_QC' = as.double(pressure_QC), 'temp_QC' = as.double(temp_QC),
    'psal_QC' = as.double(psal_QC), 'oxy_QC' = as.double(oxy_QC),
    'pdens_QC' = as.double(pdens_QC), 
    'oxy_sat_QC' = as.double(oxy_sat_QC),
    'nitrate_QC' = as.double(nitrate_QC),
    'pH_QC' = as.double(pH_QC),
    'chl_QC' = as.double(chl),
    'poc_QC' = as.double(poc_QC))
  data_list[[j]] <- mat_test
  print(file[['dim']][['N_PROF']][['len']])
  nc_close(file)
}
library(dplyr)
soccom <- dplyr::bind_rows(data_list)
save(soccom, file = paste0('paper/data/soccom_all_variables_',Sys.getenv('soccom_date'),
                           '.RData'))
