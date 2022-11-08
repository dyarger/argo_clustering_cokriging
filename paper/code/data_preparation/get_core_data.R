file_folder <- Sys.getenv('argo_location') 

library(tidyverse)
library(ncdf4)
traj_info <- read.delim(
  file = paste0(file_folder, "dac/ar_index_global_traj.txt"),
  skip = 8, header = T,
  sep = ",", stringsAsFactors = FALSE
)

# subset to floats with profiles < -30 latitude
traj_SO <- traj_info %>%
  filter(latitude_min < -30) %>%
  mutate(jurisdiction = sapply(strsplit(file, '/'), function(x) x[1]),
         float =  sapply(strsplit(file, '/'), function(x) x[2]))
prof_info <- read.delim(
  file = paste0(file_folder, "dac/ar_index_global_prof.txt"),
  skip = 8, header = T,
  sep = ",", stringsAsFactors = FALSE
)

prof_info <- prof_info %>%
  mutate(jurisdiction = sapply(strsplit(file, '/'), function(x) x[1]),
         float =  sapply(strsplit(file, '/'), function(x) x[2]))

prof_SO <- right_join(prof_info, traj_SO, by = c('jurisdiction', 'float'))
# subset to profiles after Jan 1, 2014
prof_SO_time <- prof_SO %>%
  mutate(Date = as.Date(substr(as.character(date), 1, 8), format ='%Y%m%d')) %>%
  filter(Date > as.Date('2014-01-01', format = '%Y-%m-%d'))

# Load files in to R
prof_files <- prof_SO_time[['file.x']]
prof_info <- sapply(strsplit(prof_files, '/'), function(x) x[4])
mode <- substr(prof_info, 1, 1)
float <- sapply(strsplit(substr(prof_info, 2, nchar(prof_info)), '_'), function(x) x[1])
profile <- 
  sapply(strsplit(
    sapply(strsplit(substr(prof_info, 2, nchar(prof_info)), '_'), function(x) x[2]),
    '\\.'
  ), function(x) x[1])

# consider only delayed mode profiles
prof_files <- prof_files[mode == 'D']
prof_files <- prof_files[!duplicated(prof_files)]
system('mkdir paper/data/core_profiles')
# Load information from each file
library(parallel)
return_vec <- mclapply(1:length(prof_files), mc.cores = 4, mc.preschedule = T,
                       prof_files = prof_files, file_folder = file_folder,
                       function(i, prof_files, file_folder) {
  file <- nc_open(paste0(file_folder, 'dac/', prof_files[i]))
  lat <- ncvar_get(file, "LATITUDE")[1]
  long <- ncvar_get(file, "LONGITUDE")[1]
  pos_qc <- strsplit(ncvar_get(file, "POSITION_QC"), "")[[1]][1]
  day <- ncvar_get(file, "JULD")[1]
  day_ref <- ncvar_get(file, "REFERENCE_DATE_TIME")[1]
  day_qc <- strsplit(ncvar_get(file, "JULD_QC"), "")[[1]][1] # quality control on day
  float <- as.numeric(ncvar_get(file, "PLATFORM_NUMBER"))[1] # float number
  cycle <- ncvar_get(file, "CYCLE_NUMBER")[1] # cycle number increments by 1 for each profile float collects
  mode <- substr(ncvar_get(file, "DATA_MODE")[1], 1, 1)
  if (mode %in% c("D", "A")) {
    if (length(dim(ncvar_get(file, "PRES_ADJUSTED"))) == 2) {
      pressure <- ncvar_get(file, "PRES_ADJUSTED")[, 1]
      temperature <- ncvar_get(file, "TEMP_ADJUSTED")[, 1]
      salinity <- ncvar_get(file, "PSAL_ADJUSTED")[, 1]
      pressure_qc <- strsplit(ncvar_get(file, "PRES_ADJUSTED_QC")[1], "")[[1]]
      temperature_qc <- strsplit(ncvar_get(file, "TEMP_ADJUSTED_QC")[1], "")[[1]]
      salinity_qc <- strsplit(ncvar_get(file, "PSAL_ADJUSTED_QC")[1], "")[[1]]
      pres_adj_error <- ifelse(is.na(ncvar_get(file, "PRES_ADJUSTED_ERROR")[, 1]),
                               Inf, ncvar_get(file, "PRES_ADJUSTED_ERROR")[, 1]
      )
    } else {
      pressure <- ncvar_get(file, "PRES_ADJUSTED")
      temperature <- ncvar_get(file, "TEMP_ADJUSTED")
      salinity <- ncvar_get(file, "PSAL_ADJUSTED")
      pressure_qc <- strsplit(ncvar_get(file, "PRES_ADJUSTED_QC"), "")[[1]]
      temperature_qc <- strsplit(ncvar_get(file, "TEMP_ADJUSTED_QC"), "")[[1]]
      salinity_qc <- strsplit(ncvar_get(file, "PSAL_ADJUSTED_QC"), "")[[1]]
      pres_adj_error <- ifelse(is.na(ncvar_get(file, "PRES_ADJUSTED_ERROR")[1]),
                               Inf, ncvar_get(file, "PRES_ADJUSTED_ERROR")[1]
      )
    }
  } else {
    nc_close(file)
    return('R')
  }
  nc_close(file)
  
  if (cycle == 0 | sum(pres_adj_error >= 20) > 0) {
    return('rejected')
  }
  year <- substr(day_ref, 1, 4)
  month <- substr(day_ref, 5, 6)
  day_ref_final <- substr(day_ref, 7, 8)
  hours <- substr(day_ref, 9, 10)
  minutes <- substr(day_ref, 11, 12)
  seconds <- substr(day_ref, 13, 14)
  date_time <- as.POSIXct(day * 60 * 60 * 24,
                          origin = as.POSIXct(ISOdatetime(
                            year = year, month = month, day = day_ref_final,
                            hour = hours, min = minutes,
                            sec = seconds, tz = "GMT"
                          )),
                          tz = "GMT"
  )
  date_use <- as.Date(date_time)
  
  df_return <- data.frame(pressure, temperature, salinity,
                          float, cycle, lat, long, pos_qc, day, day_qc,
                          pressure_qc, temperature_qc, salinity_qc, pres_adj_error,
                          reference = day_ref, date = date_use, date_time
  )
  single_prof <- df_return[df_return[['day_qc']] != 4 & 
                             df_return[['pressure_qc']] %in% c(1,2) &
                             df_return[['temperature_qc']] %in% c(1,2) & 
                             df_return[['salinity_qc']] %in% c(1,2), ]
  profile_unique <- paste0(single_prof[['float']][1], '_',
                           single_prof[['cycle']][1])
  save(single_prof, file = paste0('paper/data/core_profiles/profile_', profile_unique, '.RData'))
  
  if (i %% 1000 == 0) {
    print(i)
  }
  return('saved')
})
