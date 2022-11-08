library(dplyr)
date <- Sys.getenv('soccom_date')
load(paste0('paper/data/soccom_all_variables_', date, '.RData'))
variable <- Sys.getenv('variable')
soccom[['variable_use_QC']] <- soccom[[paste0(variable, '_QC')]]
soccom[['variable_use']] <- soccom[[variable]]
# removals to BGC data
df_TS <- soccom %>%
  filter(psal > 28, latitude < -25, pressure < 2000, pressure > 0,
         variable_use_QC != 8) %>%
  filter(psal_QC == 0 |  temp_QC == 0 |  pressure_QC == 0 | 
           !is.na(pressure) | !is.na(longitude)) %>%
  mutate(profile_unique = paste0(float, '_', profile)) %>%
  group_by(profile_unique) %>%
  arrange(pressure) %>%
  mutate(bad2 = sum(psal < 33, na.rm = T), bad1= sum(psal > 37.2, na.rm = T),
         bad3 = max(abs(lead(pressure) - pressure), na.rm = T) > 200,
         bad4 = max(pressure, na.rm = T) < 1450,
         bad5 = min(pressure, na.rm = T) > 100,
         bad6 = n() < 15, 
         bad7 = sum( diff(pdens[pressure > 400])/diff(pressure[pressure > 400]) < -.05,
                     na.rm = T),
         bad8 = (sum(!is.na(variable_use_QC), na.rm = T) < 5) |
           (sum(variable_use_QC %in% c(0, 4), na.rm = T) < 5),
         bad9 = max(oxy > 400, na.rm = T)) %>%
  ungroup() %>%
  filter((bad1 + bad2 + bad3 + bad4 + bad5 + bad6 + bad7 + bad8 + bad9) ==0) %>%
  dplyr::select(-starts_with('bad'))


df_TS <- df_TS %>%
  mutate(day = as.Date(day, format = '%m/%d/%Y'),
         dayofyear = julian(day, origin = as.Date('2000-01-01')) %% 365.25)

df_TS_unique <- df_TS %>%
  filter(!duplicated(profile_unique)) %>%
  filter(!duplicated(cbind(longitude, latitude))) %>%
  dplyr::select(profile_unique)

df_TS <- df_TS %>%
  right_join(df_TS_unique)
# furthur removal to BGC data
df <- df_TS %>%
  filter_at(all_of(paste0(variable, '_QC')) , any_vars(. %in% c(0, 4))) %>%
  filter_at(all_of(variable) , any_vars(!is.na(.)))

if (variable != 'oxy') {
  df <- dplyr::select(df, -oxy, -oxy_sat)
}
if (variable != 'nitrate') {
  df <- dplyr::select(df, -nitrate)
}
if (variable != 'pH') {
  df <- dplyr::select(df, -pH)
}
if (variable != 'chl') {
  df <- dplyr::select(df, -chl)
}
if (variable != 'poc') {
  df <- dplyr::select(df, -poc)
}

df <- dplyr::select(df, -time, -variable_use)

df_list <- list(dplyr::select(df, -ends_with('_QC'),  # BGC info
                              - pdens)%>%
                  mutate(longitude = ifelse(longitude > 180, longitude - 360, longitude)), 
                dplyr::select(df_TS, -ends_with('_QC'), -pdens)%>% # info with all T/S measurements
                  mutate(longitude = ifelse(longitude > 180, longitude - 360, longitude))) 
save(df_list, file = paste0('paper/data/soccom_processed_', variable, '_', date, '.RData'))
