library(tidyverse)

load('paper/data/soccom_processed_oxy_05_05_21.RData')
theme_set(theme_bw())
soccom <- df_list[[2]][,c('longitude', 'latitude', 'day', 'pressure', 'temp', 'psal', 'oxy', 
                          'profile_unique')]
load('paper/data/core_processed_06_21.RData')
core_data <- core_data[, c('longitude', 'latitude', 'date', 'pressure', 'temp', 'psal', 'profile_unique')]

soccom_unique <- soccom[!duplicated(soccom[['profile_unique']]),]
core_unique <- core_data[!duplicated(core_data[['profile_unique']]),]

core_unique[['day']] <- core_unique[['date']]
all_unique <- rbind(cbind(soccom_unique[,c('longitude', 'latitude', 'day')], type = 'BGC'),
                    cbind(core_unique[,c('longitude', 'latitude', 'day')], type = 'Core'))
source('paper/code/src/plot_src.R')
theme_set(theme_bw())
continents <-   geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),
                             fill = 'white', color = 'gray25', size = .45)
ggplot(data = all_unique %>% filter(substr(day, 1, 4) == '2020'), aes(x = longitude, y = latitude))+
  geom_point(size = .2) + 
  facet_wrap(~type, ncol = 1) + 
  SO_theme + SO_coord + continents + 
  theme(text = element_text(size = 25))
ggsave('paper/images/data_2020.png', height = 8.86, width = 4.1)

floats <- sapply(strsplit(soccom$profile_unique, '_'), function(x) x[1])
unique_floats <- unique(floats)
cycle <- sapply(strsplit(soccom$profile_unique, '_'), function(x) x[2])

soccom_one <- soccom[floats == unique_floats[100] & as.numeric(cycle) %in% 1:10 ,] %>%
  rename(`Temperature (°C)` = temp, 
         `Salinity (Practical\nSalinity Units)` = psal,
         `Oxygen (μmol/kg)` = oxy)

soccom_one_long <- tidyr::pivot_longer(soccom_one, cols = c(`Temperature (°C)`, `Salinity (Practical\nSalinity Units)`,
                                                            `Oxygen (μmol/kg)`), 
                                       names_to = 'variable', values_to = 'value'
                                       ) %>%
  filter(!is.na(value))
soccom_one_long[['variable']] <- factor(soccom_one_long[['variable']],
                                   levels = c('Oxygen (μmol/kg)',
                                              'Temperature (°C)', 'Salinity (Practical\nSalinity Units)'))
ggplot(data =soccom_one_long, aes(x = pressure, y = value, group = profile_unique)) +
  geom_point(size = .1)+
  geom_line(alpha = .2)+
  facet_wrap(~variable, scales = 'free_x') + 
  coord_flip(xlim = c(2000, 0)) + 
  labs(x = 'Pressure (decibars)') + 
  theme(axis.title.x = element_blank())+ 
  theme(text = element_text(size = 18), axis.text.x = element_text(hjust = .75))
ggsave('paper/images/data_example.png', width = 9, height = 7.5)



ggplot(data =soccom_one_long %>%
         filter(variable != 'Oxygen (μmol/kg)'), aes(x = pressure, y = value, group = profile_unique)) +
  geom_point(size = .1)+
  geom_line(alpha = .2)+
  facet_wrap(~variable, scales = 'free_x') + 
  coord_flip(xlim = c(2000, 0)) + 
  labs(x = 'Pressure (decibars)') + 
  theme(axis.title.x = element_blank())+ 
  theme(text = element_text(size = 18), axis.text.x = element_text(hjust = .75))
ggsave('paper/images/data_example2.png', width = 7, height = 7.5)


soccom_one_unique <- soccom_one_long[!duplicated(soccom_one_long$profile_unique),]
library(ggrepel)
ggplot() +
  geom_point(data =soccom_one_unique %>% arrange(day), aes(x = longitude, y = latitude))+
  geom_path(data =soccom_one_unique %>% arrange(day), aes(x = longitude, y = latitude))+
  geom_label_repel(data = soccom_one_unique %>% arrange(day) %>% mutate(plot_it = c(T, rep(F, 8), T)),
                   aes(x = longitude, y = latitude, label = as.character(day),
                       alpha = plot_it),
                   max.overlaps = 1000, force_pull = 5, force = 5)+
  labs(x = 'Longitude', y = 'Latitude')+
  scale_alpha_discrete(range = c(0, 1)) + 
  theme(text = element_text(size = 18), axis.text.x = element_text(hjust = .75),
        legend.position = 'none')+
  continents + 
  #coord_quickmap(xlim = c(33, 120), ylim = c(-68, -30))
  coord_quickmap(xlim = c(63, 100), ylim = c(-57, -43))
ggsave('paper/images/data_example3.png', width = 7, height = 4)








