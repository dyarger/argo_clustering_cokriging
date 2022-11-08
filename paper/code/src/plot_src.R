library(ggplot2)

continents <-   geom_polygon(data = map_data('world'), aes(x = long, y = lat, group = group),
                             fill = 'white', color = 'gray25', size = .2)
SO_coord <- coord_map('ortho', orientation = c(-90, 0, 0))
SO_theme <-   theme(axis.ticks = element_blank(), axis.text = element_blank(),
                    axis.title = element_blank(),
                    legend.position = 'bottom')
plot_lines <- rbind( data.frame(expand.grid('latitude' = seq(-70, -30, by = 10),
                                            'longitude' = seq(-180, 180, by = 2)),
                                type = 'latitude'),
                     data.frame(expand.grid('latitude' = seq(-70, -30, by = 2),
                                            'longitude' = seq(-180, 180, by = 30)),
                                type = 'longitude'))
latitude_lines <- geom_line(data = filter(plot_lines, type == 'latitude'), aes(x = longitude, y = latitude, group = latitude), 
                            color = 'gray35', size = .1)
longitude_lines <- geom_line(data = filter(plot_lines, type == 'longitude'), aes(x = longitude, y = latitude, group = longitude), 
                             color = 'gray35', size = .1)

labels_df <- data.frame(type = c('means', 'upper_pred', 'lower_pred', 'upper_total', 'lower_total',
                                 'sd_pred', 'sd_total', 'sd_meas'),
                        Label = c('Mean', '', '', '2SD bound', '2SD bound', 'Pred SD', 'Marg SD',
                                  'Meas Error SD'))

load('paper/data/fronts.RData')
fronts_dark <-   geom_contour( data = all_fronts,
                               aes(x = ifelse(longitude > 180, longitude - 360, longitude), 
                                   y = latitude, group = front, z = as.numeric(South)),bins = 1,
                               color = 'gray34', size = .24)

fronts_light <-   geom_contour( data = all_fronts,
                                aes(x = ifelse(longitude > 180, longitude - 360, longitude), 
                                    y = latitude, group = front, z = as.numeric(South)),bins = 1,
                                color = 'gray85', size = .24)

