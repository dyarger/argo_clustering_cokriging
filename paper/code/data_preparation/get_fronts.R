.libPaths()
library(tidyverse)
library(raster)
library(rgdal)
library(fields)
library(ncdf4)
library(gsw)

year <-  2014
month <- '09'
month_name <- month.abb[as.numeric(month)]
file_start <- "https://masie_web.apps.nsidc.org/pub/DATASETS/NOAA/G02135/south/monthly/geotiff/"
imported_raster <- raster(paste0(
  file_start, "/", month, "_", month_name, "/",
  "S_", year, month, "_concentration_v3.0.tif"
))
projection(imported_raster) <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs "
new_xy <- raster::projectRaster(
  from = imported_raster, res = c(.5/5, .25/5),
  crs = "+proj=longlat +datum=WGS84 +no_defs "#, method = "bilinear"
)
new_xy <- as.data.frame(new_xy, xy = T)
colnames(new_xy) <- c("long", "lat", "concentration")
RG_fine_grid <- nc_open(Sys.getenv('RG_location'))
temp_mean <- ncvar_get(RG_fine_grid, 'ARGO_TEMPERATURE_MEAN')
psal_mean <- ncvar_get(RG_fine_grid, 'ARGO_SALINITY_MEAN')
long <- RG_fine_grid[['dim']][['LONGITUDE']][['vals']] 
lat <- RG_fine_grid[['dim']][['LATITUDE']][['vals']]
pressure <- RG_fine_grid[['dim']][['PRESSURE']][['vals']]

SA <- array(gsw::gsw_SA_from_SP(SP = psal_mean,p = pressure, longitude = long, latitude = lat),
            dim = dim(psal_mean))
CT <- array(gsw::gsw_pt_from_t(SA = SA, t = temp_mean, p = pressure),
            dim = dim(psal_mean))

CT_100 <- CT[,,pressure == 100]
CT_400 <- CT[,,pressure == 400]

subtropical_front <- data.frame(expand.grid('longitude' = long, 'latitude' = lat), 
                                subtropical =  as.double(CT_100) < 11) %>%
  filter(latitude < -20)
subantarctic_front <- data.frame(expand.grid('longitude' = long, 'latitude' = lat), 
                                 subantarctic =  as.double(CT_400) < 5) %>%
  filter(latitude < -20)

new_xy_ord <- new_xy %>%
  filter(long > -180, long < 180) %>%
  mutate(long2 = ifelse(long <= 20, long + 360, long)) %>%
  arrange(long2, lat)

x_mat <- unique(new_xy_ord[['long2']])
y_mat <- unique(new_xy_ord[['lat']])
z_mat <- matrix(nrow = length(x_mat), ncol = length(y_mat), 
                new_xy_ord[['concentration']] > 150 &
                  new_xy_ord[['concentration']] <= 1000, byrow = T)
fields_obj <- list('x' = x_mat, 'y' = y_mat,
                   'z' = z_mat)

sea_ice_RG <- fields::interp.surface(obj = fields_obj, loc = expand.grid('longitude' = long, 'latitude' = c(seq(-75.5, -65.5, by = .25), lat)))

seaice_extent <- data.frame(expand.grid('longitude' = long, 'latitude' = c(seq(-75.5, -65.5, by = .25), lat)),
                            seaice = sea_ice_RG) %>%
  filter(latitude < -20) %>%
  mutate(seaice = ifelse(is.na(seaice), 0, seaice),
         seaice = seaice > .5)

for (long_i in unique(seaice_extent$longitude)) {
  max_lat <- max(seaice_extent$latitude[seaice_extent$longitude == long_i & 
                                      seaice_extent$seaice &
                                        seaice_extent$latitude < -56])
  seaice_extent$seaice[seaice_extent$longitude == long_i & 
                           seaice_extent$latitude < max_lat] <- T
}

all_fronts <- full_join(seaice_extent, subtropical_front) %>%
  left_join(subantarctic_front) %>%
  mutate(subantarctic = ifelse(is.na(subantarctic), T, subantarctic),
         subtropical = ifelse(is.na(subtropical), T, subtropical)) %>%
  pivot_longer(c(seaice, subtropical, subantarctic),names_to = 'front', values_to = 'South')
save(all_fronts, file = 'paper/data/fronts.RData')


source('paper/code/src/plot_source.R')
ggplot(data = all_fronts, aes(x = ifelse(longitude > 180, longitude - 360, longitude), 
                              y = latitude, color = front))+
  geom_contour(aes(z = as.numeric(South)),bins = 1) +
  SO_theme  +  SO_coord +  continents


