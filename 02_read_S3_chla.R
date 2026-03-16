source('remote-scripts/00_plotting.R')
library(tidyverse)
library(ncdf4)

nc <- nc_open("data/seychelles_chl_300m/cmems_obs-oc_glo_bgc-plankton_my_l4-olci-300m_P1M_CHL_50.00E-57.00E_9.00S-2.00S_2016-04-01-2025-10-01.nc")

# inspect structure
print(nc)

# extract
lon  <- ncvar_get(nc, "longitude")
lat  <- ncvar_get(nc, "latitude")
time <- ncvar_get(nc, "time") # Days Since 1900-01-01
chl  <- ncvar_get(nc, "CHL") # mass_concentration_of_chlorophyll_a_in_sea_water

nc_close(nc)

origin <- as.Date('1900-01-01')
dates <- origin + time

chla_mean <- apply(chl, c(1, 2), mean, na.rm = TRUE)

df <- expand.grid(lon = lon, lat = lat) %>%
  as_tibble() %>%
  mutate(chla_mean = as.vector(chla_mean))

# add sample size (ie cloud free days)
chla_n <- apply(chl, c(1, 2), function(x) sum(!is.na(x)))

df <- df %>%
  mutate(chla_n = as.vector(chla_n))

hist(df$chla_n)

# Plot
amira_raw<-ggplot() +
  geom_raster(data = df %>% filter(!is.na(chla_mean)),
              aes(x = lon, y = lat, fill = chla_mean)) +
  scale_fill_gradientn(
    colors = grad_cols,
    transform = 'log10',
    name = 'chl-a',
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  coord_sf(xlim = c(52.5, 54), ylim = c(-6.5, -4.7), expand = FALSE) +
  theme_minimal() +
  labs(x = NULL, y = NULL, subtitle = 'Sentinel-3 (~600m reso, raw)')


## load Sentinel as raster
chl_rast <- rast("data/seychelles_chl_300m/cmems_obs-oc_glo_bgc-plankton_my_l4-olci-300m_P1M_CHL_50.00E-57.00E_9.00S-2.00S_2016-04-01-2025-10-01.nc")
reso<-res(chl_rast)[1]

# Load bathy once outside the loop
bat_nc <- rast('data/gebco_2025_sub_ice_topo/GEBCO_2025_sub_ice.nc')
bat_nc <- terra::crop(bat_nc, ext(chl_rast))
bat_land <- bat_nc %>% mutate(elevation = ifelse(elevation > 0, 1, NA))

# Pre-allocate list
chl_masked_list <- vector("list", nlyr(chl_rast))

for (i in seq_len(nlyr(chl_rast))) {
  
  cat("Processing layer", i, "/", nlyr(chl_rast), "-", format(dates[i]), "\n")
  
  layer <- chl_rast[[i]]  # single SpatRaster layer
  
  chl_masked_list[[i]] <- masker(
    extent  = ext(chl_rast),
    reso    = reso, # 0.005555554 from chl S3
    chl     = layer,
    bathy  = bat_nc,
    cropper = FALSE,
    interpolate = FALSE
  )
  
}

# Stack results
chl_masked_stack <- rast(chl_masked_list)

chla_mean <- app(chl_masked_stack, mean, na.rm = TRUE)

df_mean <- as.data.frame(chla_mean, xy = TRUE) %>%
  as_tibble() %>%
  rename(lon = x, lat = y, chla_mean = mean)


# Plot
amira<-ggplot() +
  geom_raster(data = df_mean %>% filter(!is.na(chla_mean)),
              aes(x = lon, y = lat, fill = chla_mean)) +
  scale_fill_gradientn(
    colors = grad_cols,
    transform = 'log10',
    name = 'chl-a',
    limits=c(10^-1, 10^1),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  coord_sf(xlim = c(52.5, 54), ylim = c(-6.5, -4.7), expand = FALSE) +
  theme_minimal() +
  labs(x = NULL, y = NULL, subtitle = 'Sentinel-3 (~600m reso, masked to 30m)')


ggsave("maps/Sentinel3_Amirantes_mask.tiff",
       plot   = amira,
       width  = 10,
       height = 8,
       dpi    = 600)

ggsave("maps/Sentinel3_Amirantes_obs.tiff",
       plot   = amira_raw,
       width  = 10,
       height = 8,
       dpi    = 600)
