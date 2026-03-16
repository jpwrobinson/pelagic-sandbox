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

  # extract variable: dimensions = island x time
df <- chl %>% as.data.frame()
chla_mean <- apply(chl, c(1, 2), mean, na.rm = TRUE)

df <- expand.grid(lon = lon, lat = lat) %>%
  as_tibble() %>%
  mutate(chla_mean = as.vector(chla_mean))

df %>%
  ggplot(aes(x = lon, y = lat, fill = chla_mean)) +
  geom_raster() +
  coord_quickmap() +
  theme_minimal()


# add sample size (ie cloud free days)
chla_n <- apply(chl, c(1, 2), function(x) sum(!is.na(x)))

df <- df %>%
  mutate(chla_n = as.vector(chla_n))

hist(df$chla_n)


land <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")

# Log-transform for display (essential for chl-a)
df <- df %>%
  mutate(chla_log = log10(chla_mean))

# Plot
ggplot() +
  geom_raster(data = df %>% filter(!is.na(chla_mean)),
              aes(x = lon, y = lat, fill = chla)) +
  geom_sf(data = land, fill = "grey40", colour = NA) +
  scale_fill_viridis_c(
    option = "turbo",           # or "viridis", "plasma"
    name = expression(log[10]~Chl-a~(mg~m^{-3})),
    na.value = NA
  ) +
  coord_sf(xlim = c(50, 57), ylim = c(-9, -2), expand = FALSE) +
  theme_minimal() +
  labs(x = NULL, y = NULL)