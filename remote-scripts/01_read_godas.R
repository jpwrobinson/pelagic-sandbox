library(sf)
library(terra)


# https://psl.noaa.gov/data/gridded/data.godas.html
# Behringer, D.W., M. Ji, and A. Leetmaa, 1998: An improved coupled model for ENSO prediction and implications for ocean initialization. Part I: The ocean data assimilation system. Mon. Wea. Rev., 126, 1013-1021.

# Long-term monthly means for SSH relative to geoid for 1991-2020
# 0.333 degree latitude x 1.0 degree longitude global grid (418x360)
# 74.5S - 64.5N, 0.5E - 359.5E

ssh<-rast('data/godas/sshg.mon.ltm.1991-2020.nc')

# plot(ssh) # this is a layer for each month
# names(ssh)
# crs(ssh)

# work with January to start
# jan<-ssh$sshg_1
# plot(jan)

# overlay reefs
source('00_islands.R')
latlon<-data.frame(lon = island$longitude, lat = island$latitude, island_group=island$island_group)
# align projection systems - longitude wrap
latlon$lon<-ifelse(latlon$lon < 0, 360+latlon$lon, latlon$lon)
# points(latlon)

ssh_vals_m<-terra::extract(ssh, latlon[,1:2]) %>% select(starts_with('sshg')) %>% 
  mutate(mean_annual_ssh =rowMeans(across(everything())), island = island$island, island_group=latlon$island_group) %>% 
  pivot_longer(-c(mean_annual_ssh, island, island_group), names_to = 'month', values_to = 'ssh') %>% 
  mutate(month_num = as.numeric(str_replace_all(month, 'sshg_', '')),
         month = month.abb[month_num])

ssh_vals_m_C<-ssh_vals_m %>% group_by(island_group, month_num, month) %>% 
  summarise(ssh = mean(ssh))

ssh_vals<-ssh_vals_m %>% group_by(island) %>% summarise(ssh = mean(ssh))
ssh_vals_C<-ssh_vals_m_C %>% group_by(island_group) %>% summarise(ssh = mean(ssh))
