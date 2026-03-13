library(tidyverse)
library(ncdf4)

nc <- nc_open("data/seychelles_chl_300m/cmems_obs-oc_glo_bgc-plankton_my_l4-olci-300m_P1M_CHL_50.00E-57.00E_9.00S-2.00S_2016-04-01-2025-10-01.nc")

# inspect structure
print(nc)

# extract
lon  <- ncvar_get(nc, "longitude")
lat  <- ncvar_get(nc, "latitude")
time <- ncvar_get(nc, "time")
chl  <- ncvar_get(nc, "CHL")

nc_close(nc)

