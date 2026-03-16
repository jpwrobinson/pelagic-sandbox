## Extracting and masking chl-a datasets
library(sf)
library(s2)
library(terra)
library(tidyverse)
library(tidyterra)
source('remote-scripts/00_plotting.R')


#----------------------------------------------------------------#
##  Global bathymetry dataset from GEBCO at the 15-arcsec res   ##
## 	         (approx. 0.45 km of resolution).                   ##
## 	gebco.net/data_and_products/gridded_bathymetry_data/#global ##
#----------------------------------------------------------------#

# bat_nc <- rast('data/GEBCO_2024_sub_ice_topo.nc')

masker<-function(bathy, extent, reso, chl, cropper=FALSE, interpolate = TRUE){
  
  # Arguments
  # bathy = bathymetry file (use bat_nc above)
  # extent = to crop to focal islands / region
  # reso = chl-a resolution for matching with bathymetry (in degrees, e.g. for MODIS, 0.04166667 = 1/24 = 2.5' = ~4.6 km)
  # chl = netCDF of chl-a concentrations
  
  # Crop to focal extent
  # set to TRUE if using in single operation (ie not looping through nc layers)
  if(cropper == TRUE)
    {
  chl <- terra::crop(chl, extent)
  bat_nc <- terra::crop(bat_nc, extent)
  }

  # 1. Resample bathymetric layer for chl-a resolution
  x <- rast(bat_nc)
  res(x) <- reso

  ## And resample
  bat_nc1 <- resample(bat_nc, x, method = 'cubic', threads = TRUE)
  crs(bat_nc1) <- crs(chl)

  ## 2. Mask the shallow pixels (<30 m depth) by setting to NA
  bat_sh_mask<-bat_nc1
  bat_sh_mask[values(bat_sh_mask) > -30] <- NA

  ## 3. Now mask the chlorophyll data layer
  chl_mask <- mask(chl, bat_sh_mask)

  # 4. Interpolate shallow cells (but allow fails for NA in both layers)
  if(interpolate == TRUE){
  chl_approx <- tryCatch(
    approximate(c(chl_mask, chl)),
    error = function(e) chl_mask  # fall back to masked layer as-is
  )
  chl_mask<-chl_approx[[1]]
  }

  # 5. re-mask land values
  bat_la_mask<-bat_nc1
  bat_la_mask[values(bat_la_mask) >= 0] <- NA
  chl_mask <- mask(chl_mask,bat_la_mask)

  return(chl_mask)
  
  # Output = chl-a with shallow (<30m) pixels removed
}
