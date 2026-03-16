## Extracting and masking chl-a 
## Packages
library(sf)
library(s2)
library(terra)
library(tidyverse)
library(tidyterra)
source('remote-scripts/00_plotting.R')
sey_extent<-c(52.5,57.5, -8, -3)
sey_extent2<-c(54,57, -5.5, -3)
ind_extent<-c(52,75, -8, 7.5)
mald_extent<-c(71.5, 75, 0, 7.5)

#-------------------------------------------------------------#
## Global chlorophyl distribution NASA Ocean Colour Level 3  ##
## 	  Aqua-MODIS 4 km resolution composite (2002-2023).      ##
## 		https://oceancolor.gsfc.nasa.gov/l3/order/ 			 ##
#-------------------------------------------------------------#

## Loads the layer using terra and
## interpolates it to a higher resolution
chl_nc4 <- terra::rast('data/AQUA_MODIS.20020704_20230731.L3m.CU.CHL.chlor_a.4km.nc')
crs(chl_nc4) <- 'EPSG:4326'
chl_nc4[values(chl_nc4) > 1] <- NA

#----------------------------------------------------------------#
##  Global bathymetry dataset from GEBCO at the 15-arcsec res   ##
## 	         (approx. 0.45 km of resolution).                   ##
## 	gebco.net/data_and_products/gridded_bathymetry_data/#global ##
#----------------------------------------------------------------#

bat_nc <- rast('data/GEBCO_2024_sub_ice_topo.nc')

## Sets the new resolution to 4km
x <- rast(bat_nc)
res(x) <- 0.04166667

## And resampling
bat_nc1 <- resample(bat_nc, x, method = 'cubic', threads = TRUE)
crs(bat_nc1) <- crs(chl_nc4)

## Masking shallow pixels (<30 m depth)
bat_mask<-bat_nc1
bat_mask[values(bat_mask) > -30] <- NA

## And masking chlorophyll data layer
chl_nc4_crop <- terra::crop(chl_nc4, ind_extent)
bat_mask_crop <- terra::crop(bat_mask, ind_extent)
chl_nc4_mask <- mask(chl_nc4_crop,bat_mask_crop)

# Now interpolate shallow cells
chl_approx<-approximate(c(chl_nc4_mask, chl_nc4_crop))
chl_approx<-chl_approx[[1]]

# re-mask land values
bat_mask2<-bat_nc1
bat_mask2[values(bat_mask2) >= 0] <- NA
bat_mask2_crop <- terra::crop(bat_mask2, ind_extent)
chl_approx_sey <- mask(chl_approx,bat_mask2_crop)

bat_crop <- terra::crop(bat_nc1, ind_extent)
bat_crop[values(bat_crop) >= 0] <- 0


#---------------------------------------------------------------#
## Extracts the chlorophyll-a values near  reefs (0-5, 5-10km) ##
## according to the Allen Coral Atlas benthic cat dataset      ##
## 	as per above (coral/algae + rock categories only)				   ##
#---------------------------------------------------------------#

## _mf = shallow-water masked 5-10 km from reefs
## _mc = shallow-water masked 0-5 km from reefs
## _nf = not masked 5-10 km from reefs
## _nc = not masked 0-5 km from reefs


## And for chlorophyll raster masked by bathymetry
## (all pixels shallower than 30m excluded)
# chl_reef_mf <- ext_filt(chl_nc1_mask, reef_pt_f, var = chlor_a) 

pdf(file = 'maps/chla_maps.pdf', height=7, width=8)

# Plot Seychelles.
ggplot() + 
  # geom_tile(data = chl_df, aes(x=x, y = y, fill=chlor_a)) +
  geom_spatraster(data = chl_approx, aes(fill = chlor_a)) +
  scale_x_continuous(expand=c(0,0), limits=sey_extent2[c(1,2)]) +
  scale_y_continuous(expand=c(0,0), limits=sey_extent2[c(3,4)]) +
  scale_fill_gradientn(
    colors = grad_cols,
    transform = 'log10',
    name = 'chl-a',
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  white_axis_theme

ggplot() + 
  # geom_tile(data = chl_df, aes(x=x, y = y, fill=chlor_a)) +
  geom_spatraster(data = chl_approx, aes(fill = chlor_a)) +
  scale_x_continuous(expand=c(0,0), limits=mald_extent[c(1,2)]) +
  scale_y_continuous(expand=c(0,0), limits=mald_extent[c(3,4)]) +
  scale_fill_gradientn(
    colors = grad_cols,
    transform = 'log10',
    name = 'chl-a',
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    white_axis_theme
dev.off()


