#########################################################
######## Rethinking Darwin's coral reef paradox #########
############################### R Morais Mar-2023 #######


## Options 
options(
	scipen=10^6, 
	dplyr.summarise.inform=FALSE)

## Packages
library(sf)
library(s2)
library(terra)
library(maps)
library(tidyverse)
library(units)
library(here)
library(rnaturalearth)


## Directories
reefdir <- here::here('Allen_coral_atlas')
outdir <- here::here('Outputs')
meowdir <- here::here('MEOW')
rastdir <- here::here('Rasters')
batdir <- here::here('GEBCO_bat')

## Functions
source('Code/functions.R')
detach('package:tm')
detach('package:NLP')


####################################################################################
#####  Loading coral reef distribution polygons (UNEP-WCMC) and chlorophyll   ######
#####   concentration from satellite (Nasa Ocean Color - AQUA MODIS 4km res)  ######
####################################################################################

#---------------------------------------------------------#
## Coral reef distribution polygons (Allen Coral Atlas)  ##
## It was rasterised and includes only the benthic cover ##
## dataset, and the categories coral/algae and rock      ##
#---------------------------------------------------------#

## Allen Coral Atlas coral/algae + rock distribution raster
reef_f <- rast(file.path(reefdir,"aca_coral_algae_rock_buffer_5-10km_res_1km.tiff"))
reef_c <- rast(file.path(reefdir,"aca_coral_algae_rock_buffer_0-5km_res_1km.tiff"))

## And a dataframe with coordinates and values
## filtering non-NAs (meaning, coral reef cells)
xyval_f <- cbind(
  xyFromCell(reef_f,cell=seq_len(ncell(reef_f))),
  val = values(reef_f)) |>
  as.data.frame() |>
  rename(val = `aca_coral_algae_rock_buffer_5-10km_res_1km`) |>
  filter(!is.na(val)) |>
  select(-val)

xyval_c <- cbind(
  xyFromCell(reef_c,cell=seq_len(ncell(reef_c))),
  val = values(reef_c)) |>
  as.data.frame() |>
  rename(val = `aca_coral_algae_rock_buffer_0-5km_res_1km`) |>
  filter(!is.na(val)) |>
  select(-val)


## And a SpatVec for operations
reef_pt_f <- vect(as.matrix(xyval_f))
reef_pt_c <- vect(as.matrix(xyval_c))
  

## Finally, saving a raster combining both buffers (0-5 and 
## 5-10 km) for non-remote sensing data extraction (filtering  
## in-situ chlorophyll, and nutrient data). First, resampling
## rasters to equalise dimensions, extent and value storage
reef_c <- resample(reef_c,reef_f)
reef_f <- resample(reef_f,reef_f)

## A plot shows it works
plot(reef_c, xlim = c(120,130),ylim = c(5,-5),col = 'red')
plot(reef_f, xlim = c(120,130),ylim = c(5,-5),add = TRUE)


reef <- reef_f | reef_c
reef <- app(reef, function(x) ifelse(!is.na(x),1,NA))
plot(reef, xlim = c(120,130),ylim = c(5,-5),col = 'blue')
## It works

writeRaster(reef,file.path(reefdir,"aca_coral_algae_rock_buffer_0-10km_res_1km.tiff"), overwrite=TRUE)

#-------------------------------------------------------------#
## Global chlorophyl distribution NASA Ocean Colour Level 3  ##
## 	  Aqua-MODIS 4 km resolution composite (2002-2023).      ##
## 		https://oceancolor.gsfc.nasa.gov/l3/order/ 			 ##
#-------------------------------------------------------------#

## Loads the layer using terra and
## interpolates it to a higher resolution
chl_nc4 <- rast(file.path(rastdir, 'AQUA_MODIS.20020704_20230731.L3m.CU.CHL.chlor_a.4km.nc')) 

## Sets the new resolution to 1km
x <- rast(chl_nc4)
res(x) <- 0.01

## And resampling (i.e., increasing resolution)
chl_nc1 <- resample(chl_nc4, x, method = 'cubic', threads = TRUE)


#----------------------------------------------------------------#
##  Global bathymetry dataset from GEBCO at the 15-arcsec res   ##
## 	         (approx. 0.45 km of resolution).                   ##
## 	gebco.net/data_and_products/gridded_bathymetry_data/#global ##
#----------------------------------------------------------------#

bat_nc <- rast(file.path(batdir, 'GEBCO_2024_sub_ice_topo.nc'))

## Sets the new resolution to 1km
x <- rast(bat_nc)
res(x) <- 0.01

## And resampling
bat_nc1 <- resample(bat_nc, x, method = 'cubic', threads = TRUE)

## Masking shallow pixels (<30 m depth)
bat_nc1[values(bat_nc1) > -30] <- NA
crs(bat_nc1) <- crs(chl_nc1)

## And masking chlorophyll data layer
chl_nc1_mask <- mask(chl_nc1,bat_nc1)


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
chl_reef_mf <- ext_filt(chl_nc1_mask, reef_pt_f, var = chlor_a) 
chl_reef_mc <- ext_filt(chl_nc1_mask, reef_pt_c, var = chlor_a)

## And extracting values for non-masked chlorophyll raster
chl_reef_nf <- ext_filt(chl_nc1, reef_pt_f, var = chlor_a)
chl_reef_nc <- ext_filt(chl_nc1, reef_pt_c, var = chlor_a)



## Converting into a vector that excludes newly made
## NA values, probably due to mismatch in the resolution
## of land between the NASA Ocean Colouur (original 4km res)
## and the ACA (original 5m res - 800x higher res



## Building an indicator table with the different strategies above
sens_list <- list(chl_reef_mc,chl_reef_nc,chl_reef_mf,chl_reef_nf)

sens_df <- data.frame(
  buffer = c("0-5km","0-5km","5-10km","5-10km"),
  bathym_mask = c('<30m masked','not-masked','<30m masked','not-masked'),
  median_chl_mgm3 = sapply(sens_list, function(x) median(x$chlor_a)),
  mean_chl_mgm3 = sapply(sens_list, function(x) mean(x$chlor_a)),
  prop_oligo =  sapply(sens_list, function(x) sum(x$chlor_a > 0.1) / nrow(x))
  ) |>
  mutate(across(median_chl_mgm3:prop_oligo, \(x) format(x, digits=3)))

write.csv(sens_df, file.path(outdir,'aca_sens_table.csv'), row.names=FALSE)    


#------------------------------------------------------#
##  Merges chlo conc data with MEOW biogeography info ##
##  including realm, province and ecoregion. It was   ##
##  downloaded directly from the UNEP-WCMC website:   ##
##  https://data.unep-wcmc.org/datasets/38            ##
#------------------------------------------------------#

## Continuing only with the masked and far from reef data
chl_reef <- chl_reef_mf

## Getting Spalding et al/s 2007 MEOW biogeographic data
all_meow <- read_sf(file.path(meowdir, 'meow_ecos.shp')) |>
  select(realm = REALM, lat_zone = Lat_Zone, province = PROVINCE, ecoregion = ECOREGION)

## A SpatVector object
meow <- vect(all_meow)
crs(meow) <- crs(chl_nc1_mask)

## And a vect form of the chl_reef_mf data
chl_reef_v <- vect(chl_reef, geom = c('x','y'), crs = crs(chl_nc1_mask))

## Extracting info from meow
chl_reef_meow <- terra::extract(meow,chl_reef_v)


## It works well but there are NAs
table(is.na(chl_reef_meow$realm))
wNA <- which(is.na(chl_reef_meow$realm))

coordNA <- chl_reef[wNA,c('x','y')]|> mutate(
  lat = trunc(y),
  lon = trunc(x)) |>
group_by(lon, lat) |>
  mutate(
    lat = unique(lat),
    lon = unique(lon))


coordNA_sf <- coordNA|>
  group_by(lon, lat) |>
  reframe(
    lat = unique(lat),
    lon = unique(lon),
    x = median(x),
    y = median(y)) |> 
  st_as_sf(coords = c('lon','lat'), crs = st_crs(all_meow))

## They all seem to be in 4 groups falling between the cracks

## Determining the closest ecoregion
dist <- st_distance(all_meow, coordNA_sf)
na_meow <- all_meow[apply(dist, 2, which.min),]
  
na_info <- bind_cols(
  coordNA_sf, 
  as.data.frame(na_meow) |> 
    select(-geometry)) |>
    mutate(
      lat = trunc(y),
      lon = trunc(x))

## Plotting: yep, it works
maps::map('world', fill=TRUE, border=NA, col='grey30')
plot(all_meow[,1], add=TRUE, col='transparent')
points(as.data.frame(coordNA)[,1], as.data.frame(coordNA)[,2],col='blue')
plot(na_meow[,1], add=TRUE, col=alpha('red',0.2))

## Merging it back 
coordNA_proc <- left_join(coordNA, na_info |> select(-x,-y), by = c('lat','lon')) |> 
  ungroup() |>
  select(realm, lat_zone, province,ecoregion)

chl_reef_meow[wNA,-1] <- coordNA_proc 

## and adjusting some contentious zones
reef_chlo <- chl_reef_meow |>
  select(realm, lat_zone, province, ecoregion) |>
  mutate(
    realm = case_when(
      realm == 'Temperate Australasia' | realm == 'Temperate Northern Pacific' ~ 'Central Indo-Pacific',
      realm == 'Temperate Southern Africa' ~ 'Western Indo-Pacific',
      TRUE ~ realm),
    x = chl_reef$x,
    y = chl_reef$y,
    chlor_a = chl_reef$chlor_a)

#------------------------------------------------------#
##  And finally determining distance from the coast   ##
##  so that reefs in inshore and offshore locations   ##
##  can be determined based on data distribution      ##
#------------------------------------------------------#

latlong <- st_crs("+proj=longlat")

## Complete data ##
reef_chlo_sf <- reef_chlo |>
  st_as_sf(
    coords = c('x', 'y'),
    crs=latlong)  


## World data normal latlong proj
w <- ne_countries(
  returnclass = 'sf') |>
  st_transform(crs=latlong)

## 3 polygon groups are not valid
## so making them valid
all(s2_is_valid(w))
w[which(!s2_is_valid(w)),] <- st_make_valid(w[which(!s2_is_valid(w)),])
all(s2_is_valid(w))

## If needed ##
sf_use_s2(FALSE)
st_is_valid(w)
w <- st_make_valid(w)

## Attention: this can take many hours, but not because of apply,
## of st_distance. Paralellising with mclapply can't help
dist_aff = apply(st_distance(reef_chlo_sf, w),MARGIN=1,FUN=min)/1000 # convert to km

## A quick look. It seems there may be a break between 120 and 240 km
ggplot(data=data.frame(d=dist_aff)) +
  geom_histogram(aes(x=d),bins=100) +
  scale_x_sqrt(breaks=c(0,30,60,120,240,480,960,1920,3840))

## Defining a split threshold
median(dist_aff[dist_aff > 0])
exp(mean(log(dist_aff[dist_aff > 0])))
thre <- exp(mean(log(c(60,360))))


## So defining inshore and offshore systems based on that
chl_reef_full <- reef_chlo |>
  mutate(habitat = ifelse(dist_aff < thre, 'inshore','offshore'),
         habitat = factor(habitat, c('inshore','offshore')))


## Saving this data
saveRDS(chl_reef_full, file.path(outdir,'reef_chlo_full.RDS'))


# After running message
cat("#---> ---> ---> Second script (Extracting chl-a data for near-reef cells) done.", sep='\n')