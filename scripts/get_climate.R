# Script to get climate variables using CRU TS for botanical countries. Averaging from 1979 - 2019. 
# Saves a "data/climate.rds" file: dataframe containing mean, sd and sample size (=n) for each variable.
# Sample size refers to the number of raster layer cells the values are calculated from.
# Note that sd = 0 for all variables with n=1, this should also apply to PET values however across the oceans values are just set to 0 in the raster file. fix?

## RUNTIME: approx 10 minutes on 16 gb RAM


#### Data sources #######################################################
# CRU TS
# https://crudata.uea.ac.uk/cru/data/hrg/v4announcement.htm
# http://dap.ceda.ac.uk/thredds/fileServer/badc/cru/doc/File_Format_CRU_TS_v4.pdf
#  0.5° latitude by 0.5° longitude grid = approx 50km


# VARIABLES:
# PRE: annual precipitation (mm per year)
# SEA: precipitation seasonality (coefficient of variation of monthly values)
# MAT: mean annual temperature
# TRA: temperature annual range (maximum temperature of warmest month minus minimum temperature of coldest month)
# PET: potential evapotranspiration in mm/yr (sum each year, take average). 

# explanation for mean and sd calculations:
## MAT mean annual temperature is the mean of all temperatures in one place. 
## For each pixel, we calculate the mean of all annual temperature values. As every pixel has the same number 
## of temperature values, taking the overall mean works. 

# not implemented yet because simple multiplication:
# "Note that PET values are mean mm/day for each month (with a scaling factor of 10 applied to the PET ascii (*.dat) files, but NOT the PET netcdf files (*.nc). The pet values in the datafiles therefore need to be multiplied by the number of days for each month to get the mean pet for that month."

# I will put a buffer around small islands to actually capture all values. sometimes the shapefile and the raster layers are slightly off, this will adjust for it.



#### RUN ###############################################################

library(cruts)
library(raster)
library(vegan)


#### GET CLIMATOLOGICAL VALUES ################################################################

#### PET ####
rasta <- cruts2raster("data/cru/cru_ts4.04.1901.2019.pet.dat.nc",
                      timeRange=c("1979-01-01","2019-01-01"))

# plot(rasta[[c(1:12)]]) # plots first 12 months
# plot(calc(rasta[[c(1:12)]], mean, na.rm=TRUE)) # mean of the 12 months

# get sum PET per year, then average over the years
year_start <- seq(1,480,12)
year_end <- seq(12,480,12)
annual <- c()
annual_stack <- c()
for(i in 1:length(year_start)){ # 
    annual <- calc(rasta[[year_start[i]:year_end[i]]], sum) # this step somehow eliminates NAs
    if(i==1){annual_stack <- annual
    }else{
        annual_stack <- stack(annual_stack, annual)
    }
    if(!i%%1)cat(i,"\r")
}
rm(rasta)

annual_average_pet <- calc(annual_stack, mean)
rm(annual, annual_stack)
writeRaster(annual_average_pet, "data/cru/annual_average_pet.grd", overwrite=TRUE)
rm(annual_average_pet)



#### TEMPERATURE ####
# get average per year, then average of all years
rasta <- cruts2raster("data/cru/cru_ts4.04.1901.2019.tmp.dat.nc",
                      timeRange=c("1979-01-01","2019-01-01"))

# MAT
# The average of averages is only equal to the average of all values in two cases: if the number of elements of all groups is the same; or. the trivial case when all the group averages are zero. 
# --> case one applies to us, all years have 12 values
annual_average_temp <- calc(rasta, mean)
writeRaster(annual_average_temp, "data/cru/annual_average_temp.grd" ,overwrite=TRUE)
rm(annual_average_temp)

# temperature annual range (maximum temperature of warmest month minus minimum temperature of coldest month)
annual_stack <- c()
for(i in 1:length(year_start)){ # 
  annual_max <- calc(rasta[[year_start[i]:year_end[i]]], max)
  annual_min <- calc(rasta[[year_start[i]:year_end[i]]], min)
  range <- annual_max - annual_min
  if(i==1){annual_stack <- range
  }else{
    annual_stack <- stack(annual_stack, range)
  }
  if(!i%%1)cat(i,"\r")
}
rm(rasta)

annual_temp_range <- calc(annual_stack, mean)
rm(annual_stack, annual_max, annual_min, range)
writeRaster(annual_temp_range, "data/cru/annual_temp_range.grd", overwrite=TRUE)
rm(annual_temp_range)




#### PRECIPITATION ####
# get sum PRE per year, then average over the years
# get cv per year, then average over the years
rasta <- cruts2raster("data/cru/cru_ts4.04.1901.2019.pre.dat.nc",
                      timeRange=c("1979-01-01","2019-01-01"))
annual_sum <- c()
annual_stack_sum <- c()
annual_cv <- c()
annual_stack_cv <- c()

cv.fun <- function(x, ...){sd(x)/mean(x)}

for(i in 1:length(year_start)){ # coefficient of variation of monthly precipitation
  annual_sum <- calc(rasta[[year_start[i]:year_end[i]]], sum)
  annual_cv<- calc(rasta[[year_start[i]:year_end[i]]], cv.fun)
  if(i==1){annual_stack_sum <- annual_sum
          annual_stack_cv <- annual_cv
  }else{
    annual_stack_sum <- stack(annual_stack_sum, annual_sum)
    annual_stack_cv <- stack(annual_stack_cv, annual_cv)
  }
  if(!i%%1)cat(i,"\r")
}
rm(rasta)

annual_average_pre <- calc(annual_stack_sum, mean)
rm(annual_sum, annual_stack_sum)
writeRaster(annual_average_pre, "data/cru/annual_pre.grd", overwrite=TRUE)
rm(annual_average_pre)

seasonality_pre <- calc(annual_stack_cv, mean)
rm(annual_cv, annual_stack_cv)
writeRaster(seasonality_pre, "data/cru/seasonality_pre.grd", overwrite=TRUE)
rm(seasonality_pre)






#### EXTRACT FOR BOTANICAL COUNTRIES ###############################################

## In the latter step we took mean values for all variables over the years 1979 - 2019 (serving as an averaged current climate)
## Mean values now are spatial means, they are averaged over a botanical country. Same goes for sd
library(rgdal)
library(geosphere)
shape <- readOGR("shapefile/level3_mod.shp")



pet <- raster("data/cru/annual_average_pet.grd")
mat <- raster("data/cru/annual_average_temp.grd")
tra <- raster("data/cru/annual_temp_range.grd")
pre <- raster("data/cru/annual_pre.grd")
sea <- raster("data/cru/seasonality_pre.grd")
proj4string(shape) <- proj4string(pet)

par(mfrow=c(3,2), mar=c(3, 2, 2, 2))
plot(pet, main="annual PET")
plot(mat, main="annual temp")
plot(tra, main="annual temp range")
plot(pre, main="annual precipitation")
plot(sea, main="seasonality")

climate_vars <- c("pet", "mat", "tra", "pre", "sea")
stat_vars <- c("mean", "sd", "n")
combs <- nrow(expand.grid(climate_vars, stat_vars))
m <- matrix(seq(1:combs), ncol=3, byrow = TRUE)
num_list <- split(m, rep(1:nrow(m)))

res <- matrix(nrow=nrow(shape@data), ncol=combs)
rownames(res) <- shape@data$LEVEL_3_CO
disag_id <- c()
upsale_count <- c()

Sys.time()
for(i in 1:nrow(shape@data)){
  # loop over botanical countries
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]])
  for(j in 1:length(climate_vars)){
    # loop over each climate layer
    lay <- get(climate_vars[j])
    rest <- raster::extract(lay, shape_sub)
    rest <- na.omit(rest[[1]])
    # increase resolution necessary?
    if(all(is.na(rest))==TRUE){
      print("disaggregate to increase resolution")
      #disag <- c(disag, 1)
      upsale_count <- c(upsale_count, 1)
      disag_id <- c(disag_id, paste(i,j))
      # to avoid huge raster files crop to extent of shapefile sub * 20
      newExtent <- extent(bbox(shape_sub))
      lay2 <- crop(lay, newExtent*20)
      lay2 <- disaggregate(lay2, 10)
      rest <- raster::extract(lay2, shape_sub)
      rest <- na.omit(rest[[1]])
      
    }else{}
    
    # get mean and sd
    res[i,num_list[[j]][1]] <- mean(rest)
    # set SD to zero for n=1 cases?
      #ifelse(is.na(sd(rest))==TRUE, stdiv <- 0, stdiv <- sd(rest))
      #res[i,num_list[[j]][2]] <- stdiv
    res[i,num_list[[j]][2]] <- sd(rest)
    # sample size
    res[i,num_list[[j]][3]] <- length(rest)
  }
  if(!i%%1)cat(i,"\r")
}
Sys.time()

disag_id

temp <- c(paste(climate_vars[1], stat_vars, sep="_"),
          paste(climate_vars[2], stat_vars, sep="_"),
          paste(climate_vars[3], stat_vars, sep="_"),
          paste(climate_vars[4], stat_vars, sep="_"),
          paste(climate_vars[5], stat_vars, sep="_"))
res <- as.data.frame(res)
names(res) <- temp


saveRDS(res, "data/climate.rds")





#### missing data ########################################################### 
library(VIM)
aggr_plot <- aggr(res, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(data), cex.axis=.9, gap=3, ylab=c("Histogram of missing data","Pattern"))

# no data available is usually tiny islands
# missing SDs are due to just one data point available. hence no SD
# 10 botanical countries missing, 36 for PET. 0.83% of data is complete

table(is.na(res$pet_mean))
table(is.na(res$mat_mean))

# 
# 
# #### CAPTURING MORE DATA WITH BUFFER #####################################
# library(sf)
# # getting buffer for shapes missing data, here: PET
# pet <- raster("data/cru/annual_average_pet.grd")
# 
# miss.pet <- rownames(res)[is.na(res$pet_mean)]
# plot(shape[shape$LEVEL_3_CO%in% miss.pet,], lwd=1, border=shape$LEVEL_3_CO)
# shape2 <- st_as_sf(shape) # create sf object
# buff <- st_buffer(shape2[shape2$LEVEL_3_CO %in% miss.pet,], dist = 0.02)
# plot(st_geometry(buff), add=TRUE)
# 
# 
# # use buffer
# value <- c()
# stand_dev <- c()
# n <- c()
# for(i in 1:nrow(buff)){
#   rest <- extract(pet, buff[i,])
#   if(i>1){print(rest)}
#   rest <- na.omit(rest[[1]])
#   value <- c(value, mean(rest))
#   stand_dev <- c(stand_dev, sd(rest))
#   n <- c(n, length(rest))
#   if(!i%%1)cat(i,"\r")
# }
# 
# # does not do anything for PET
# 
# 
# # getting buffer for shapes missing data, here: MAT
# mat <- raster("data/cru/annual_average_temp.grd")
# 
# miss.mat <- rownames(res)[is.na(res$mat_mean)]
# plot(shape[shape$LEVEL_3_CO%in% miss.mat,], lwd=1, border=shape$LEVEL_3_CO)
# shape2 <- st_as_sf(shape) # create sf object
# buff <- st_buffer(shape2[shape2$LEVEL_3_CO %in% miss.mat,], dist = 0.2)
# plot(st_geometry(buff), add=TRUE)
# 
# 
# # use buffer
# value <- c()
# stand_dev <- c()
# n <- c()
# for(i in 1:nrow(buff)){
#   rest <- extract(mat, buff[i,])
#   rest <- na.omit(rest[[1]])
#   value <- c(value, mean(rest))
#   stand_dev <- c(stand_dev, sd(rest))
#   n <- c(n, length(rest))
#   if(!i%%1)cat(i,"\r")
# }

# captures 2 more, but rather pointless as they will be excluded no matter what as PET is missing


###########################################################################


# # dealing with small polygons
# newExtent <- extent(bbox(subset(shape, shape$LEVEL_3_CO=="CAY")))
# mat2 <- crop(mat, newExtent*20)
# mat_high <- disaggregate(mat2, 10)
# 
# plot(subset(shape, shape$LEVEL_3_CO=="CAY"))
# #plot(mat, add=TRUE)
# plot(rasterToPolygons(mat_high), add=TRUE, border="black")
# 
# extract(mat, subset(shape, shape$LEVEL_3_CO=="CAY"))
# extract(mat_high, subset(shape, shape$LEVEL_3_CO=="CAY"))
# 







# 
# 
# 
# # #### CHELSA DATA ###############################################################
# 
# # CHELSA is free and comes as tif - has no evapotranspiration though
# # download on 25th June from http://chelsa-climate.org/downloads/
# # For climatological values, a file contains the means of a given period (usually 1979-2013)
# # advantage: climatological values are already calculated
# # 30 arc sec = approx 1km
# 
# 
# library(rgdal)
# 
# shape <- readOGR("shapefile/level3.shp")
# 
# mat <- raster("data/chelsa/CHELSA_bio10_01_mean_annual_temperature.tif")
# tar <- raster("data/chelsa/CHELSA_bio10_07_temperature_annual_range.tif")
# prec <- raster("data/chelsa/CHELSA_bio10_12_annual_precipitation.tif")
# prec_season <- raster("data/chelsa/CHELSA_bio10_15_precipitation_seasonality.tif")
# 
# 
# cru <- readRDS("data/climate.rds")
# cru_missing <- rownames(cru)[which(cru$mat_n==0)]
# 
# 
# Sys.time()
# mat_mean <- c()
# mat_sd <- c()
# even_mat <- c()
# tar_mean <- c()
# tar_sd <- c()
# even_tar <- c()
# prec_mean <- c()
# prec_sd <- c()
# even_prec  <- c()
# prec_season_mean <- c()
# prec_season_sd <- c()
# even_prec_season <- c()
# 
# for(i in 2:length(cru_missing)){ #nrow(shape@data)
#   shape_sub <- subset(shape, shape$LEVEL_3_CO==cru_missing[i])
# #  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]])
# 
#   # MAT
#   rest <- extract(mat, shape_sub)
#   rest <- na.omit(rest[[1]])
#   mat_mean <- c(mat_mean, mean(rest))
#   mat_sd <- c(mat_sd, sd(rest))
#   even_mat <- c(even_mat, diversity(table(rest))/log(length(unique(rest))))
#   print("check1")
# 
#   # TAR
#   rest <- extract(tar, shape_sub)
#   rest <- na.omit(rest[[1]])
#   tar_mean <- c(tar_mean, mean(rest))
#   tar_sd <- c(tar_sd, sd(rest))
#   even_tar <- c(even_tar, diversity(table(rest))/log(length(unique(rest))))
#   print("check2")
#   
#   # PREC
#   rest <- extract(prec, shape_sub)
#   rest <- na.omit(rest[[1]])
#   prec_mean <- c(prec_mean, mean(rest))
#   prec_sd <- c(prec_sd, sd(rest))
#   even_prec <- c(even_prec, diversity(table(rest))/log(length(unique(rest))))
#   print("check3")
# 
#   # PREC_SEASON
#   rest <- extract(prec_season, shape_sub)
#   rest <- na.omit(rest[[1]])
#   prec_season_mean <- c(prec_season_mean, mean(rest))
#   prec_season_sd <- c(prec_season_sd, sd(rest))
#   even_prec_season <- c(even_prec_season, diversity(table(rest))/log(length(unique(rest))))
# 
#   if(!i%%1)cat(i,"\r")
# }
# Sys.time()
# 
# res <- data.frame(country=cru_missing[2:length(cru_missing)],
#                   mat_mean, mat_sd, even_mat,
#                   tar_mean, tar_sd, even_tar,
#                   prec_mean, prec_sd, even_prec,
#                   prec_season_mean, prec_season_sd, even_prec_season)
# # res <- data.frame(country=shape$LEVEL_3_CO[1],
# #                   mat_mean, mat_sd, even_mat,
# #                   tar_mean, tar_sd, even_tar,
# #                   prec_mean, prec_sd, even_prec,
# #                   prec_season_mean, prec_season_sd, even_prec_season)
# 
# saveRDS(res, file="chelsa.rds")
# 
# 
