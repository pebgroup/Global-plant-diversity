# Script to get climate variables using CRU TS for botanical countries. Averaging from 1979 - 2019. 
# Saves a "data/climate.rds" file: dataframe containing mean, sd and sample size (=n) for each variable.
# Sample size refers to the number of raster layer cells the values are calculated from.
# Note that sd = 0 for all variables with n=1

## RUNTIME: approx 10 minutes on 16 gb RAM


#### Data sources #######################################################
# CRU TS
# https://crudata.uea.ac.uk/cru/data/hrg/v4announcement.htm
# http://dap.ceda.ac.uk/thredds/fileServer/badc/cru/doc/File_Format_CRU_TS_v4.pdf
#  0.5° latitude by 0.5° longitude grid


# VARIABLES:
# PRE: annual precipitation (mm per year)
# PRS: precipitation seasonality (coefficient of variation of monthly values)
# MAT: mean annual temperature
# TRA: temperature annual range (maximum temperature of warmest month minus minimum temperature of coldest month)
# PET: potential evapotranspiration in mm/yr (sum each year, take average). 

# explanation for mean and sd calculations:
## MAT mean annual temperature is the mean of all temperatures in one place. 
## For each pixel, we calculate the mean of all annual temperature values. As every pixel has the same number 
## of temperature values, this is the overall mean. 

# "Note that PET values are mean mm/day for each month (with a scaling factor of 10 applied to the PET ascii (*.dat) files, but NOT the PET netcdf files (*.nc). The pet values in the datafiles therefore need to be multiplied by the number of days for each month to get the mean pet for that month."

# I put a buffer around small islands to actually capture all values. sometimes the shapefile and the raster layers are slightly off, this will adjust for it.



#### GET CLIMATe VALUES ################################################################

#### PET ####
#rasta <- cruts2raster("data/cru/cru_ts4.04.1901.2019.pet.dat.nc", timeRange=c("1979-01-01","2019-01-01"))
rasta <- stack("data/cru/cru_ts4.04.1901.2019.pet.dat.nc")
first_layer <- grep("1979", names(rasta))[1]
last_layer <- grep("2018", names(rasta))[length(grep("2018", names(rasta)))]
rasta <- subset(rasta, seq(first_layer, last_layer))



# get sum PET per year, then average over the years
year_start <- seq(1,480,12)
year_end <- seq(12,480,12)
annual <- c()
annual_stack <- c()
for(i in 1:length(year_start)){ # 
    annual <- calc(rasta[[year_start[i]:year_end[i]]], sum) # this step eliminates NAs
    if(i==1){annual_stack <- annual
    }else{
        annual_stack <- stack(annual_stack, annual)
    }
    if(!i%%1)cat(i,"\r")
}
rm(rasta)

annual_average_pet <- calc(annual_stack, mean)
rm(annual, annual_stack)
writeRaster(annual_average_pet, "processed_data/cru/annual_average_pet.grd", overwrite=TRUE)
rm(annual_average_pet)



#### TEMPERATURE ####
# get average per year, then average of all years
# rasta <- cruts2raster("data/cru/cru_ts4.04.1901.2019.tmp.dat.nc",
#                       timeRange=c("1979-01-01","2019-01-01"))
rasta <- stack("data/cru/cru_ts4.04.1901.2019.tmp.dat.nc")
rasta <- subset(rasta, seq(first_layer, last_layer))

# MAT
# The average of averages is equal to the average of all values if the number of elements of all groups is the same: just take average
annual_average_temp <- calc(rasta, mean)
writeRaster(annual_average_temp, "processed_data/cru/annual_average_temp.grd", overwrite=TRUE)
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
writeRaster(annual_temp_range, "processed_data/cru/annual_temp_range.grd", overwrite=TRUE)
rm(annual_temp_range)




#### PRECIPITATION ####
# get sum PRE per year, then average over the years
# get cv per year, then average over the years
# rasta <- cruts2raster("data/cru/cru_ts4.04.1901.2019.pre.dat.nc",
#                       timeRange=c("1979-01-01","2019-01-01"))
rasta <- stack("data/cru/cru_ts4.04.1901.2019.pre.dat.nc")
rasta <- subset(rasta, seq(first_layer, last_layer))
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
writeRaster(annual_average_pre, "processed_data/cru/annual_pre.grd", overwrite=TRUE)
rm(annual_average_pre)

seasonality_pre <- calc(annual_stack_cv, mean)
rm(annual_cv, annual_stack_cv)
writeRaster(seasonality_pre, "processed_data/cru/seasonality_pre.grd", overwrite=TRUE)
rm(seasonality_pre)






#### average for botanical countries ###############################################

## In the latter step we took mean values for all variables over the years 1979 - 2019 (serving as an averaged current climate)
## Following mean values are spatial means, they are averaged for each botanical country.
shape <- readOGR("shapefile/level3_mod.shp")

pet <- raster("processed_data/cru/annual_average_pet.grd")
mat <- raster("processed_data/cru/annual_average_temp.grd")
tra <- raster("processed_data/cru/annual_temp_range.grd")
pre <- raster("processed_data/cru/annual_pre.grd")
prs <- raster("processed_data/cru/seasonality_pre.grd")
#proj4string(pet) <- proj4string(shape)
#proj4string(shape) <- proj4string(pet)

par(mfrow=c(3,2), mar=c(3, 2, 2, 2))
plot(pet, main="annual PET")
plot(mat, main="annual temp")
plot(tra, main="annual temp range")
plot(pre, main="annual precipitation")
plot(prs, main="seasonality")

climate_vars <- c("pet", "mat", "tra", "pre", "prs")
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






# Paleoclimate ------------------------------------------------------------

# # 2021 (MAT only: 10.5281/zenodo.4568897)
# # No anomalies calculated yet.. so this is paleogeography
# rasta <- stack("data/paleoclim/MioMIP1.nc")
# names(rasta)[grep("Mid.Mio", names(rasta))]
# # first_layer <- grep("1979", names(rasta))[1]
# # last_layer <- grep("2018", names(rasta))[length(grep("2018", names(rasta)))]
# rasta.mio <- subset(rasta, names(rasta)[grep("Mid.Mio", names(rasta))])
# plot(rasta.mio)


# Bradshaw 2012 (https://doi.org/10.5194/cp-8-1257-2012)
#summary and anomalies on website:
#https://www.paleo.bristol.ac.uk/ummodel/data/tczth/standard_html/tczth.html
#download nc data here: https://www.paleo.bristol.ac.uk/cgi-bin/Page_new2.cgi #
#select experiment and present day standard to compare with
mio_pet_ano <- stack("data/paleoclim/tczth-tcztg_precip_ann_fsy1.nc")
#plot(pMAP_ano)

mio_mat_ano <- stack("data/paleoclim/tczth-tcztg_mat_ann_fsy1.nc")
#plot(pMAT_ano)

# stack them
pclimate_vars <- c("mio_pet_ano", "mio_mat_ano")
pstat_vars <- c("mean", "sd", "n")
combs <- nrow(expand.grid(pclimate_vars, pstat_vars))
m <- matrix(seq(1:combs), ncol=3, byrow = TRUE)
num_list <- split(m, rep(1:nrow(m)))

pres <- matrix(nrow=nrow(shape@data), ncol=combs)
rownames(pres) <- shape@data$LEVEL_3_CO
disag_id <- c()
upsale_count <- c()

# match projections
proj4string(mio_pet_ano) <- proj4string(shape)
proj4string(mio_mat_ano) <- proj4string(shape)

# match coordinates
bbox(mio_pet_ano)
bbox(shape)
apply(coordinates(mio_pet_ano), 2, range)
# 3.75 longitudes missing?
unique(coordinates(mio_pet_ano)[,1]) #  no, diff steps: 96 lng steps
unique(coordinates(mio_pet_ano)[,2])
#coordinates(mio_pet_ano)[,1] <- coordinates(mio_pet_ano)[,1]-180
extent(mio_pet_ano) # larger extent since pixel size 3.75 for lng + 1.5 for lat
# change coordinates
bb <- extent(-181.875, 181.875, -91.25, 91.25)
extent(mio_pet_ano) <- bb
extent(mio_mat_ano) <- bb

Sys.time()
for(i in 1:nrow(shape@data)){
  # loop over botanical countries
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]])
  for(j in 1:length(pclimate_vars)){
    # loop over each climate layer
    lay <- get(pclimate_vars[j])
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
    pres[i,num_list[[j]][1]] <- mean(rest)
    # set SD to zero for n=1 cases?
    #ifelse(is.na(sd(rest))==TRUE, stdiv <- 0, stdiv <- sd(rest))
    #res[i,num_list[[j]][2]] <- stdiv
    pres[i,num_list[[j]][2]] <- sd(rest)
    # sample size
    pres[i,num_list[[j]][3]] <- length(rest)
  }
  if(!i%%1)cat(i,"\r")
}
Sys.time()


