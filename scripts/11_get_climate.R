# Script to get climate variables using CRU TS for botanical countries.
# Averaging from 1979 - 2019. Saves a "data/climate.rds" file: dataframe
# containing mean, sd and sample size (=n) for each variable. Sample size refers
# to the number of raster layer cells the values are calculated from. Note that
# sd = 0 for all variables with n=1

## RUNTIME: approx 10 minutes


#### Data sources #######################################################
# CRU TS
# https://crudata.uea.ac.uk/cru/data/hrg/v4announcement.htm
# http://dap.ceda.ac.uk/thredds/fileServer/badc/cru/doc/File_Format_CRU_TS_v4.pdf
#  0.5° latitude by 0.5° longitude grid


# VARIABLES: PRE: annual precipitation (mm per year) PRS: precipitation
# seasonality (coefficient of variation of monthly values) MAT: mean annual
# temperature TRA: temperature annual range (maximum temperature of warmest
# month minus minimum temperature of coldest month) PET: potential
# evapotranspiration in mm/yr (sum each year, take average).

#About mean and sd calculations: # MAT mean annual temperature is the mean of
#all temperatures in one place. # For each pixel, we calculate the mean of all
#annual temperature values. As every pixel has the same number # of temperature
#values, this is the overall mean.


# I put a buffer around small islands in the shapefile to deal with small
# offsets between shapefile and raster layers, this will adjust for it without
# including inaccurate values as there are not other values around islands.



library(raster)
library(rgdal)


# Present climate ---------------------------------------------------------

#### PET ####
rasta <- stack("../data/cru/cru_ts4.04.1901.2019.pet.dat.nc")
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
writeRaster(annual_average_pet, "../processed_data/cru/annual_average_pet.grd", overwrite=TRUE)
rm(annual_average_pet)



#### TEMPERATURE ####
# get average per year, then average of all years
rasta <- stack("../data/cru/cru_ts4.04.1901.2019.tmp.dat.nc")
rasta <- subset(rasta, seq(first_layer, last_layer))

# MAT
# The average of averages is equal to the average of all values if the number of elements of all groups is the same: just take average
annual_average_temp <- calc(rasta, mean)
writeRaster(annual_average_temp, "../processed_data/cru/annual_average_temp.grd", overwrite=TRUE)
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
writeRaster(annual_temp_range, "../processed_data/cru/annual_temp_range.grd", overwrite=TRUE)
rm(annual_temp_range)




#### PRECIPITATION ####
# get sum PRE per year, then average over the years
# get cv per year, then average over the years
rasta <- stack("../data/cru/cru_ts4.04.1901.2019.pre.dat.nc")
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
writeRaster(annual_average_pre, "../processed_data/cru/annual_pre.grd", overwrite=TRUE)
rm(annual_average_pre)

seasonality_pre <- calc(annual_stack_cv, mean)
rm(annual_cv, annual_stack_cv)
writeRaster(seasonality_pre, "../processed_data/cru/seasonality_pre.grd", overwrite=TRUE)
rm(seasonality_pre)






#### Average for botanical countries ##########################################

shape <- readOGR("../data/shapefile_bot_countries/level3_mod.shp")

pet <- raster("../processed_data/cru/annual_average_pet.grd")
mat <- raster("../processed_data/cru/annual_average_temp.grd")
tra <- raster("../processed_data/cru/annual_temp_range.grd")
pre <- raster("../processed_data/cru/annual_pre.grd")
prs <- raster("../processed_data/cru/seasonality_pre.grd")

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
      upsale_count <- c(upsale_count, 1)
      disag_id <- c(disag_id, paste(i,j))
      # to avoid huge raster files crop to extent of shapefile sub * 20
      newExtent <- extent(bbox(shape_sub))
      lay2 <- crop(lay, newExtent*20)
      lay2 <- disaggregate(lay2, 10)
      rest <- raster::extract(lay2, shape_sub)
      rest <- na.omit(rest[[1]])
      
    }else{}
    
    # get mean, sd and sample size
    res[i,num_list[[j]][1]] <- mean(rest)
    res[i,num_list[[j]][2]] <- sd(rest)
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

saveRDS(res, "../processed_data/climate.rds")





# Paleoclimate ------------------------------------------------------------

# Bradshaw 2012 (https://doi.org/10.5194/cp-8-1257-2012)
#summary and anomalies on website:
#https://www.paleo.bristol.ac.uk/ummodel/data/tczth/standard_html/tczth.html
#download nc data here: https://www.paleo.bristol.ac.uk/ummodel/users/Bradshaw_et_al_2012/new2/ #
# https://www.paleo.bristol.ac.uk/cgi-bin/Page_new2.cgi
#select experiment and present day standard to compare with


# Miocene precipitation anomaly ####
path <- "../data/paleoclim/precipitation/"
lay <- paste0(path, dir(path)[grep("^tcztk-tczti", dir(path))])
lay2 <- paste0(path, dir(path)[grep("^tczti", dir(path))])
lay3 <- paste0(path, dir(path)[grep("^tcztk_", dir(path))])
months <- letters[1:12]
months2 <- letters[13:24]
months3 <- paste0("n",letters[13:24])
for(z in 1:length(lay2)){
  temp <- raster(lay[z])
  temp2 <- raster(lay2[z])
  temp3 <- raster(lay3[z])
  assign(months[z], temp)
  assign(months2[z], temp2)
  assign(months3[z], temp3)
}
mio_pre_12m_ano <- stack(a,b,c,d,e,f,g,h,i,j,k,l)
mio_pre_12m <- stack(m,n,o,p,q,r,s,t,u,v,w,x)
pi_pre_12m <- stack(nm,nn,no,np,nq,nr,ns,nt,nu,nv,nw,nx)

mio_pre_ano <- calc(mio_pre_12m_ano, sum)
mio_pre <- calc(mio_pre_12m, sum)
pi_pre <- calc(pi_pre_12m, sum)
# note: the sum(one)-sum(two) == sum(one-two)

# mio_pre_ano should be the same as pi_pre - mio_pre
mio_pre_ano2 <- calc(stack(mio_pre, pi_pre), diff)
par(mfrow=c(2,1))
plot(mio_pre_ano)
plot(mio_pre_ano2)
# all is well


# TEMPERATURE ANOMALY #####
mio_mat_ano <- stack("../data/paleoclim/tcztk-tczti_mat_ann_fsy1.nc")
mio_mat <- stack("../data/paleoclim/tczti_mat_ann_fsy1.nc")


# Recenter maps #####
recenter.map <- function(x){
  x1 <- crop(x, extent(-1.875, 178.125, -91.25, 91.25))
  x2 <- crop(x, extent(178.125, 358.125, -91.25, 91.25))   
  extent(x2) <- c(-181.875, -1.875, -91.25, 91.25)
  m <- merge(x1, x2)
  return(m)
}
mio_list <- lapply(mget(ls()[grep("mio_", ls())]), recenter.map)
# unlist items
for(i in 1:length(mio_list)){assign(names(mio_list[i]), mio_list[[i]])}
pi_pre <- recenter.map(pi_pre)

# LGM ---------------------------------------------------------------------

mat_cur <- raster("../data/paleoclim/CHELSA_cur_V1_2B_r10m/10min/bio_1.tif")
pre_cur <- raster("../data/paleoclim/CHELSA_cur_V1_2B_r10m/10min/bio_12.tif")
mat_lgm <- raster("../data/paleoclim/chelsa_LGM_v1_2B_r10m/10min/bio_1.tif")
pre_lgm <- raster("../data/paleoclim/chelsa_LGM_v1_2B_r10m/10min/bio_12.tif")

extent(mat_cur); extent(pre_cur)
extent(mat_lgm); extent(pre_lgm)

mat_lgm_crop <- crop(mat_lgm, extent(-180.0001, 179.9999, -90.00014, 83.99984))
mat_lgm_ano <- calc(stack(mat_lgm_crop, mat_cur), diff)
pre_lgm_crop <- crop(pre_lgm, extent(-180.0001, 179.9999, -90.00014, 83.99984))
pre_lgm_ano <- calc(stack(pre_lgm_crop, mat_cur), diff)


# Average for botanical countries -------------------------------------------


shape <- readOGR("../data/shapefile_bot_countries/level3_mod.shp")
pclimate_vars <- c("mio_pre_ano", "mio_mat_ano", "mio_pre", "mio_mat", "pi_pre",
                   "mat_lgm_ano", "pre_lgm_ano")
pstat_vars <- c("mean", "sd", "n")
combs <- nrow(expand.grid(pclimate_vars, pstat_vars))
m <- matrix(seq(1:combs), ncol=3, byrow = TRUE)
num_list <- split(m, rep(1:nrow(m)))

pres <- matrix(nrow=nrow(shape@data), ncol=combs)
rownames(pres) <- shape@data$LEVEL_3_CO
disag_id <- c()
upsale_count <- c()

# make sure projections match
projection(pi_pre) <- projection(shape)

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
      upsale_count <- c(upsale_count, 1)
      disag_id <- c(disag_id, paste(i,j))
      # to avoid huge raster files crop to extent of shapefile sub * 20
      newExtent <- extent(bbox(shape_sub))
      lay2 <- crop(lay, newExtent*20)
      lay2 <- disaggregate(lay2, 10)
      rest <- raster::extract(lay2, shape_sub)
      rest <- na.omit(rest[[1]])
      
    }else{}
    
    # get mean, sd and sample size
    pres[i,num_list[[j]][1]] <- mean(rest)
    pres[i,num_list[[j]][2]] <- sd(rest)
    pres[i,num_list[[j]][3]] <- length(rest)
  }
  if(!i%%1)cat(i,"\r")
}

temp <- c(paste(pclimate_vars[1], pstat_vars, sep="_"),
          paste(pclimate_vars[2], pstat_vars, sep="_"),
          paste(pclimate_vars[3], pstat_vars, sep="_"),
          paste(pclimate_vars[4], pstat_vars, sep="_"),
          paste(pclimate_vars[5], pstat_vars, sep="_"),
          paste(pclimate_vars[6], pstat_vars, sep="_"),
          paste(pclimate_vars[7], pstat_vars, sep="_"))
pres <- as.data.frame(pres)
names(pres) <- temp

saveRDS(pres, "../processed_data/paleo_climate.rds")



