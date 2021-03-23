library(classInt)
library(rgdal)
library(raster)
#library(mmSAR) not available for R 3.5.2
library(rgeos)
library(foreign)

#set working directory
setwd("~/Documents/WOLF/PROJECTS/58 World Checklist paper/analyses december 2018")

#path to WCSP download
wcsp_path = "/Users/au265104/Documents/WOLF/PROJECTS/58 World Checklist paper/October 2018"

#load("WCSP_geo.RData")

#read WCSP
read.csv(paste(wcsp_path, "/published_names_19_10_2018.csv", sep=""), header=TRUE, sep="|") -> published
read.csv(paste(wcsp_path, "/unpublished_names_19_10_2018.csv", sep=""), header=TRUE, sep="|") -> unpublished
wcsp <- rbind(published,unpublished)
rm(published, unpublished)

#read distributions 
read.csv(paste(wcsp_path, "/published_distribution_19_10_2018.csv", sep=""), header=TRUE, sep="|") -> published
read.csv(paste(wcsp_path, "/unpublished_distribution_19_10_2018.csv", sep=""), header=TRUE, sep="|") -> unpublished
dist <- rbind(published,unpublished)
rm(published, unpublished)
dist = dist[dist$introduced==0 & dist$extinct==0,] #remove introduced and extinct records
dist$area_code_l3 = toupper(dist$area_code_l3)

shape = readOGR(dsn = "/Volumes/My Book/WOLF/Aarhus PC/C/GIS workspace/projektuebergreifend",, layer = "level3")

goodspp = wcsp[wcsp$genus_hybrid_marker == "" & wcsp$species_hybrid_marker == "" & wcsp$infraspecific_rank == "" & wcsp$species != "" & wcsp$taxon_status_description == "Accepted",]

# create empty community matrix
comm = matrix(nrow=nrow(goodspp), ncol=nrow(shape@data))
colnames(comm) = as.vector(shape@data$LEVEL_3_CO)
rownames(comm) = goodspp$checklist_id
comm[is.na(comm)] <- 0

# remove distribution records that cannot be assigned to proper TDWG units
dist = dist[dist$area_code_l3 %in% colnames(comm),]

# synonymise distributions table by adding accepted_plant_id column
accids = wcsp$accepted_name_id
names(accids) = wcsp$checklist_id
dist$accepted_name_id = accids[as.vector(dist$checklist_id)]
dist$accepted_name_id = as.vector(dist$accepted_name_id)
rm(accids)

# remove distribution records that have no accepted id
dist = dist[dist$accepted_name_id != "",]

# remove distribution records of infraspecific taxa
dist = dist[dist$accepted_name_id %in% rownames(comm),]

# indexing rownames for faster processing
idx = 1:nrow(comm)
names(idx) = rownames(comm)
dist$accepted_name_idx = idx[dist$accepted_name_id]

# populating community table
for(i in 1:nrow(dist)){
  comm[dist[i,"accepted_name_idx"], dist[i,"area_code_l3"]] = 1
  if(i%%16400==0) print("*")
}

# add species richness to polygons
SR = colSums(comm)
shape@data$SR = SR[shape@data$LEVEL_3_CO]

#calculate area
#behrmann projection
shape_behr = spTransform(shape, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"))
shape_behr@data$area=area(shape_behr)

#plot area for testing
mypal <- c("aliceblue", "red3")
colscale <- findColours(classIntervals(shape_behr@data$area, style="jenks"), mypal)
plot(shape_behr, col=colscale)

#calculate area correction
plot(shape_behr@data$area, shape_behr@data$SR)
plot(log(shape_behr@data$area), log(shape_behr@data$SR))

as = list(name="data.wcsp", data=data.frame(a=shape_behr@data$area, s=shape_behr@data$SR))
as$data = as$data[!is.na(as$data$s),]

data(power)
data(expo)
data(negexpo)
data(monod) 
data(ratio) 
data(logist) 
data(lomolino) 
data(weibull)
mods = list(power,expo,negexpo,logist,ratio,lomolino,weibull) #monod doesn't work, creates some kind of singularity

AICs = vector("numeric",length(mods))

for(i in 1:length(mods)){
  res = rssoptim(mods[[i]],as) 
  AICs[i] = res$AIC
}

AICs - min(AICs)

#best model
res = rssoptim(weibull,as)

#fitted values
a = shape_behr@data$area
fitted = res$par[1]*(1-exp(-1*res$par[2]*a^res$par[3]))

#plot fit
plot(shape_behr@data$area, shape_behr@data$SR,xlab="area", ylab="species richness")
temp = data.frame(a,fitted)
temp = temp[sort(temp$a, index.return=T)$ix,]
lines(temp$a, temp$fitted, col="red")
rm(temp)

#logged
plot(log(shape_behr@data$area), log(shape_behr@data$SR),xlab="log(area)", ylab="log(species richness)")
temp = data.frame(a,fitted)
temp = temp[sort(temp$a, index.return=T)$ix,]
lines(log(temp$a), log(temp$fitted), col="red")
rm(temp)

#calculate residuals (area corrected SR)
shape_behr@data$SR_res = shape_behr@data$SR-fitted

#plot raw richness
mypal <- c("aliceblue", "red3")
colscale <- findColours(classIntervals(shape_behr@data$SR, style="jenks"), mypal)
plot(shape_behr, col=colscale, main="Raw species richness", border=NA)

#plot area corrected richness
colscale <- findColours(classIntervals(shape_behr@data$SR_res, style="jenks"), mypal)
plot(shape_behr, col=colscale, main="Area corrected species richness", border=NA)

#import essential info from match_WCSP_phylo.R

#read root distances
RD_full <- read.table("RD_full.txt")

#read link table
link <- read.table("link.txt")
rownames(link) <- link$checklist_id

#calculate MRD
MRD = vector("numeric",ncol(comm))
names(MRD) = colnames(comm)
for(i in 1:ncol(comm)){
  tips = link[rownames(comm)[comm[,i]==1],2]
  tips = tips[!is.na(tips)]
  MRD[i] = mean(RD_full[tips,])
  print(i)
}

#20.1.2019

#add MRD to shapefile

# add species richness to polygons
shape_behr@data$MRD =  MRD[shape_behr@data$LEVEL_3_CO]

#plot MRD
colscale <- findColours(classIntervals(shape_behr@data$MRD, style="jenks"), mypal)
plot(shape_behr, col=colscale, main="Mean root distance")

plot(MRD, SR)
plot(shape_behr@data$area, shape_behr@data$MRD)

#read additional information on TDWG units from 2012 PNAS paper
PNAS_info <- read.csv("TDWG_level3_Realms_2010Dec_mod2019.csv", sep=";") #data for VNA added except Lat Lon
rownames(PNAS_info) <- PNAS_info$LEVEL_3_CO

#add island info to shapefile
shape_behr@data$island = PNAS_info[as.vector(shape_behr@data$LEVEL_3_CO),"ISISLAND"]

#rescale area in order to be able to write to shapefile
#shape_behr@data$area = shape_behr@data$area/1000000 

writeOGR(shape_behr, dsn = getwd(), layer = "shape_behr", driver = "ESRI Shapefile")

#21.2.2019

#calculate MedRD
MedRD = vector("numeric",ncol(comm))
names(MedRD) = colnames(comm)
for(i in 1:ncol(comm)){
  tips = link[rownames(comm)[comm[,i]==1],2]
  tips = tips[!is.na(tips)]
  MedRD[i] = median(RD_full[tips,])
  print(i)
}

hist(MedRD)

plot(MRD, MedRD)

#load predictors from Stabfor study

predpath = "/Users/au265104/Documents/WOLF/PROJECTS/35 STABFOR 3 Large-scale/~analysis/predictors"

pred <- as.data.frame(matrix(ncol=0,nrow=nrow(shape_behr@data)))
rownames(pred) <- shape_behr@data$LEVEL_3_CO

Current_AP <- read.dbf(paste(predpath, "/Current_AP.dbf", sep=""))
rownames(Current_AP) <- Current_AP[,"LEVEL_3_CO"]
pred$Current_AP <- Current_AP[rownames(pred),"MEAN"]

Current_MAT <- read.dbf(paste(predpath, "/Current_MAT.dbf", sep=""))
rownames(Current_MAT) <- Current_MAT[,"LEVEL_3_CO"]
pred$Current_MAT <- Current_MAT[rownames(pred),"MEAN"]

Elevation <- read.dbf(paste(predpath, "/Elevation.dbf", sep=""))
rownames(Elevation) <- Elevation[,"LEVEL_3_CO"]
pred$Elevation <- Elevation[rownames(pred),"RANGE"]

LGM_ano_AP <- read.dbf(paste(predpath, "/LGM_ano_AP.dbf", sep=""))
rownames(LGM_ano_AP) <- LGM_ano_AP[,"LEVEL_3_CO"]
pred$LGM_ano_AP <- LGM_ano_AP[rownames(pred),"MEAN"]

LGM_ano_MAT <- read.dbf(paste(predpath, "/LGM_ano_MAT.dbf", sep=""))
rownames(LGM_ano_MAT) <- LGM_ano_MAT[,"LEVEL_3_CO"]
pred$LGM_ano_MAT <- LGM_ano_MAT[rownames(pred),"MEAN"]

#LGM_velo_AP <- read.dbf(paste(predpath, "/LGM_velo_AP.dbf", sep=""))
#rownames(LGM_velo_AP) <- LGM_velo_AP[,"LEVEL_3_CO"]
#pred$LGM_velo_AP <- LGM_velo_AP[rownames(pred),"MEAN"]

#LGM_velo_MAT <- read.dbf(paste(predpath, "/LGM_velo_MAT.dbf", sep=""))
#rownames(LGM_velo_MAT) <- LGM_velo_MAT[,"LEVEL_3_CO"]
#pred$LGM_velo_MAT <- LGM_velo_MAT[rownames(pred),"MEAN"]

Mio_ano_AP <- read.dbf(paste(predpath, "/Mio_ano_AP.dbf", sep=""))
rownames(Mio_ano_AP) <- Mio_ano_AP[,"LEVEL_3_CO"]
pred$Mio_ano_AP <- Mio_ano_AP[rownames(pred),"MEAN"]

Mio_ano_MAT <- read.dbf(paste(predpath, "/Mio_ano_MAT.dbf", sep=""))
rownames(Mio_ano_MAT) <- Mio_ano_MAT[,"LEVEL_3_CO"]
pred$Mio_ano_MAT <- Mio_ano_MAT[rownames(pred),"MEAN"]

Plio_ano_AP <- read.dbf(paste(predpath, "/Plio_ano_AP.dbf", sep=""))
rownames(Plio_ano_AP) <- Plio_ano_AP[,"LEVEL_3_CO"]
pred$Plio_ano_AP <- Plio_ano_AP[rownames(pred),"MEAN"]

Plio_ano_MAT <- read.dbf(paste(predpath, "/Plio_ano_MAT.dbf", sep=""))
rownames(Plio_ano_MAT) <- Plio_ano_MAT[,"LEVEL_3_CO"]
pred$Plio_ano_MAT <- Plio_ano_MAT[rownames(pred),"MEAN"]

soiltypes <- read.dbf(paste(predpath, "/soiltypes.dbf", sep=""))
rownames(soiltypes) <- soiltypes[,"LEVEL_3_CO"]
pred$soiltypes <- soiltypes[rownames(pred),"VARIETY"]

climrar <- read.table(paste(predpath, "/Climate rarity country means.csv", sep=""),sep=",",header=T)
rownames(climrar) <- climrar[,"Country"]
pred$climrar <- climrar[rownames(pred),"ClimateRarity"]

#plot all vs all scatter plots
data = cbind(SR = shape_behr@data$SR, MRD = shape_behr@data$MRD, area = shape_behr@data$area,  pred)
data = data[!shape_behr@data$LEVEL_3_CO == "ANT",] #removing Antarctica
plot(data)

#correlations
write.csv(cor(data, use="pairwise"), "correlations.csv")

#exploratory modelling
mod_SR = lm(SR ~ log(area) + Current_AP + Current_MAT + Elevation + LGM_ano_AP + LGM_ano_MAT + 
           Mio_ano_AP + Mio_ano_MAT + Plio_ano_AP + Plio_ano_MAT  +  soiltypes  +  climrar, data=data)

mod_MRD = lm(MRD ~ log(area) + Current_AP + Current_MAT + Elevation + LGM_ano_AP + LGM_ano_MAT + 
              Mio_ano_AP + Mio_ano_MAT + Plio_ano_AP + Plio_ano_MAT  +  soiltypes  +  climrar, data=data)

# 22.2.2019 + 7.3.2019

# data transformations and standardisation

#data <- data[rownames(data)!="GNL",] #exclude Greenland, which is an extreme outlier in SAR

data$logSR = log(data$SR)
data_no0 <- data[data$SR>0,]
data_no0$sqrtAP <- sqrt(data_no0$Current_AP)
data_no0$sqrtElevation <- sqrt(data_no0$Elevation)
data_no0$logarea <- log(data_no0$area)

SR_res = shape_behr@data$SR_res
names(SR_res) <- shape_behr@data$LEVEL_3_CO
data_no0$SR_res <- SR_res[rownames(data_no0)]
rm(SR_res)

data_no0$logSR_res = data_no0$SR_res+1-min(data_no0$SR_res)
data_no0$logSR_res = log(data_no0$logSR_res)

data_no0_std <- data_no0
for(i in 1:ncol(data_no0_std)){
  data_no0_std[,i] <- (data_no0_std[,i]-mean(data_no0_std[,i], na.rm = TRUE))/sd(data_no0_std[,i], na.rm = TRUE)
}

data_no0_std <- data_no0_std[,!colnames(data_no0_std) %in% c("Current_AP", "Elevation", "area")]

#summarising key heterogeneity variables

data_heterog <- data_no0_std[,c("logarea", "sqrtElevation", "soiltypes")]

pca_heterog <- prcomp(~logarea+sqrtElevation+soiltypes, data=data_heterog, na.action = na.omit, center = TRUE, scale = TRUE)

summary(pca_heterog)

data_no0_std$pca_heterog = rep(NA,nrow(data_no0_std))

data_no0_std[rownames(pca_heterog$x),"pca_heterog"] = -1*pca_heterog$x[,"PC1"] 

# ... and without area (for analyses of area-corrected richness)
data_heterog1 <- data_no0_std[,c("sqrtElevation", "soiltypes")]

pca_heterog1 <- prcomp(~sqrtElevation+soiltypes, data=data_heterog1, na.action = na.omit, center = TRUE, scale = TRUE)

summary(pca_heterog)

data_no0_std$pca_heterog1 = rep(NA,nrow(data_no0_std))

data_no0_std[rownames(pca_heterog$x),"pca_heterog1"] = -1*pca_heterog1$x[,"PC1"] 

# testing for collinearity

library(car)

vif(lm(SR~., data=data_no0_std))

kappa(data_no0_std[!is.na(rowSums(data_no0_std)),-1])

# playing with SEM

library(lavaan)

# model <- ' 
#   # measurement model
#     f1 =~ MRD
#   # regressions
#     logSR ~ sqrtAP + Current_MAT + sqrtElevation + soiltypes + logarea + f1
#     f1 ~ sqrtAP + Current_MAT + sqrtElevation + soiltypes + logarea
#   # residual correlations
#     sqrtAP ~~ Current_MAT 
#     sqrtAP ~~ sqrtElevation
#     sqrtAP ~~ soiltypes
#     sqrtAP ~~ logarea
#     Current_MAT ~~ sqrtElevation
#     Current_MAT ~~ soiltypes
#     Current_MAT ~~ logarea
#     sqrtElevation ~~ soiltypes
#     sqrtElevation ~~ logarea
#     soiltypes ~~ logarea
#   '

model <- ' 
  # measurement model
    f1 =~ MRD
  # regressions
    logSR ~ sqrtAP + Current_MAT + pca_heterog + f1
    f1 ~ sqrtAP + Current_MAT + pca_heterog
  # residual correlations
    sqrtAP ~~ Current_MAT 
    sqrtAP ~~ pca_heterog
    Current_MAT ~~ pca_heterog
'

#with area standardized SR
model <- ' 
  # measurement model
    f1 =~ MRD
  # regressions
    SR_res ~ sqrtAP + Current_MAT + pca_heterog1 + f1
    f1 ~ sqrtAP + Current_MAT + pca_heterog1
  # residual correlations
    sqrtAP ~~ Current_MAT 
    sqrtAP ~~ pca_heterog1
    Current_MAT ~~ pca_heterog1
'

fit <- sem(model, data=data_no0_std)
summary(fit, standardized=TRUE, rsquare = TRUE)

#plots

plot(data_no0_std[,c("SR_res","MRD","pca_heterog1","sqrtAP","Current_MAT")], pch=15, col="grey")




#save.image(file = "WCSP_geo.RData")

save(dist,file="dist.RData")
