#load pa matrix
pa <- read.table("../02 extract WCSPF data/checklist_pa.csv",sep=";",header=T)
rownames(pa) <- pa[,1]
pa <- pa[,2:ncol(pa)]

rownames(pa)[rownames(pa)== "Cunninghamia_lanceolata_var__konishii"] <- "Cunninghamia_lanceolata_var_konishii"
rownames(pa)[rownames(pa)== "Pinus_caribaea_var__hondurensis"] <- "Pinus_caribaea_var_hondurensis"
rownames(pa)[rownames(pa)== "Pinus_caribaea_var__bahamensis"] <- "Pinus_caribaea_var_bahamensis"
rownames(pa)[rownames(pa)== "Pinus_caribaea_var__caribaea"] <- "Pinus_caribaea_var_caribaea"

#eliminate hybrids
temp <- strsplit(rownames(pa),"_")
hybrid <- c()
for(i in temp){
  if(i[2] == "x") hybrid <- c(hybrid,1) else hybrid <- c(hybrid,0)
}
pa <- pa[hybrid == 0,]
rm(temp, i, hybrid)

#map genera
temp <- strsplit(rownames(pa),"_")
genus <- c()
for(i in temp){
  genus <- c(genus,i[1])
}
genus <- as.data.frame(cbind(genus,rownames(pa)))
rm(temp, i)

# map families
fam <- read.table("genera-families.csv",sep=";")
rownames(fam) <- fam[,1]
genus <- data.frame(fam=fam[genus[,1],2],genus)
rm(fam)

# map JCs major clades
genus <- data.frame(group=numeric(nrow(genus)),genus)
genus[genus[,2]=="Pinaceae",1] <- 1
genus[genus[,2]=="Podocarpaceae" | genus[,2]=="Araucariaceae",1] <- 2
genus[genus[,2]=="Cupressaceae" | genus[,2]=="Taxaceae" | genus[,2]=="Sciadopityaceae",1] <- 3

# realms
realm_we2 <- read.table("realm_we2.txt",sep=";",header=T)
rownames(realm_we2) <- realm_we2[,"LEVEL_3_CO"]

colmap <- c("#038e3b", "#9966cc", "#999999", "#ff9933", "#03428e", "#ff0000")
names(colmap) <- c("AFRICA", "ASIA-TROPICAL", "AUSTRALASIA", "NORTHERN AMERICA", "PALAEARCTIC", "SOUTHERN AMERICA")  

lev <- c("AFRICA" , "ASIA-TROPICAL" , "AUSTRALASIA" , "NORTHERN AMERICA" , "PALAEARCTIC" , "SOUTHERN AMERICA")

# load coordinates
TDWG <- read.table("C:/GIS workspace/STABFOR 3/TDWG_level3_Realms_2010Dec_modified from Daniel.csv",sep=";",header=T)
rownames(TDWG) <- TDWG[,1]

#load indices
load("data/Ic_gamma_clade1.Rdata")
load("data/Ic_gamma_clade2.Rdata")
load("data/Ic_gamma_clade3.Rdata")
load("data/NRI_NTI_clade1.Rdata")
load("data/NRI_NTI_clade2.Rdata")
load("data/NRI_NTI_clade3.Rdata")

data.clade1 <- data.frame(NRI = apply(NRI.world.clade1,1,mean)[rownames(Ic.world.clade1)],NTI = apply(NTI.world.clade1,1,mean)[rownames(Ic.world.clade1)],Ic = apply(Ic.world.clade1,1,mean),gamma=apply(gamma.world.clade1,1,mean))
data.clade2 <- data.frame(NRI = apply(NRI.world.clade2,1,mean)[rownames(Ic.world.clade2)],NTI = apply(NTI.world.clade2,1,mean)[rownames(Ic.world.clade2)],Ic = apply(Ic.world.clade2,1,mean),gamma=apply(gamma.world.clade2,1,mean))
data.clade3 <- data.frame(NRI = apply(NRI.world.clade3,1,mean)[rownames(Ic.world.clade3)],NTI = apply(NTI.world.clade3,1,mean)[rownames(Ic.world.clade3)],Ic = apply(Ic.world.clade3,1,mean),gamma=apply(gamma.world.clade3,1,mean))

data <- list()
data[[1]] <- data.clade1
data[[2]] <- data.clade2
data[[3]] <- data.clade3

# load global PCS

load("data/globalPCS20130816.RData")

for(i in 1:3){
  colnames(RES[[i]]) <- c("global.NRI", "global.NTI", "global.Ic", "global.gamma")
  data[[i]] <- cbind(data[[i]], RES[[i]][rownames(data[[i]]),])
  
}

# load data

pred <- as.data.frame(matrix(ncol=0,nrow=ncol(pa)))
rownames(pred) <- colnames(pa)

Current_AP <- read.dbf("predictors/Current_AP.dbf")
rownames(Current_AP) <- Current_AP[,"LEVEL_3_CO"]
pred$Current_AP <- Current_AP[rownames(pred),"MEAN"]

Current_MAT <- read.dbf("predictors/Current_MAT.dbf")
rownames(Current_MAT) <- Current_MAT[,"LEVEL_3_CO"]
pred$Current_MAT <- Current_MAT[rownames(pred),"MEAN"]

Elevation <- read.dbf("predictors/Elevation.dbf")
rownames(Elevation) <- Elevation[,"LEVEL_3_CO"]
pred$Elevation <- Elevation[rownames(pred),"RANGE"]

LGM_ano_AP <- read.dbf("predictors/LGM_ano_AP.dbf")
rownames(LGM_ano_AP) <- LGM_ano_AP[,"LEVEL_3_CO"]
pred$LGM_ano_AP <- LGM_ano_AP[rownames(pred),"MEAN"]

LGM_ano_MAT <- read.dbf("predictors/LGM_ano_MAT.dbf")
rownames(LGM_ano_MAT) <- LGM_ano_MAT[,"LEVEL_3_CO"]
pred$LGM_ano_MAT <- LGM_ano_MAT[rownames(pred),"MEAN"]

#LGM_velo_AP <- read.dbf("predictors/LGM_velo_AP.dbf")
#rownames(LGM_velo_AP) <- LGM_velo_AP[,"LEVEL_3_CO"]
#pred$LGM_velo_AP <- LGM_velo_AP[rownames(pred),"MEAN"]

#LGM_velo_MAT <- read.dbf("predictors/LGM_velo_MAT.dbf")
#rownames(LGM_velo_MAT) <- LGM_velo_MAT[,"LEVEL_3_CO"]
#pred$LGM_velo_MAT <- LGM_velo_MAT[rownames(pred),"MEAN"]

Mio_ano_AP <- read.dbf("predictors/Mio_ano_AP.dbf")
rownames(Mio_ano_AP) <- Mio_ano_AP[,"LEVEL_3_CO"]
pred$Mio_ano_AP <- Mio_ano_AP[rownames(pred),"MEAN"]

Mio_ano_MAT <- read.dbf("predictors/Mio_ano_MAT.dbf")
rownames(Mio_ano_MAT) <- Mio_ano_MAT[,"LEVEL_3_CO"]
pred$Mio_ano_MAT <- Mio_ano_MAT[rownames(pred),"MEAN"]

Plio_ano_AP <- read.dbf("predictors/Plio_ano_AP.dbf")
rownames(Plio_ano_AP) <- Plio_ano_AP[,"LEVEL_3_CO"]
pred$Plio_ano_AP <- Plio_ano_AP[rownames(pred),"MEAN"]

Plio_ano_MAT <- read.dbf("predictors/Plio_ano_MAT.dbf")
rownames(Plio_ano_MAT) <- Plio_ano_MAT[,"LEVEL_3_CO"]
pred$Plio_ano_MAT <- Plio_ano_MAT[rownames(pred),"MEAN"]

soiltypes <- read.dbf("predictors/soiltypes.dbf")
rownames(soiltypes) <- soiltypes[,"LEVEL_3_CO"]
pred$soiltypes <- soiltypes[rownames(pred),"VARIETY"]

climrar <- read.table("predictors/Climate rarity country means.csv",sep=",",header=T)
rownames(climrar) <- climrar[,"Country"]
pred$climrar <- climrar[rownames(pred),"ClimateRarity"]

# AET <- read.dbf("predictors/AET.dbf")
# rownames(AET) <- AET[,"LEVEL_3_CO"]
# pred$AET <- AET[rownames(pred),"MEAN"]
# 

# load climate suitability

global_raw <- read.table("suitability/global_raw_poly.txt", sep=" ", header=T)
rownames(global_raw) <- global_raw$LVL3
suit <- data.frame(current.suit = global_raw[,"Current"], min.suit = global_raw[,"min.suit"], mean.suit = global_raw[,"mean.suit"])
rownames(suit) <- rownames(global_raw)

# clade1_raw <- read.table("suitability/clade1_raw.txt", sep=" ", header=T)
# rownames(clade1_raw) <- clade1_raw$LVL3
# suit <- cbind(suit, clade1.min.suit = clade1_raw[rownames(suit),"min.suit"], clade1.mean.suit = clade1_raw[rownames(suit),"mean.suit"])
# 
# clade2_raw <- read.table("suitability/clade2_raw.txt", sep=" ", header=T)
# rownames(clade2_raw) <- clade2_raw$LVL3
# suit <- cbind(suit, clade2.min.suit = clade2_raw[rownames(suit),"min.suit"], clade2.mean.suit = clade2_raw[rownames(suit),"mean.suit"])
# 
# clade3_raw <- read.table("suitability/clade3_raw.txt", sep=" ", header=T)
# rownames(clade3_raw) <- clade3_raw$LVL3
# suit <- cbind(suit, clade3.min.suit = clade3_raw[rownames(suit),"min.suit"], clade3.mean.suit = clade3_raw[rownames(suit),"mean.suit"])

# current.suit <- data.frame(global = global_raw[,"Current"], clade1 = clade1_raw[rownames(global_raw),"Current"], clade2 = clade2_raw[rownames(global_raw),"Current"], clade3 = clade3_raw[rownames(global_raw),"Current"])
# rownames(current.suit) <- rownames(global_raw)

# clear redundant vars
rm(list=ls()[!ls() %in% c("data", "lev", "pa", "genus", "pred", "realm_we2", "colmap", "TDWG", "suit", "current.suit")])
 


