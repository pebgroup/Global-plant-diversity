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