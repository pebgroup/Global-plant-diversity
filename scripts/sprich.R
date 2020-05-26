setwd("F:/PROJECTS/38 Tansley review/analysis")

# Load ecoregion data from Kier et al. 2005 J Biogeogr 32, 1107-1116
# Downloaded from http://databasin.org/datasets/43478f840ac84173979b22631c2ed672
# Processed in ArcGIS in G:\WOLF\Aarhus PC\C\GIS workspace\38 Tansley
# This dataset contains multiple polygons per ecoregion
sprich = read.table("sprich2.txt", sep=",", header=T, dec=".", stringsAsFactors=F)
sprich$AREA = as.numeric(sprich$AREA)
sprich$BIOME = as.numeric(sprich$BIOME)
sprich$plant_spcs = as.numeric(sprich$plant_spcs)

# Eliminate invalid ecoregions
sprich = sprich[sprich$BIOME<50,]
sprich = sprich[sprich$plant_spcs>0,]

# Merge the different polygons of each ecoregion
sprich_red = cbind(aggregate(sprich$AREA, by=list(sprich$ECO_CODE), FUN=sum),
                   BIOME=aggregate(sprich$BIOME, by=list(sprich$ECO_CODE), FUN=max)$x,
                   plant_spcs=aggregate(sprich$plant_spcs, by=list(sprich$ECO_CODE), FUN=max)$x)

# Calculate log(area) and log(species richness) for area correction
# Area is in square meters
sprich_red$logAREA = log(sprich_red$x)
sprich_red$logplant_spcs = log(sprich_red$plant_spcs)

# Check species-area relationship
plot(logplant_spcs~logAREA, data = sprich_red)

# Global standardisation to 10000 km2
mod = lm(logplant_spcs~logAREA, data = sprich_red)
sprich_red$plant_spcs_std_glob = sprich_red[,"plant_spcs"] * (10000000000/sprich_red[,"x"])^coefficients(mod)[2]

# Per-biome standardisation to 10000 km2
STD_BIOM = NULL
for(r in 1:14){
  sprich_temp = sprich_red[sprich_red$BIOME==r,]
  mod = lm(logplant_spcs~logAREA, data = sprich_temp)
  std_biom = sprich_temp[,"plant_spcs"] * (10000000000/sprich_temp[,"x"])^coefficients(mod)[2]
  names(std_biom) = sprich_temp$Group.1
  STD_BIOM = c(STD_BIOM, std_biom)
}
sprich_red$plant_spcs_std_biom = STD_BIOM[sprich_red$Group.1]

# Compare the two standardizations
plot(plant_spcs_std_biom ~ plant_spcs_std_glob, data = sprich_red)

# Plot area-standardized richness (standardized within each biome) against biome
boxplot(plant_spcs_std_biom~BIOME, data = sprich_red)
abline(h=mean(sprich_red$plant_spcs_std_biom), col="red")

# Plot area-standardized richness (standardized within each biome) for one biome against all others
boxplot(plant_spcs_std_biom~(sprich_red$BIOME == 3), data = sprich_red)


sprich_red[sort(sprich_red$plant_spcs_std_biom, decreasing=T, index.return=T)$ix,]




plot(logplant_spcs~logAREA, data = sprich_red[sprich_red$BIOME==2,])
mod = lm(logplant_spcs~logAREA, data = sprich_red[sprich_red$BIOME==2,])

abline(mod, col="red")

r = sample(rownames(sprich_red[sprich_red$BIOME==2,]),1)

points(logplant_spcs~logAREA, data = sprich_red[r,], col="red") 

points(log(10000000000), log(sprich_red[r,"plant_spcs"] * (10000000000/sprich_red[r,"x"])^coefficients(mod)[2]), col="red", pch=8)

sprich_red[,"plant_spcs"] * (10000000000/sprich_red[,"x"])^coefficients(mod)[2]


hist(sprich_red$x)

acorr = lm(plant_spcs~AREA)

boxplot(acorr$residuals~BIOME)

boxplot(plant_spcs~BIOME)

str(sprich)

sprich1 = sprich[BIOME==1,]

levels(sprich$BIOME)

lm(plant_spcs~as.factor(BIOME)+AREA)

sprich_red = cbind(aggregate(sprich$AREA, by=list(sprich$ECO_CODE), FUN=sum),
BIOME=aggregate(sprich$BIOME, by=list(sprich$ECO_CODE), FUN=max)$x,
plant_spcs=aggregate(sprich$plant_spcs, by=list(sprich$ECO_CODE), FUN=max)$x)

boxplot(plant_spcs~BIOME, data=sprich_red)
abline(h=mean(sprich_red$plant_spcs), col="red")

plot(plant_spcs~x, data = sprich_red)
