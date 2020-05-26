# check Pinto-Sanchez et al. 2014 for their species-area connection
# 
#
# SPECIES RICHNESS - AREA CONNECTION
area <- read.csv("pinto-sanchez_SR_area_data.csv")
# exclude areas they excluded
area <- area[!area$too_big=="yes",]

hist(area$Area..hectares., breaks=20)
hist(area$Species.richness)

# they state 39 sites with relevant data: I get only 36 sites with data, rest is NA
library(ggplot2)
theme_set(theme_classic())
ggplot(data= area, aes(Area..hectares.,Species.richness))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

cor.test(area$Area..hectares., area$Species.richness, method="s")
# rho=-0.06, p=0.71


# SPECIES RICHNESS - SAMPLING EFFORT
eff <- read.csv("pinto-sanchez_SR_area_data_search_hours.csv")
hist(eff$Total.sampling.effort..hours.invested.)

# NOTE: ref for the third excluded area is not right in one of the sets - differs
eff <- eff[!eff$too_big=="yes",]

dat <- merge(area, eff, all.x=TRUE)

ggplot(data= dat, aes(Species.richness, Total.sampling.effort..hours.invested.))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

cor.test(dat$Species.richness, dat$Total.sampling.effort..hours.invested., method="s")
# rho=0.10, p=0.54
# 


 
