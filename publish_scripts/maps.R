# map figures

# SR MAP ##########################################################################

temp <- readRDS("data/data_for_SEM.rds")
temp2 <- shp[shp$LEVEL_3_CO %in% rownames(temp),]

thicc_lines <-temp2[which(temp2$area<1200000000),]
lcol <- min(thicc_lines$sr)/max(temp2$sr)
ucol <- max(thicc_lines$sr)/max(temp2$sr)

(sr_map <- ggplot(temp2) + 
    geom_sf(aes(fill=sr),lwd=0.1) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
    scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
                           begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("SR", option = "plasma", trans = "sqrt")+ #, 
    #  scale_fill_gradient("Species richness", guide="colorbar", low="white", high="darkred",
    #                      limits=c(min(sr_df$species_richness),max(sr_df$species_richness)))+
    #guides(fill = guide_colourbar(barwidth = 20, direction="vertical"))+ # stretch that colorbar
    theme_void()+
    theme(legend.position = c(0.21, 0.3), 
          legend.background = element_blank(),
          legend.key = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"))
)
ggsave(sr_map, filename=paste0("results/sr_map_red_", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=10, height=7)


# MRD MAP ##########################################################################
lcol.mrd <- (min(thicc_lines$mrd)-min(temp2$mrd))/(max(temp2$mrd)-min(temp2$mrd))
ucol.mrd <- (max(thicc_lines$mrd)-min(temp2$mrd))/(max(temp2$mrd)-min(temp2$mrd))
ggplot(temp2) + 
  geom_sf(aes(fill = mrd), lwd=0.1) + 
  geom_sf(data=thicc_lines, lwd=1.5, aes(col=mrd), show.legend=F)+
  scale_colour_viridis_c("SR", option = "plasma", 
                         begin = lcol.mrd, end = ucol.mrd)+
  scale_fill_viridis_c(option = "plasma")+ #, 
  theme_void()+
  theme(legend.position = c(0.21, 0.3), 
        legend.background = element_blank(),
        legend.key = element_blank())
ggsave(filename=paste0("results/mrd_map_red_", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=10, height=7)






## uncertainty maps ###########################################################################3
grid.arrange(ncol=1,
             ggplot(shp) + 
               geom_sf(aes(fill = MRD.sd), lwd=0.1) + 
               scale_fill_viridis_c("MRD sd", option = "plasma")+ #, 
               #  guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
               theme_void()
             ,
             ggplot(shp) + 
               geom_sf(aes(fill = MDR.sd), lwd=0.1) + 
               scale_fill_viridis_c("Equal splits SD", option = "plasma")+ #, 
               # guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
               theme_void()
)
cor.test(shp$MRD.sd, shp$MDR.sd)




range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}



## LDG stuff reduced data that has been used in SEMs ####
#library(scales)
shpbins <- temp2
shpbins$sr_scaled <- range01(shpbins$sr)
shpbins$mrd_scaled <- range01(shpbins$mrd)
range(shpbins$mrd_scaled, na.rm=T)
range(shpbins$sr_scaled, na.rm=T)
shpbins <- st_drop_geometry(shpbins)
temp <- pivot_longer(shpbins[,c("lat","mrd_scaled", "sr_scaled", "CONTINENT")], cols = c("mrd_scaled", "sr_scaled"))
temp <- as.data.frame(temp)

(continent_scatterplot <- ggplot(temp, aes(lat, value, col=name))+
    geom_point(alpha=0.5)+
    scale_y_continuous("scaled MRD and SR")+#,
    #sec.axis = sec_axis(~.*1, name="mean root distance"))+ 
    scale_x_continuous("latitude centroid botanical country")+
    facet_wrap(~CONTINENT, nrow = 3, scale="free_x")+ #, scale="free"
    geom_smooth(data=temp[temp$lat<0,], method="lm", se=F)+
    geom_smooth(data=temp[temp$lat>0,], method="lm", se=F)+
    #geom_smooth(method="loess", se=F)+
    scale_colour_startrek(labels=c("MRD", "SR"), name="")+
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8))
)
ggsave("scatterplot_SR_MRD_latitude_continents_red.png", dpi=600, width=6, height=4)

(global_scatterplot_red <- ggplot(temp2, aes(lat, sr))+
    geom_point(alpha=0.5)+
    scale_y_continuous("species richness", trans = "sqrt")+
    scale_x_continuous("latitude centroid botanical country")+
    geom_smooth(data=shp[shp$lat<0,], method="lm", col="red")+
    geom_smooth(data=shp[shp$lat>0,], method="lm", col="blue")
)
(global_scatterplot_mrd_red <- ggplot(temp2, aes(lat, mrd))+
    geom_point(alpha=0.5)+
    scale_y_continuous("mean root distance")+
    scale_x_continuous("latitude centroid botanical country")+    
    geom_smooth(data=shp[shp$lat<0,], method="lm", col="red")+
    geom_smooth(data=shp[shp$lat>0,], method="lm", col="blue")
)
(global_scatterplot_MRD_SR_red <- ggplot(temp2, aes(mrd, sr))+
    geom_point(alpha=0.5)+
    scale_y_continuous("species richness", trans = "sqrt")+
    scale_x_continuous("mean root distance")+
    geom_smooth(method="lm", col="black")
)
cor.test(temp2$sr, temp2$mrd)
grid.arrange(nrow=1, global_scatterplot_red, global_scatterplot_mrd_red, global_scatterplot_MRD_SR_red)


# All in one map SR #####
library(ggthemes)
library(ggpubr)
(sr_map2 <- ggplot(temp2) + 
    geom_sf(aes(fill=sr),lwd=0.1) + 
    geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
    scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
                           begin = lcol, end = sqrt(ucol))+
    scale_fill_viridis_c("SR", option = "plasma", trans = "sqrt")+ #, 
    theme(legend.position = c(0.21, 0.3),
          legend.background = element_blank(),
          legend.key = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm")
    )+
    coord_sf(xlim = c(-180, 180), ylim = c(-70, 85), expand = F)+
    xlab(" ")
)
(sr_ldg_map <- ggplot(temp2, aes(lat,sr))+
    geom_point(alpha=0.5)+
    scale_y_continuous("SR/1000", trans = "sqrt", 
                       labels = c("0", "5", "10", "15", "20"),
                       breaks = c(0, 5000, 10000, 15000, 20000))+
    scale_x_continuous("", labels = c("40°S","20°S", "0°", "20°N","40°N","60°N", "80°N"),
                       breaks = c(-40, -20, 0, 20,40,60, 80), 
                       position="top")+
    geom_smooth(se=T)+
    #geom_smooth(data=shp[shp$lat<0,], method="lm", col="black")+
    #geom_smooth(data=shp[shp$lat>0,], method="lm", col="black")+
    coord_flip(xlim =c(-70, 85), ylim =c(0, max(temp2$sr)+500),expand = F)+
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.margin = margin(0, 0, 0, 0.2, "cm"))
)
ggarrange(sr_map2, sr_ldg_map, 
          #labels = c("A", "B"),
          ncol = 2, nrow = 1, 
          widths = c(4,1))
ggsave("sr_map_AIO.png", dpi=600, width=10, height=3.7)





# All in one map MRD ####
lcol.mrd <- (min(thicc_lines$mrd)-min(temp2$mrd))/(max(temp2$mrd)-min(temp2$mrd))
ucol.mrd <- (max(thicc_lines$mrd)-min(temp2$mrd))/(max(temp2$mrd)-min(temp2$mrd))
mrd_map2 <- ggplot(temp2) + 
  geom_sf(aes(fill = mrd), lwd=0.1) + 
  geom_sf(data=thicc_lines, lwd=1.5, aes(col=mrd), show.legend=F)+
  scale_colour_viridis_c("SR", option = "plasma", 
                         begin = lcol.mrd, end = ucol.mrd)+
  scale_fill_viridis_c(option = "plasma")+ #, 
  theme(legend.position = c(0.21, 0.3),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm")
  )+
  coord_sf(xlim = c(-180, 180), ylim = c(-70, 85), expand = F)+
  xlab(" ")
(mrd_ldg_map <- ggplot(temp2, aes(lat,mrd))+
    geom_point(alpha=0.5)+
    scale_y_continuous("MRD")+
    scale_x_continuous("", labels = c("40°S","20°S", "0°", "20°N","40°N","60°N", "80°N"),
                       breaks = c(-40, -20, 0, 20,40,60, 80), 
                       position="top")+
    geom_smooth(se=T)+
    #geom_smooth(data=shp[shp$lat<0,], method="lm", col="black")+
    #geom_smooth(data=shp[shp$lat>0,], method="lm", col="black")+
    coord_flip(xlim =c(-70, 85),ylim =c(min(temp2$mrd)-1, max(temp2$mrd)+1),expand = F)+
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.margin = margin(0, 0, 0, 0.2, "cm"))
)
ggarrange(mrd_map2, mrd_ldg_map, 
          #labels = c("A", "B"),
          ncol = 2, nrow = 1, 
          widths = c(4,1))
ggsave("mrd_map_AIO.png", dpi=600, width=10, height=3.7)



# goal: coord flip + group americas, asias + europe+africa
temp$g3 <- rep(NA, nrow(temp))
temp$g3[grep("AMERICA", temp$CONTINENT)] <- "Americas"
temp$g3[grep("ASIA", temp$CONTINENT)] <- "Australasia"
temp$g3[grep("AFRICA|EUROPE", temp$CONTINENT)] <- "Afrope"
temp$g3 <- factor(temp$g3, levels = c("Americas","Afrope","Australasia"))
t2 <- subset(temp2, CONTINENT=="NORTHERN AMERICA")
(sr_map_bw <- ggplot(temp2[temp2$CONTINENT %in% c("NORTHERN AMERICA", "SOUTHERN AMERICA"),]) + 
    geom_sf(lwd=0.1, fill="white") + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm")
    )+
    coord_sf(xlim = c(-180, -10), ylim = c(-70, 85), expand = F)+
    xlab(" ")
)
ggsave("map_bw_americas.png", width=6, height=6)
library(ggimage)
#library(ggpubr)
img <- png::readPNG("map_bw.png")
(c1 <- ggplot(temp[temp$g3=="Americas",], aes(lat, value, col=name))+
    geom_point(alpha=0.5)+
    scale_y_continuous("scaled MRD and SR")+#,
    scale_x_continuous("latitude")+
    # facet_wrap(~g3, ncol = 3)+ #, scale="free"
    geom_smooth(method="loess", se=T)+
    coord_flip(ylim = c(0,1))+
    scale_colour_startrek(labels=c("MRD", "SR"), name="")+
    theme_void()+
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8),
          #panel.grid.major.y =element_line(),
          legend.position = "none",
          plot.margin = margin(1, 1, 2, 1, "cm"))
)
ggbackground(c1,"map_bw_americas.png")
(c2 <- ggplot(temp[temp$g3=="Afrope",], aes(lat, value, col=name))+
    geom_point(alpha=0.5)+
    scale_y_continuous("scaled MRD and SR")+#,
    scale_x_continuous("latitude")+
    # facet_wrap(~g3, ncol = 3)+ #, scale="free"
    geom_smooth(method="loess", se=T)+
    coord_flip(ylim = c(0,1))+
    scale_colour_startrek(labels=c("MRD", "SR"), name="")+
    theme_void()+
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8),
          #panel.grid = element_blank(),
          legend.position = "none")
)
(c3 <- ggplot(temp[temp$g3=="Australasia",], aes(lat, value, col=name))+
    geom_point(alpha=0.5)+
    scale_y_continuous("scaled MRD and SR")+#,
    scale_x_continuous("latitude")+
    # facet_wrap(~g3, ncol = 3)+ #, scale="free"
    geom_smooth(method="loess", se=T)+
    coord_flip(ylim = c(0,1))+
    scale_colour_startrek(labels=c("MRD", "SR"), name="")+
    theme_void()+
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8),
          panel.background =  element_blank(),
          legend.position = "none")
)

p1 <- ggarrange(c1,c2,c3, nrow=1)
ggbackground(p1, "map_bw.png")
ggsave("map_with_gradients.png", width=6, height=3)
#ggbackground(continent_scatterplot_comb, img)
(continent_scatterplot_comb <- ggplot(na.omit(temp), aes(lat, value, col=name))+
    geom_point(alpha=0.5)+
    scale_y_continuous("scaled MRD and SR")+#,
    scale_x_continuous("latitude")+
    facet_wrap(~g3, ncol = 3)+ #, scale="free"
    geom_smooth(method="loess", se=T)+
    coord_flip(ylim = c(0,1))+
    scale_colour_startrek(labels=c("MRD", "SR"), name="")+
    #theme_void()+
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8))
)

#ggsave("scatterplot_SR_MRD_latitude_continents_red.png", dpi=600, width=6, height=4)



# area + lat ####
# discovering the area neg MRD connection
(global_scatterplot_lat_area_red <- ggplot(temp2, aes(area, lat))+
   geom_point(alpha=0.5)+
   scale_y_continuous("lat")+
   scale_x_continuous("area", trans="sqrt")+
   geom_smooth(method="lm", col="black")
)
cor.test(temp2$lat, temp2$area)

#, label = unit_format(unit = "K"), limits=c(0,max(shp$sr)/1000)
#  geom_smooth(data=shp[shp$lat<0,], method="lm", col="red")+
#  geom_smooth(data=shp[shp$lat>0,], method="lm", col="blue")
#geom_text(aes(label=LEVEL_3_CO, col=factor(mrd_z_score_binary)),hjust=0, vjust=0)+
cor.test(abs(temp2$lat), temp2$sr, method = "p") # pearson correlation = -0.38
cor.test(temp2$lat[temp2$lat<0], temp2$sr[temp2$lat<0], method = "p")
cor.test(temp2$lat[temp2$lat>0], temp2$sr[temp2$lat>0], method = "p")
# summary(lm(sr~lat, data=shp[shp$lat<0,]))
# summary(lm(sr~lat, data=shp[shp$lat>0,]))
# summary(lm(sr~abs(lat), data=shp))

cor.test(abs(temp2$lat), temp2$mrd, method = "p")
cor.test(temp2$lat[temp2$lat<0], temp2$mrd[temp2$lat<0], method = "p")
cor.test(temp2$lat[temp2$lat>0], temp2$mrd[temp2$lat>0], method = "p")
summary(lm(mrd~lat, data=temp2[temp2$lat<0,]))
summary(lm(mrd~lat, data=temp2[temp2$lat>0,]))
summary(lm(mrd~abs(lat), data=temp2))


# Robustness of lat patterns #####
names(temp2)
table(temp2$lat<0) # much more data points in the northern hemisphere
reps <- 500
samp.size <- 50
res.hemisphere <- data.frame(cor.nh.mrd=NA, pval.nh.mrd=NA,
                             cor.sh.mrd=NA, pval.sh.mrd=NA,
                             cor.nh.sr=NA, pval.nh.sr=NA,
                             cor.sh.sr=NA, pval.sh.sr=NA,
                             n=c(1:reps))
set.seed(6348957)
for(i in 1:reps){
  nh <- temp2[temp2$lat>0,]
  temp <- nh[sample(1:nrow(nh), samp.size, replace = F),]
  res.hemisphere$cor.nh.sr[i] <- cor.test(temp$lat, temp$sr)$estimate
  res.hemisphere$pval.nh.sr[i] <- cor.test(temp$lat, temp$sr)$p.value
  res.hemisphere$cor.nh.mrd[i] <- cor.test(temp$lat, temp$mrd)$estimate
  res.hemisphere$pval.nh.mrd[i] <- cor.test(temp$lat, temp$mrd)$p.value
  
  sh <- temp2[temp2$lat<0,]
  temp <- sh[sample(1:nrow(sh), samp.size, replace = F),]
  res.hemisphere$cor.sh.sr[i] <- cor.test(temp$lat, temp$sr)$estimate
  res.hemisphere$pval.sh.sr[i] <- cor.test(temp$lat, temp$sr)$p.value
  res.hemisphere$cor.sh.mrd[i] <- cor.test(temp$lat, temp$mrd)$estimate
  res.hemisphere$pval.sh.mrd[i] <- cor.test(temp$lat, temp$mrd)$p.value
  
  if(!i%%1)cat(i,"\r")
}
# barplots mit significances filling
pres <- pivot_longer(res.hemisphere[,c(grep("cor|^n", names(res.hemisphere)))], cols=contains(c("cor")), names_to = "name", values_to = "cor")
pres2 <- pivot_longer(res.hemisphere[,c(grep("pval|^n", names(res.hemisphere)))], cols=contains(c("pval")), names_to = "name", values_to = "pval")
pres$name <- gsub("cor.", "", pres$name)
pres2$name <- gsub("pval.", "", pres2$name)
pres3 <- merge(pres,pres2,all=TRUE)
rm(pres,pres2)

pres3$sig <- pres3$pval<0.05
aggregate(cor~name,data=pres3,mean)
aggregate(pval~name,data=pres3,mean)

ggplot(pres3, aes(cor, fill=sig))+
  geom_histogram()+
  geom_vline(data=ddply(pres3, "name", summarize, cormean = mean(cor)), aes(xintercept=cormean),
             lty=2)+
  geom_vline(xintercept=0)+
  facet_wrap(~name)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 8))+
  scale_fill_startrek(name="p-value<0.05")



