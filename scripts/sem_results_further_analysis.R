# get SEM path coefficient estimates table, etc
dat_no.na$level3 <- row.names(dat_no.na)

# Table of estimates ###########################################################################
## add indirect effect names to variables
cat(max.var.mod_mod7)

fm <- "
sr_trans ~ s*soil + trf*sub_trop_mbf + a*area + mrd + ms*mont_gs + t*tra_m + m*mat_m + prs_m + tr0*tri
mrd ~ pre_m + sub_trop_mbf + prs_m + tra_m + s1*soil + tri + mat_m + a3*area 
soil ~ a1*area + ms1*mont_gs + trf1*sub_trop_mbf
sub_trop_mbf ~ p1*pre_m + t1*tra_m + m1*mat_m + a2*area + tr1*tri

# indirect paths on SR
sr_area_via_trf := a1*s + a2*trf
sr_area_total := a + a1*s + a2*trf

sr_subtrop_mbf_via_soil := trf1*trf
sr_subtrop_mbf_total := trf1*trf +trf

sr_pre_m_via_trf := p1*trf

sr_tra_m_via_trf := t1*trf
sr_tra_m_total := t1*trf + t

sr_mat_m_via_trf := m1*trf
sr_mat_m_total :=  m1*trf + m

# indirect paths on MRD
mrd_area_via_soil := a1*s1
mrd_area_total := a1*s1 + a3
"


fm.fit <- sem(fm, data = dat_no.na, estimator="MLM")
summary(fm.fit, standardized = TRUE, rsq=T, fit.measures=TRUE)
fitmeasures(fm.fit, c("cfi.robust", "rmsea.robust", "aic", "pvalue.scaled", "chisq.scaled", "df", "chisq.scaling.factor"))

# res <- parameterestimates(fm.fit, standardized = TRUE, boot.ci.type = "bca.simple")
# reg.coef <- res[res$op=="~",]
# reg.coef.sig <- reg.coef[reg.coef$pvalue<=0.5,]
# reg.coef.ns <- reg.coef[reg.coef$pvalue>0.5,]
# plot(sort(reg.coef.sig$std.all))
# abline(h=c(0,0.2,0.4,-0.2))
# k5 <- kmeans(sort(reg.coef.sig$std.all), centers = 5)
# plot(sort(reg.coef.sig$std.all), col=k5$cluster, pch=20)
# k3 <- kmeans(sort(abs(reg.coef.sig$std.all)), centers = 3)
# plot(sort(abs(reg.coef.sig$std.all)), col=k3$cluster, pch=20)
# abline(h=c(0,0.25,0.5))
# #limits: <-0.2, -0.2-0, 0-0.2, 0.2-0.4, <0.5
# # == 3 strength groups: 0 - 0.2 , 0.25 - 0.5, >0.5 that
# 
# # sorting path coefficients
# reg.coef.sig[order(abs(reg.coef.sig$std.all)),]
# reg.coef.ns[order(abs(reg.coef.ns$std.all)),]
# 
# reg.coef[order(reg.coef$lhs, reg.coef$rhs),] # grouped for regression + alphabetical variables
# 
# 
# 
# pres <- reg.coef[order(reg.coef$lhs, reg.coef$rhs),c("lhs", "op", "rhs", "std.all", "se", "ci.lower", "ci.upper", "pvalue")] # grouped for regression + alphabetical variables
# pres[,c(4:7)] <- round(pres[,c(4:7)],2)
# pres[,c(8)] <- round(pres[,c(8)],3)
# pres[order(pres$lhs, abs(pres$std.all)),]
# 


# documentation notes that standardized = TRUE only displays standardized estimates. "Note that SEs and tests are still based on unstandardized estimates. Use standardizedSolution() to obtain SEs and test statistics for standardized estimates"

pres2 <- standardizedSolution(fm.fit)
# filter for direct and indirect paths
pres2 <- pres2[pres2$op=="~" |pres2$op==":=" ,-6]

# round estimates and p-values separately
pres2[,c(4,5,7,8)] <- round(pres2[,c(c(4,5,7,8))],2)
pres2[,c(6)] <- round(pres2[,c(c(6))],3)

write_xlsx(pres2[order(pres2$op, pres2$lhs, abs(pres2$est.std), decreasing = T),],
           "processed_data/sem_model_coefficients.xls") 
# open that table in google sheets and copy paste it from there to avoid manual work and typos!



# Model path coefficients and correlations ############################################################
lavCor(fm.fit) # observed correlations
inspect(fm.fit, what="cor.all") # model-implied correlations among variables
resid(fm.fit, "cor")
# Large positive values indicate the model under-predicts the correlation; large negative values suggest over-prediction of correlation. Values around |r>.1| are worth closer consideration.
# Note: adding a soil~mat_m path improves model fit, has no further effect though and is not causation



# Univariate VS model regressions ##############################################################
## Scatterplot showing all significant effects for SR and MRD
pres2 <- standardizedSolution(fm.fit)
pres3 <- pres2[pres2$pvalue<0.05 & pres2$op=="~",]
sig.sr <- pres3[pres3$lhs %in% c("sr_trans"),]
plot_sr <- dat_no.na[,c("sr_trans", sig.sr$rhs)]
plot_sr <- pivot_longer(plot_sr, col=-sr_trans)
plot_sr$std.all <- rep(sig.sr$est.std, nrow(dat_no.na))
plot_sr$name[plot_sr$name=="sea_m"] <- "prs_m"
ggplot(plot_sr, aes(x=value, y=sr_trans))+
  geom_point(alpha=0.3)+
  ylab("species richness")+
  facet_wrap(~name, scales = "free")+
  geom_smooth(method="lm", se=F)+
  geom_abline(aes(slope=std.all, intercept=0), col="grey20")+
  theme(strip.background = element_blank())
ggsave("figures/scatterplot_sig_effects_SR.png",  width=7, height=7, units = "in", dpi = 600)

sig.mrd <- pres3[pres3$lhs %in% c("mrd"),]
plot_mrd <- dat_no.na[,c("mrd", sig.mrd$rhs)]
plot_mrd <- pivot_longer(plot_mrd, col=-mrd)
plot_mrd$std.all <- rep(sig.mrd$est.std, nrow(dat_no.na))
plot_mrd$name[plot_mrd$name=="sea_m"] <- "prs_m"
ggplot(plot_mrd, aes(x=value, y=mrd))+
  geom_point(alpha=0.3)+
  ylab("species richness")+
  facet_wrap(~name, scales = "free")+
  geom_smooth(method="lm", se=F)+
  geom_abline(aes(slope=std.all, intercept=0), col="grey20")+
  theme(strip.background = element_blank())
ggsave("figures/scatterplot_sig_effects_MRD.png",  width=5, height=5, units = "in", dpi = 600)


# --> Might be pointing at a Simpsons paradox?

## SR tra_m and tra_m 
cor(dat_no.na[,c("sr_trans", "mat_m", "prs_m", "tra_m")])
n=10
dat_no.na$mat_breaks <- c(1:n)[cut(dat_no.na$mat_m, breaks=n)]
# dat_no.na$tra_breaks <- c(1:n)[cut(dat_no.na$tra_m, breaks=n)]
# dat_no.na$prs_breaks <- c(1:n)[cut(dat_no.na$prs_m, breaks=n)]

simp_tra1<- Simpsons(tra_m, sr_trans, clusterid = mat_breaks, data=dat_no.na, nreps=500) 
simp_sea1<-Simpsons(prs_m, sr_trans, clusterid = mat_breaks, data=dat_no.na, nreps=500) 
# signs for simpsons paradox in both when clustering for mat_m



tra_sr <- ggplot(dat_no.na, aes(x=tra_m, y=sr_trans, col=factor(mat_breaks)))+
    geom_point()+
    scale_color_viridis_d(option = "plasma", 'MAT bins')+
    #scale_color_gradient2('MAT intervals', low = "blue", mid = "yellow", high = "red", midpoint = 5)
    #scale_color_continuous("MAT intervals")+
    guides(colour = guide_legend(reverse=T))+
    geom_smooth(method = "lm", se=F, col="grey")+
    geom_smooth(method = "lm", aes(group=mat_breaks), se=F)+
    theme(legend.position = "none")+
    ylab("species richness")
prs_sr <- ggplot(dat_no.na, aes(x=prs_m, y=sr_trans, col=factor(mat_breaks)))+
    geom_point()+
    scale_color_viridis_d(option = "plasma", 'MAT bins')+
    guides(colour = guide_legend(reverse=T))+
    geom_smooth(method = "lm", se=F, col="grey")+
    geom_smooth(method = "lm", aes(group=mat_breaks), se=F)+
    ylab("species richness")


## MRD area
n=10
dat_no.na$soil_breaks <- c(1:n)[cut(dat_no.na$soil, breaks=n)]
table(dat_no.na$soil_breaks)

simp_area1 <- Simpsons(area, mrd, clusterid = soil_breaks, data=dat_no.na, nreps=500) 

area_mrd <- ggplot(dat_no.na, aes(x=area+abs(min(area)), y=mrd, col=factor(soil_breaks)))+ 
  geom_point()+
  scale_x_continuous("area", trans="sqrt")+
  scale_color_viridis_d(option = "plasma", 'soil bins')+
  geom_smooth(method = "lm", se=F, col="grey")+
  geom_smooth(method = "lm", aes(group=soil_breaks), se=F)+
  guides(colour = guide_legend(reverse=T))#+
empty <- ggplot(dat_no.na, aes(x=area+abs(min(area)), y=mrd, col=factor(soil_breaks)))+ 
  geom_blank()+
  theme_void()
top_row <- plot_grid(ncol = 2, rel_widths = c(0.8,1),
          labels = c("a)", "b)"),
          tra_sr, prs_sr)
bottom_row <- plot_grid(ncol = 2, rel_widths = c(1, 0.77),
                        labels = c("c)", ""),
                        area_mrd, empty)
plot_grid(top_row, bottom_row, nrow=2)
ggsave("figures/simpsons_paradox.png", width=7, height=7, units = "in", dpi = 600)





## Global patterns for tra_m map
### get slopes and p-values for each mat bin
(res <- t(sapply(split(dat_no.na, list(dat_no.na$mat_breaks)),
         function(subd){summary(lm(sr_trans ~ tra_m, data = subd))$coefficients[2,c(1,4)]}))) 
res <- as.data.frame(res)

# get mean mat_m per group
mat_mean <- tapply(dat_no.na$mat_m, dat_no.na$mat_breaks, mean)
res$mat <- as.numeric(mat_mean)
names(res)[2] <- "pvalue"
ggplot(res, aes(x=mat, y=Estimate))+
  geom_point(aes(shape=pvalue<0.05, col=Estimate), size=1.7)+
  #scale_color())+
  scale_y_continuous("SR~tra coefficient")+
  theme(legend.position = "null", 
        plot.background = element_blank(),
    panel.background = element_blank())+
  geom_hline(yintercept = 0, lty=2)
ggsave("figures/slope_mat.png",width=2.5, height=2,units="in")

# check countries with mat_m >0
# library(sf)
# shp <- st_as_sf(shape)
shp$tra_dynamics <- "positive"
shp$tra_dynamics[shp$LEVEL_3_CO %in% dat_no.na$level3[dat_no.na$mat_breaks %in% c(9,10)]] <- "negative"
# add slopes for all bins
res$mat_breaks <- c(1:10)
dat_no.na <- merge(dat_no.na, res, all.x=TRUE)
shp <- merge(shp, dat_no.na[,c("level3","Estimate", "pvalue", "mat_m")], by.x="LEVEL_3_CO", by.y="level3", all.x=TRUE)

shp$Estimate2 <- shp$Estimate
shp$Estimate2[shp$pvalue>0.05] <- NA
(my_plot3 <- ggplot(shp) + 
  geom_sf(aes(fill = Estimate2), lwd=0.1) + 
 scale_fill_continuous("", na.value="grey70")+
  theme_void()+
theme(legend.position = c(0.265, 0.383), # c(0.265, 0.322) excluding ANT
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.key.width = unit(3, "mm"),
      legend.key.height = unit(6, "mm"),
      legend.text = element_blank())
)

logo_file <- readPNG("figures/slope_mat.png")
my_plot_4 <- ggdraw() +
  draw_image(logo_file,  x = -0.35, y = -0.12, scale = .25) + #x = -0.35, y = -0.14, scale = .25
  draw_plot(my_plot3)
my_plot_4
ggsave("figures/tra_scale_map.png", dpi=600, width=10, height=7)


#save.image("processed_data/results.RData")





# Spatial autocorrelation ########################################################################
## Uses code from https://github.com/jebyrnes/spatial_correction_lavaan to get residuals 

load("processed_data/results.RData")
dat_no.na <- dat_no.na[order(dat_no.na$level3),]

# get coordinates 
shp$centroids <- st_centroid(shp) %>% 
  st_coordinates()
shp$y <- shp$centroids[,1]
shp$x <- shp$centroids[,2]

# subset sf object to regions included in SEM
shp <- shp[shp$LEVEL_3_CO %in% dat_no.na$level3,]
shp <- merge(shp, dat_no.na[,c("sr_trans", "level3")], by.x="LEVEL_3_CO", by.y="level3")

# Get model residuals
## Fit the model without indirect effects and meanstructure=TRUE as required for 
## residual function and without path names to for later correction
fm.npn <- "
sr_trans ~ soil + sub_trop_mbf + area + mrd + mont_gs + tra_m + mat_m + prs_m + tri
mrd ~ pre_m + sub_trop_mbf + prs_m + tra_m + soil + tri + mat_m + area 
soil ~ area + mont_gs + sub_trop_mbf
sub_trop_mbf ~ pre_m + tra_m + mat_m + area + tri
"
fm.p <- sem(fm.npn, data = dat_no.na, estimator="MLM", meanstructure=TRUE)
resids <- as.data.frame(my.residuals_lavaan(fm.p)) 
rawdata <- data.frame(inspect(fm.p, "data")) # get raw data from the model
names(resids) <- paste0(names(resids),"_residuals")
plot(rawdata$sr_trans, dat_no.na$sr_trans) # check order for attaching level3 ID
resids$level3 <- dat_no.na$level3

# add residuals to spatial object
shp <- merge(shp, resids, by.x="LEVEL_3_CO", by.y="level3", all.x=TRUE) 

# add fitted SR values to spatial object
fitted_vals <- my.fitted_lavaan(fm.p)
names(fitted_vals) <- paste0(names(fitted_vals),"_fitted")
fitted_vals$level3 <- dat_no.na$level3 
shp <- merge(shp, fitted_vals[,c("sr_trans_fitted", "level3")], by.x="LEVEL_3_CO", by.y="level3", all.x=TRUE)

# plot original vs fitted and residuals
plot_grid(labels = c("a)", "b)", "c)"), ncol = 2,
          ggplot(shp, aes(x=sr_trans, y=sr_trans_fitted))+
            geom_point()+
            xlab("species richness")+
            ylab("species richness fitted")+
            geom_abline(slope=1, intercept=0)
          ,
          ggplot(shp, aes(x=sr_trans, y=sr_trans_residuals))+
            geom_point()+
            xlab("species richness")+
            ylab("species richness SEM residuals")+
            geom_abline(slope=0, intercept=0)
          ,
          ggplot(shp, aes(x=abs(lat), y=sr_trans_residuals, col=sr_trans))+
            geom_point()+
            ylab("species richness SEM residuals")+
            xlab("absolute latitude")+
            scale_color_continuous("SR")+
            geom_smooth(method="lm")
          )
ggsave(file="figures/SEM_residuals_SR.png",
       width=7, height=7, units = "in", dpi = 600)



## Residuals spatial autocorrelation ###################################################

# Weighted distance matrix
distMat <- as.matrix(dist(cbind(shp$x, shp$y)))
distsInv <- 1/distMat # invert matrix for weights
diag(distsInv) <- 0

## GLOBAL autocorrelation
MI.sr <- moran.mc(shp$sr_trans_residuals[,1], mat2listw(distsInv, style = "B"), nsim=599, zero.policy=T)
MI.mrd  <-  moran.mc(shp$mrd_residuals[,1], mat2listw(distsInv, style = "B"), nsim=599, zero.policy=T)
MI.sr$statistic
MI.mrd$statistic

## DISTANCE bands
# distances are in kilometers
moran <- data.frame(dist.class = seq(100, 10000, 100),
                    sr.moransI = NA,
                    sr.moransp = NA,
                    mrd.moransI = NA,
                    mrd.moransp = NA)
coo <- cbind(shp$y, shp$x)
for(i in 1:length(moran$dist.class)){
  S.dist  <-  dnearneigh(coo, 0, moran$dist.class[i], longlat = TRUE)
  lw <- nb2listw(S.dist, style="W",zero.policy=T) 
  
  MI <- moran.mc(shp$sr_trans_residuals[,1], lw, nsim=599,zero.policy=T) 
  moran$sr.moransI[i] <- MI$statistic
  moran$sr.moransp[i] <- MI$p.value
  
  MI.mrd <- moran.mc(shp$mrd_residuals[,1], lw, nsim=599,zero.policy=T) 
  moran$mrd.moransI[i] <- MI.mrd$statistic
  moran$mrd.moransp[i] <- MI.mrd$p.value
  #moran$avg.n[i] <- median(card(lw$neighbours))
  if(!i%%1)cat(i,"\r")
}

temp <- pivot_longer(moran, cols = c("sr.moransI", "mrd.moransI"))
ggplot(temp, aes(x=dist.class, y=value, col=name)) +
  geom_line()+
  scale_x_continuous("Distance class (km)")+
  scale_color_discrete("SEM residuals", labels=c("mrd", "sr"))+
  ylab("Moran's I")
ggsave("figures/SAC_sem_residuals.png", width=5, height=, units = "in", dpi = 600)


## Correct estimates for spatial autocorrelation

# merge coordinates into original model data
coo.df <- data.frame(level3=shp$LEVEL_3_CO, x=shp$x, y=shp$y)
dat_no.na <- merge(dat_no.na, coo.df, all.x=TRUE)

fm.fit.SA <- lavSpatialCorrect_correct(fm.p, xvar=dat_no.na$x, yvar=dat_no.na$y)
# observed autocorrelation matches global calculation

# compare to original
res <- rbind(
  round(fm.fit.SA$parameters$sr_trans[,c(2:6)],3),
  round(fm.fit.SA$parameters$mrd[,c(2:6)],3),
  round(fm.fit.SA$parameters$sub_trop_mbf[,c(2:6)],3),
  round(fm.fit.SA$parameters$soil[,c(2:6)],3)
)
# remove intercepts etc
res <- res[-which(grepl("~~|~1",row.names(res))),]
res$lhs <- gsub("~.*", "", row.names(res))
res$rhs <- gsub(".*~", "", row.names(res))
names(res)[1:5] <- paste0(names(res)[1:5],"_sac")

# merge into original estimate table
all <- merge(pres2, res)
all <- all[,c("lhs", "rhs", "est.std","se", "pvalue",
              "Estimate_sac", "Std.err_sac", "P(>|z|)_sac")]

# write_xlsx(all[order(all$lhs, abs(all$est.std), decreasing = T),],
#            "processed_data/sem_model_coefficients_SAC_correction.xls")



# Maps and  figures ################################################################################

# Global SR map ------------------------------------------------------------------------------------

## get small countries that need a buffer / thicker lines
thicc_lines <-shp[which(shp$area<1200000000),]
lcol <- min(thicc_lines$sr)/max(shp$sr)
ucol <- max(thicc_lines$sr)/max(shp$sr)


(sr_map2 <- ggplot(shp) + 
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
          plot.margin = margin(0, 0, 0, 0, "cm"),
          panel.border = element_blank()
    )+
    coord_sf(xlim = c(-180, 180), ylim = c(-70, 85), expand = F, 
             label_axes = "-NE-")+
    xlab(" ")
)
(sr_ldg_map <- ggplot(shp, aes(lat,sr,col=sr))+
    geom_point(alpha=0.5)+
    scale_color_viridis_c("SR", option = "plasma", trans = "sqrt")+
    scale_y_continuous("SR/1000", trans = "sqrt", 
                       labels = c("0", "5", "10", "15", "20"),
                       breaks = c(0, 5000, 10000, 15000, 20000))+
    scale_x_continuous("", labels = c("","","", "", "","","", ""), #c("60°S","40°S","20°S", "0°", "20°N","40°N","60°N", "80°N")
                       breaks = c(-60,-40, -20, 0, 20,40,60, 80), 
                       position="bottom")+
    geom_smooth(col="grey40", lwd=0.5)+
    #geom_smooth(data=shp[shp$lat<0,], method="lm", col="black")+
    #geom_smooth(data=shp[shp$lat>0,], method="lm", col="black")+
    coord_flip(xlim =c(-70, 85), ylim =c(0, max(shp$sr)+1500),expand = F)+
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.margin = margin(0, 0, 0, -0.50, "cm"),
           legend.position = "none", 
           panel.border = element_blank())
)
plot_grid(sr_map2, sr_ldg_map, 
          ncol = 2, nrow = 1, 
          rel_widths = c(4,1), rel_heights = c(1,1), labels = c("a)", ""))
ggsave("figures/sr_map_AIO.png", dpi=600, width=10, height=3.7)

# Global MRD map ------------------------------------------------------------------------------------

lcol.mrd <- (min(thicc_lines$mrd)-min(shp$mrd))/(max(shp$mrd)-min(shp$mrd))
ucol.mrd <- (max(thicc_lines$mrd)-min(shp$mrd))/(max(shp$mrd)-min(shp$mrd))
mrd_map2 <- ggplot(shp) + 
  geom_sf(aes(fill = mrd), lwd=0.1) + 
  geom_sf(data=thicc_lines, lwd=1.5, aes(col=mrd), show.legend=F)+
  scale_colour_viridis_c(option = "plasma", 
                         begin = lcol.mrd, end = ucol.mrd)+
  scale_fill_viridis_c("MRD", option = "plasma")+ #, 
  theme(legend.position = c(0.21, 0.3),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.border = element_blank()
  )+
  coord_sf(xlim = c(-180, 180), ylim = c(-70, 85), expand = F, 
           label_axes = "-NE-")+
  xlab(" ")
(mrd_ldg_map <- ggplot(shp, aes(lat,mrd,col=mrd))+
    geom_point(alpha=0.5)+
    scale_color_viridis_c("MRD", option = "plasma")+ 
    scale_y_continuous("MRD")+
    scale_x_continuous("", labels = c("","","", "", "","","", ""),
                       breaks = c(-60,-40, -20, 0, 20,40,60, 80), 
                       position="bottom")+
    geom_smooth(col="grey40", lwd=0.5)+
    #geom_smooth(data=shp[shp$lat<0,], method="lm", col="black")+
    #geom_smooth(data=shp[shp$lat>0,], method="lm", col="black")+
    coord_flip(xlim =c(-70, 85),ylim =c(min(shp$mrd)-1, max(shp$mrd)+2),expand = F)+
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.margin = margin(0, 0, 0, -0.50, "cm"),
           legend.position = "none",
           panel.border = element_blank())
)
plot_grid(mrd_map2, mrd_ldg_map, 
          ncol = 2, nrow = 1, 
          rel_widths = c(4,1), labels = c("b)", ""))
ggsave("figures/mrd_map_AIO.png", dpi=600, width=10, height=3.7)



# Latitudinal patterns SR + MRD  -------------------------------------------------------------------

shpbins <- shp
# scale data to values between 0 and 1 for comparison
shpbins$sr_scaled <- range01(shpbins$sr)
shpbins$mrd_scaled <- range01(shpbins$mrd)

shpbins <- st_drop_geometry(shpbins)
temp <- pivot_longer(shpbins[,c("lat","mrd_scaled", "sr_scaled", "CONTINENT")], cols = c("mrd_scaled", "sr_scaled"))
temp <- as.data.frame(temp)

temp$g3 <- rep(NA, nrow(temp))
temp$g3[grep("AMERICA", temp$CONTINENT)] <- "Americas"
temp$g3[grep("ASIA", temp$CONTINENT)] <- "Australasia"
temp$g3[grep("AFRICA|EUROPE", temp$CONTINENT)] <- "Afrope"
temp$g3 <- factor(temp$g3, levels = c("Americas","Afrope","Australasia"))

(continent_scatterplot_comb <- ggplot(na.omit(temp), aes(lat, value, col=name))+
    geom_point(alpha=0.5)+
    scale_y_continuous("scaled MRD and SR", labels =c("0", "0.25", "0.50", "0.75", "1"))+#,
    scale_x_continuous("latitude")+
    facet_wrap(~g3, ncol = 3, labeller=
                 labeller(g3 = c(Americas = "Northern & Southern America",
                                 Afrope = "Africa & Europe",
                                 Australasia="Australasia, Asia-Temperate & -Tropical")))+
    geom_smooth(method="loess", se=T)+
    coord_flip(ylim = c(0,1))+
    scale_colour_startrek(labels=c("MRD", "SR"), name="")+
    #theme_void()+
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line.x = element_line(colour = "grey40"),
          axis.ticks = element_line(colour = "grey40"))+
    guides(color=guide_legend(override.aes=list(fill=NA))) # override geom_smooth() default legend key grey background
)
ggsave("figures/map_lat_patterns.pdf", width=7.6, height=4)
ggsave("figures/map_lat_patterns.png", width=7.6, height=4)

(continent_scatterplot_comb <- ggplot(na.omit(temp), aes(lat, value, col=name))+
    geom_point(alpha=0.5)+
    scale_y_continuous("scaled MRD and SR", labels =c("0", "0.25", "0.50", "0.75", "1"))+#,
    scale_x_continuous("latitude")+
    facet_wrap(~g3, ncol = 3, labeller=
                 labeller(g3 = c(Americas = "Northern & Southern America",
                                 Afrope = "Africa & Europe",
                                 Australasia="Australasia, Asia-Temperate & -Tropical")))+
    geom_smooth(method="loess", se=F)+
    coord_flip(ylim = c(0,1))+
    scale_colour_startrek(labels=c("MRD", "SR"), name="")+
    #theme_void()+
    theme(strip.background = element_blank(),
          strip.text.x = element_text(size = 8),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line.x = element_line(colour = "grey40"),
          axis.ticks = element_line(colour = "grey40"),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA))+
    guides(color=guide_legend(override.aes=list(fill=NA))) # override geom_smooth() default legend key grey background
)
ggsave("figures/map_lat_patterns_alt.pdf", width=7.6, height=4, bg = "transparent")
ggsave("figures/map_lat_patterns_alt.png", width=7.6, height=4, bg = "transparent")


# save region maps
# shp_2 = st_shift_longitude(shp)
# ggplot(shp_2[shp_2$CONTINENT %in% c("NORTHERN AMERICA", "SOUTHERN AMERICA"),]) + 
#   geom_sf(lwd=0, fill = "grey30")+
#   theme(panel.grid = element_blank())
# ggsave("figures/americas.png", width=2, height=1.8, units = "in", dpi = 600)
# 
# ggplot(shp[shp$CONTINENT %in% c("AFRICA", "EUROPE"),]) + 
#   geom_sf(lwd=0, fill = "grey30")+
#   theme(panel.grid = element_blank())
# ggsave("figures/afrope.png", width=1.6, height=2, units = "in", dpi = 600)
# 
ggplot(shp_2[shp_2$CONTINENT %in% c("ASIA-TEMPERATE", "ASIA-TROPICAL", "AUSTRALASIA"),]) +
  geom_sf(lwd=0, fill = "grey30")+
  theme(panel.grid = element_blank())
ggsave("figures/asia.png", width=2, height=2, units = "in", dpi = 600)

ggplot(shp) + 
  geom_sf(lwd=0, fill = "grey50")+
  coord_sf(expand = FALSE)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank())+
  geom_sf(data=shp_2[shp_2$CONTINENT %in% c("ASIA-TEMPERATE", "ASIA-TROPICAL", "AUSTRALASIA"),], lwd=0, fill = "grey90")+
  theme(panel.grid = element_blank())
ggsave("figures/world2.png", width=5, height=3, units = "in", dpi = 600)


## alternative latitude pattern plot
# (global_scatterplot_red <- ggplot(shp, aes(lat, sr))+
#     geom_point(alpha=0.5)+
#     scale_y_continuous("species richness", trans = "sqrt")+
#     scale_x_continuous("latitude centroid botanical country")+
#     geom_smooth(data=shp[shp$lat<0,], method="lm", col="red")+
#     geom_smooth(data=shp[shp$lat>0,], method="lm", col="blue")
# )
# (global_scatterplot_mrd_red <- ggplot(shp, aes(lat, mrd))+
#     geom_point(alpha=0.5)+
#     scale_y_continuous("mean root distance")+
#     scale_x_continuous("latitude centroid botanical country")+    
#     geom_smooth(data=shp[shp$lat<0,], method="lm", col="red")+
#     geom_smooth(data=shp[shp$lat>0,], method="lm", col="blue")
# )
# (global_scatterplot_MRD_SR_red <- ggplot(shp, aes(mrd, sr))+
#     geom_point(alpha=0.5)+
#     scale_y_continuous("species richness", trans = "sqrt")+
#     scale_x_continuous("mean root distance")+
#     geom_smooth(method="lm", col="black")
# )
# plot_grid(nrow=1, global_scatterplot_red, global_scatterplot_mrd_red, global_scatterplot_MRD_SR_red)

## latitudinal stats
cor.test(abs(shp$lat), shp$sr, method = "p") # pearson correlation = -0.38
cor.test(shp$lat[shp$lat<0], shp$sr[shp$lat<0], method = "p")
cor.test(shp$lat[shp$lat>0], shp$sr[shp$lat>0], method = "p")

cor.test(abs(shp$lat), shp$mrd, method = "p")
cor.test(shp$lat[shp$lat<0], shp$mrd[shp$lat<0], method = "p")
cor.test(shp$lat[shp$lat>0], shp$mrd[shp$lat>0], method = "p")
summary(lm(mrd~lat, data=shp[shp$lat<0,]))
summary(lm(mrd~lat, data=shp[shp$lat>0,]))
summary(lm(mrd~abs(lat), data=shp))


# Lat patterns robustness --------------------------------------------------------------------------------
## bootstrap for different sample size differs between hemispheres 
table(shp$lat<0) 
reps <- 500
samp.size <- 50
res.hemisphere <- data.frame(cor.nh.mrd=NA, pval.nh.mrd=NA,
                             cor.sh.mrd=NA, pval.sh.mrd=NA,
                             cor.nh.sr=NA, pval.nh.sr=NA,
                             cor.sh.sr=NA, pval.sh.sr=NA,
                             n=c(1:reps))
set.seed(6348957)
for(i in 1:reps){
  nh <- shp[shp$lat>0,]
  temp <- nh[sample(1:nrow(nh), samp.size, replace = F),]
  res.hemisphere$cor.nh.sr[i] <- cor.test(temp$lat, temp$sr)$estimate
  res.hemisphere$pval.nh.sr[i] <- cor.test(temp$lat, temp$sr)$p.value
  res.hemisphere$cor.nh.mrd[i] <- cor.test(temp$lat, temp$mrd)$estimate
  res.hemisphere$pval.nh.mrd[i] <- cor.test(temp$lat, temp$mrd)$p.value
  
  sh <- shp[shp$lat<0,]
  temp <- sh[sample(1:nrow(sh), samp.size, replace = F),]
  res.hemisphere$cor.sh.sr[i] <- cor.test(temp$lat, temp$sr)$estimate
  res.hemisphere$pval.sh.sr[i] <- cor.test(temp$lat, temp$sr)$p.value
  res.hemisphere$cor.sh.mrd[i] <- cor.test(temp$lat, temp$mrd)$estimate
  res.hemisphere$pval.sh.mrd[i] <- cor.test(temp$lat, temp$mrd)$p.value
  
  if(!i%%1)cat(i,"\r")
}
# barplots mit significances filling
pres <- pivot_longer(res.hemisphere[,c(grep("cor|^n", names(res.hemisphere)))], 
                     cols=contains(c("cor")), names_to = "name", values_to = "cor")
pres2 <- pivot_longer(res.hemisphere[,c(grep("pval|^n", names(res.hemisphere)))], 
                      cols=contains(c("pval")), names_to = "name", values_to = "pval")
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
ggsave("figures/lat_patterns_robustness.png", width=6, height=4.5)







