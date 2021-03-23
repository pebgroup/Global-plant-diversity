#setwd("~/Documents/WOLF/PROJECTS/58 World Checklist paper/analyses 2019")

source("functions.R")

readRDS("../data/comm_Nov2020.rds")
load("plant_sr_data/comm.RData")

  bs_spp <- betasim(comm)
  
  save(bs_spp, file = "bs_spp.RData")
  
  save.image("sim_dis.RData")
