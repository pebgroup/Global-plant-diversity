# APG iV family lookup

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
library(data.table)
library(magrittr)
library(stringr)
library(dplyr)

# clear workspace, keep functions
rm(list = setdiff(ls(), lsf.str())) 

# process for quicker access later
wcp <- fread("../data/WCVP/world_checklist_names_and_distribution_JUL_21/checklist_names.txt",
             quote="", na.strings = "")
saveRDS(wcp, "../processed_data/wcp_jul_21.rds")


# Choose database, possible values: WCP, NCBI, GBIF 
# You can chose more than 1 (e.g. c("NCBI", "GBIF")) 
db <- "GBIF" 
data_folder_path <- "../processed_data/"
WCP_input_filename <- "wcp_jul_21.rds"


# APG edits ###################################################################################
f.apg <- read.csv("../data/apgweb_parsed.csv", stringsAsFactors = F)
#grep("near_|fossil_", f.apg$Clade)

# remove some clade names don't have a valida name ("near_XXX") or a fossil name ("fossil")

f.apg <- f.apg %>% 
  filter(!(str_detect(Clade, "near_") | str_detect(Clade, "fossil"))) %>% 
  dplyr::select(Syn_Fam, Acc_Fam, Clade)

# Interesting note: "Adoxaceae" is treated as "syn" name of "Viburnaceae" in APWeb, which is conflict with APG IV; need to swap

# > f.apg[grep("Adoxaceae", f.apg$Syn_Fam),]
# Syn_Fam     Acc_Fam      Clade
# 25 Adoxaceae Viburnaceae Dipsacales
f.apg[grep("Adoxaceae", f.apg$Syn_Fam),]$Acc_Fam <- "Adoxaceae"
f.apg[grep("Viburnaceae", f.apg$Syn_Fam),]$Acc_Fam <- "Adoxaceae"

# also one more family names "Rhipogonaceae" was not recorded in APWeb
# so Manual added in "Rhipogonaceae" == "Ripogonaceae"  
# manually adding Osmundaceae as its accepted according to http://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:77126775-1
# adding Tiganophytaceae = Brassicales, appears now on WebAPG

f.apg <- rbind.data.frame(f.apg, c("Rhipogonaceae", "Ripogonaceae", "Liliales"),
                          c("Ripogonaceae", "Ripogonaceae", "Liliales"),
                          c("Osmundaceae", "Osmundaceae", ""),
                          c("Tiganophytaceae", "Tiganophytaceae", "Brassicales")) %>%  
  arrange(Syn_Fam) %>% 
  unique()

# remove special character from Isoetaceae 
f.apg$Syn_Fam <- gsub("Isoëtaceae", "Isoetaceae", f.apg$Syn_Fam)
f.apg$Acc_Fam <- gsub("Isoëtaceae", "Isoetaceae", f.apg$Acc_Fam)

# caution serval invalid strings in fern clade names (check before do)
# grep("_", f.apg$Clade)
#f.apg <- f.apg %>% filter(!(str_detect(Clade, "_"))) # uncomment this line if want to remove ferns
names(f.apg) <- c("family", "family.apg", "order")


# NCBI ###################################################################################
if("NCBI" %in% db){
  NCBI <- read.csv(paste0(data_folder_path, ncbi_input_filename), header=T, stringsAsFactors = F)
  

  #Ripogonaceae
  f.ncbi <- NCBI %>% dplyr::select(family) %>% unique()
  fam.apg <-f.apg%>% dplyr::select(family) %>% unique()
  
  setdiff(f.ncbi$family, fam.apg$family)
  # should be empty.   # [1] "Ripogonaceae" has been acounted for in line 33

  NCBI.tmp <- NCBI %>% 
    dplyr::select(-order) %>% 
    arrange(family) %>% 
    unique()
  NCBI.apg <- left_join(NCBI.tmp, f.apg, by="family")
  saveRDS(NCBI.apg, paste0(data_folder_path, ncbi_output_filename))
  write.csv(NCBI.apg, "./results/Spermatophyta_NCBI_APG_checked.csv", row.names = F, quote = F)
  # checklist <- NULL
  # for (i in NCBI.apg$ncbi_id[duplicated(NCBI.apg$ncbi_id)]){
  #   dd <- NCBI.apg[NCBI.apg$ncbi_id==i,] %>% dplyr::select(family, family.apg, order)
  #   checklist <- rbind.data.frame(checklist, dd)
  # }
  # 
  # checklist <- checklist %>% arrange() %>% unique()  
}

# WCP ###################################################################################
if("WCP" %in% db){
  WCP <- readRDS(paste0(data_folder_path, WCP_input_filename))
  
  # check if all families are represented and add if necessary manually
  setdiff(unique(WCP$family), unique(f.apg$family))
  #  "Incertae_sedis"  "Pseudotubulare"  "Gigaspermaceae" 
  
 
  # Pseudotubulare = not a family name, Synonym entry
  # Gigaspermaceae = moss family
  # Incertae_sedis:
  unique(WCP[grep("Incertae_sedis", WCP$family),]$genus)
  # [1] "Cipum"       "Anonymos"    "Theodoricea" "Thuraria"    "Urceola"     "Pouslowia"   "Euphrona"    "Ivonia"
  
  # remove Gigaspermaceae & incertae_sedis. uncomment this line below if want to keep ferns
  ## Gigaspermatacea: 2 cases, both synonyms (415628-az, 411899-az)
  ## Invertae_sedis: 20 taxa, all synonyms or unplaced. Could be matched as synonym when ignoring the family (step3), so keep them
  moss.inc <- c("Gigaspermaceae")
  WCP1 <- WCP %>% 
    filter(!(family %in% moss.inc)) 
  setdiff(unique(WCP1$family), unique(f.apg$family))
  
  # manual fix
  WCP1 <- WCP1[!WCP1$family=="Pseudotubulare",]
  WCP.apg <- left_join(WCP1, f.apg, by="family")
  
  saveRDS(WCP.apg, paste0(data_folder_path, "apg_", WCP_input_filename))
}


# GBIF #######################################################################################
if("GBIF" %in% db){
  gbif <- readRDS("../processed_data/input_tip_labels.rds") 
  
  (sort(fams <- setdiff(unique(gbif$family), unique(f.apg$family))))
  # manual research
  # Cyclostigmataceae: source = PBDB, probably doubtful: remove
  # Hookeriaceae = moss: Elharveya visicularioides
  # Hymenochaetaceae = fungus, remove
  # Sphaerocystidaceae = green algea, should not be in the tree: Planochloris pyrenoidifera
  # Stixaceae = typo of Stixidaceae (included in APG) --> change family name manually here 
  # Ulotrichaceae = green algae Hormidium ... 181911

  
  gbif$family <- as.character(gbif$family)
  gbif$family[which(gbif$family=="Stixaceae")] <- "Stixidaceae"
  todrop <-  c("Sphaerocystidaceae",
               "Hookeriaceae",
               "Cyclostigmataceae", 
               "Ulotrichaceae",
               "Hymenochaetaceae")
  if(length(which(gbif$family %in% todrop))>0){
    gbif <- gbif[-which(gbif$family %in% todrop),]
  }
  
  setdiff(unique(gbif$family), unique(f.apg$family))
  # should just be NAs for missing family information
  
  gbif.apg <- left_join(gbif, f.apg, by="family")
  
  saveRDS(gbif.apg, paste0(data_folder_path, "apg_input_tip_labels.rds"))
}
  
  