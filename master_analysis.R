# March 2022
# 
#
# Most parts written by Melanie Tietje, contributions by Wolf L. Eiserhardt
#  
# 
# This script lists (and sources) all analysis scripts with a brief description.
# 
# Please note that some scripts have a long runtime (hours) and require more
# computational resources than the average computer (e.g. RAM and multicore
# parallel processing). To save time and make results accessible without
# requiring access to a computer cluster, we provide the results of these steps.
# It is best to run those on a computational cluster, as some simply require
# excessive amounts of memory (e.g. TACT: 140GB). Bash scripts included are set
# up to run on a cluster with Slurm Workload Manager. Scripts with a long
# runtime are marked as such. Since internal parallelization in R proves
# difficult on some server setups, parallel processes are performed by cutting
# data in chunks, starting several processes.
# 
# 
# Scripts 1-18 are data compilation, prep and model selection, for working with
# the readily assembled data see 19_sem_results_further_analysis.R
#
#
# Tested with R version 4.1.2., RStudio 2022.02.0+443 "Prairie Trillium" Release
# (9f7969398b90468440a501cf065295d9050bb776, 2022-02-16) for Ubuntu Bionic



# WCVP data prep #########################################################


# *** Check family names following APG IV ---------------------------------------
source("1_APG_family_lookup.R")


# *** Clean WCVP --------------------------------------------------------------
source("2_WCP_cleanup.R")




# PHYLOGENY  ##############################################################

## 1) match tip labels Smith & Brown phylogeny with WCVP names
## 2) add missing species taxonomically
## 3) calculate root distances and DR for all tips
## 4) get verage per TDWG level 3 unit


# *** Match phylogeny tip labels to WCVP, part 1 -----------------------------
source("3_match_data.R")


# *** Prep tip labels for taxonomy matching -----------------------------------
source("4_common_format_creator_gbif.R")



# *** Check family names following APG IV ---------------------------------------
#!! set the database to GBIF
source("1_APG_family_lookup.R")



# *** Run taxonomy matching for GBIF tip labels -------------------------------

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
DB.name <- "GBIF"
data_folder_path <- "processed_data/" 
results_folder_path <- "processed_data/"
wcvp_input_filename <- "apg_wcp_jul_21_clean.rds"
gbif_input_filename <- "apg_input_tip_labels.rds"
fin_file_species <- paste0("fin_species_match_", DB.name, ".rds")

source("5_taxonomic_matcher.v.1.4.R")
# output: fin_species_match_GBIF.rds


# *** Match phyologeny tip labels to WCVP, part 2 -----------------------------
source("6_match_data_part2.R")
# writes allmb_matched_clean.tre



# *** Adding species ----------------------------------------------------------
## LONG RUNTIME ! ca 11 hours

source("7_add_species_server_version.R", print.eval=T)
# writes allmb_matched_added_species_Sep21_clean.tre



# *** Resolve polytomys, get diversification measures -------------------------
## LONG RUNTIME 
## PARALLEL PROCESSING

## Root Distance  (ca 1h on 16 cores)
source("8_resolve_polytomies.R", print.eval = TRUE)


## prep TACT runs
source("8.1.1_TACT_and_diversification_rate.R")

## TACT and DR
# 
# Submit 50 jobs using the job array bash script
# "processed_data/tact/tact_job_array.sh". Each of the started jobs needs 145GB
# memory and has a runtime of 90h. Take care files are all located within the
# same folder. Required files are 'goodsp.csv' and backbone tree
# 'gbmb_matched_no_misplaced.tre', as well as a TACT installation (in this case
# we call a .sif file). If all goes well, you will have 50 TACTed trees. We
# provide all TACTed trees in the supplement.

## process TACTed trees, get DR
source("8.1.2_TACT_and_diversification_rate.R")


# GEOGRAPHY ##################################################################


# *** Build presence matrix -------------------------------------------------
source("9_process_geography.R", print.eval = TRUE)



# *** Mean root distance for each TDGW level 3 unit -------------------------
source("10_phylostruct.R", print.eval = TRUE)



# ENVIRONMENT  ##############################################################
## LONG RUNTIME ! ca 6 hours

# *** Climate and climate anomalies -------------------------------------------
source("11_get_climate.R")


# *** Get number of soil types ------------------------------------------------
source("12_get_soil.R")


# *** Get terrain ruggedness and elevational range ----------------------------
source("13_get_topography.R")


# *** Get biome data ----------------------------------------------------------
source("14_get_biomes.R")



# Assemble final dataset --------------------------------------------------
source("15_assemble_dataset.R")


# Data prep and checks -----------------------------------------------------
## LONG RUNTIME; PARALLEL PROCESSING (ca 2h)
# 
# Parallelized via bash script. Run 'var_selection_master.sh' for 100 GBMs per
# metric, chose metric in the bash script. Results are being called in the
# following R script.

source("16_data_prep_and_checks.R")
# saves the final dataset for SEM: sem_input_data.rds



# STRUCTURAL EQUATION MODELING ###########################################

# *** SEM model selection -----------------------------------------------------
## LONG RUNTIME ! ca 1 hour on 32 cores
## PARALLEL PROCESSING
#
# Parallel processing as achieved externally by dividing the data into 32
# chunks, each being processed in a different job. Jobs are started running bash
# script 'mod_selection_master.sh', which calls either
# 'mod_selection_sequential.sh' or 'mod_selection_sequential_DivRate.sh',
# calling their respective R scripts.

# run '../processed_data/sem/mod_selection_master.sh'



# *** Process SEMs ------------------------------------------------------------ 
## Get stats from model selection runs and manually adjust models for better fit.
## Explores models with interaction effects and island vs mainland groups.

source("18_mod_selection_SEM_analysis.R", print.eval = TRUE)
source("18.1_mod_selection_DicRate_SEM_analysis.R", print.eval = TRUE)


# Maps, Figures, Output --------------------------------------------------
## model parameters, local effects, spatial autocorrelation assessment and
## correction, latitudinal SR and MRD patterns, robustness

source("19_sem_results_further_analysis.R")





# List R packages & Data --------------------------------------------------

## Called R packages
p <- system2(command = "grep", args='-r library ../processed_data', stdout=TRUE)
p2 <- system2(command = "grep", args='-r library ', stdout=TRUE)
p <- c(p,p2)
p <- gsub(".*library\\(", "", p)
p <- gsub("\\)", "", p)
p <- unique(p)
p

## R packages dependencies
pack <- available.packages(fields="Depends")
pack <- pack[row.names(pack) %in% p,]
depends <- as.character(pack[,c("Depends")])
depends <- gsub("R \\(.*?\\)|,", "", depends)
depends <- strsplit(depends, " ")
depends <- unique(unlist(depends))
depends_f <- unique(gsub("[0-9]|\\(.*|\\.|\\)|\\n|-", "", depends))
depends_f

# loaded data files
ld <- system2(command = "egrep", args="-r readRDS\\|read\\.csv\\|fread\\|read\\.tree", stdout=TRUE)
ld <- unique(gsub(".*?\\.\\./processed_data/|.*?\\.\\./data/", "", ld))
ld <- ld[-grepl("function", ld)] 
ld <- ld[-which(grepl("^[0-9]", ld))]
ld <- gsub("\\\"\\)", "", ld)
ld




# FIN #####
