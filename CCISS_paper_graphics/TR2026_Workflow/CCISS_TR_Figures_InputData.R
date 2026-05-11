

library(terra)
library(ranger)
library(ccissr)
library(data.table)
library(RColorBrewer)
library(bcmaps)
library(sf)
library(climr)
library(here)
library(ggplot2)
library(duckdb)


# -------------------------------------------
# Input data

# dir <- "//objectstore2.nrs.bcgov/ffec/BGC_models/" #permanent location
dir <- "C:/Users/CMAHONY/Data/BGC_models/" #local copy, for speed
BGCmodel <- readRDS(paste0(dir, "BGCmodel_WNA_V4.2gini.rds"))

# dem2 <- rast("//objectstore2.nrs.bcgov/ffec/DEM/dem_BC2kmGrid.tif") #permanent location
dem <- rast("dem_BC2kmGrid.tif")
dem_table <- climr::dem_to_table(dem)

# dir <- "../Common_Files/" # Kiri's location
# dir <- "./tree_summaries_code/" #Mia's location
dir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/" # local copy
bgcs <- vect(paste0(dir, "WNA_BGC_v13_13Apr2026.gpkg"))
bgc_template <- make_bgc_template(dem, bgcs)

# -------------------------------------------
# Setup parameters for model
gcms_cciss <- c("ACCESS-ESM1-5", "CNRM-ESM2-1", "EC-Earth3", "GFDL-ESM4",
                "GISS-E2-1-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")

vars_needed <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
                 "MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", "SHM", "Tave_sm", 
                 "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", "Tmin_at", "Tmin_sm", 
                 "Tmin_sp", "Tmin_wt","CMI_an", "PPT_05","PPT_06","PPT_07","PPT_08",
                 "PPT_09","PPT_at","PPT_wt","CMD_07","CMD_an"
)

# -------------------------------------------
# Setup database

con <- dbCon_cciss("./bc_2km.duckdb", threads = 8)
dbPopulate(con, bgc_template)

# Run Raw BGC Predictions
predict_bgc(con, dem_table, BGCmodel, vars_needed, gcms_cciss, periods_use = list_gcm_periods(), max_runs_use = 3L, obs_2001_2020 = TRUE)

# Run novelty detection
analog_pts <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_simple200_v13_26.csv")
target_pts <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_v4a.csv")
# DBI::dbExecute(con, "DROP TABLE IF EXISTS novelty_raw")
cciss_novelty(con, target_pts, analog_pts, ssps = c("ssp126","ssp245","ssp370", append = TRUE))

## Calculate BGC Persistance/Expansion
bgc_perexp <- bgc_persist_expand(con, by_zone = TRUE)

# compile the species list
data("S1", package = "ccissr")
data("T1", package = "ccissr")
spps <- unique(S1$spp)
spps.candidate <- T1$TreeCode[-which(T1$Exclude=="x" | T1$TreeCode == "X")]
spps <- spps[which(spps%in%spps.candidate)]
spps.native <- spps[which(spps%in%T1[Native=="N", TreeCode])]

## Calculate Species Persistance and Expansion
perexp <- spp_persist_expand(con, spp_list = c("Pl", "Fd", "Cw", "Sx", "At", "Py", "Ba", "Bl", "Bg", "Yc", "Hm", "Lw", "Hw","Dr", "Ep", "Act", "Sb", "Mb", "Ss"), fractional = TRUE)

## Calculate Relative Feasible Area
sa <- spp_suit_area(con, spp_list = spps)

