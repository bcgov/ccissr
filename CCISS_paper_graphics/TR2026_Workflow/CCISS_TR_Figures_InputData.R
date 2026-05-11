

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
bgc_template <- make_bgc_template(dem2, bgcs)

# -------------------------------------------
# Setup database

con <- dbCon_cciss("./bc_2km.duckdb", threads = 8)
dbPopulate(con, bgc_template)
