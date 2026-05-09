#devtools::load_all("/sapho/kdaust/FFEC/ccissr")
library(ccissr)
library(duckdb)
library(terra)
library(data.table)
library(climr)
setwd("/sapho/kdaust/FFEC")
con <- dbCon_cciss("Common_Files/cciss_spatial_4_2.duckdb")
# bgcs <- vect("Common_Files/BGC_BC_v13_Jan21.gdb/")
# bgcs <- bgcs["MAP_LABEL"]
# names(bgcs) <- "BGC"

# ##load dem
# dem <- rast("Common_Files/BC_DEM_200m.tif")
# dem_table <- climr::dem_to_table(dem)
# # bgcs <- project(bgcs, crs(dem))
# # bgc_rast <- make_bgc_template(dem, bgcs)
# # plot(bgc_rast$bgc_rast)

# # dbPopulate(con, bgc_rast)
# vars_needed <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
#                  "MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", "SHM", "Tave_sm", "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", "Tmin_at", "Tmin_sm", 
#                  "Tmin_sp", "Tmin_wt","CMI_an", "PPT_05","PPT_06","PPT_07","PPT_08","PPT_09","PPT_at","PPT_wt","CMD_07","CMD_an"
# )

# BGCmodel <- readRDS("Common_Files/BGCmodel_WNA_V4.2gini.rds")

# gcms_cciss <- c("ACCESS-ESM1-5", "CNRM-ESM2-1", "EC-Earth3", "GFDL-ESM4",
#                 "GISS-E2-1-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")

# # dbExecute(con, "DROP TABLE IF EXISTS bgc_raw")
# # dbExecute(con, "DROP TABLE IF EXISTS clim_raw")
# predict_bgc(con, dem_table, BGCmodel, vars_needed, gcms_cciss, periods_use = list_gcm_periods(), 
#             obs_2001_2020 = TRUE, start_tile = 1715, nthread = 4)

# dbGetQuery(con, "SELECT * FROM bgc_raw LIMIT 10")
# dat <- dbGetQuery(con, "SELECT
#                         cellnum, ssp, gcm, run, period,
#                         COUNT(*) AS n
#                         FROM bgc_raw
#                         GROUP BY cellnum, ssp, gcm, run, period
#                         HAVING COUNT(*) > 1;")
# dbExecute(con, "CREATE OR REPLACE TABLE bgc_raw_dedup AS
# SELECT DISTINCT *
# FROM bgc_raw;")
# dbExecute(con, "DROP TABLE bgc_raw;")
# dbExecute(con, "ALTER TABLE bgc_raw_dedup RENAME TO bgc_raw;")

# eda_table <- copy(ccissr::E1) ##Edatopic table
# eda_table <- eda_table[is.na(SpecialCode),]
# eda_table <- eda_table[Edatopic %in% c("B2","C4","D6"),]
# eda_table <- unique(eda_table[,.(BGC,SS_NoSpace,Edatopic)])
# dbWriteTable(con, "edatopic", eda_table, row.names = FALSE, overwrite = TRUE)

# suit <- copy(ccissr::S1) ##Suitability table
# suit <- na.omit(suit, cols = "spp")
# suit <- suit[,.(bgc,spp,ss_nospace,newfeas)]
# dbWriteTable(con, "suitability", suit, row.names = FALSE, overwrite = TRUE)

# summarise_preds(con)
# siteseries_preds(con)
dbExecute(con, "SET memory_limit = '60GB'")
dbExecute(con, "SET threads = 32")

spp_list <- c("Pl","Sx","Fd","Cw","Hw","Py", "Bl","At", "Ac", "Ep", "Yc", "Pw", "Ss", "Bg", "Lw", "Mb")
for(spp in spp_list){
  message("Processing ", spp)
  cciss_full_species(con, spp)
}