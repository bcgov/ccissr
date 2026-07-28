devtools::load_all("/sapho/kdaust/FFEC/ccissr")
#library(ccissr)
library(duckdb)
library(terra)
library(data.table)
library(climr)
setwd("/sapho/kdaust/FFEC")
con <- dbCon_cciss("Common_Files/cciss_spatial_4_2.duckdb")

flp <- vect("flp_2026.geojson")
extent_dt <- rbindlist(lapply(seq_len(nrow(flp)), function(i) {
  e <- ext(flp[i, ])
  
  data.table(
    ORG_UNIT = flp$ORG_UNIT[i],
    xmin = e$xmin,
    xmax = e$xmax,
    ymin = e$ymin,
    ymax = e$ymax
  )
}))

fwrite(extent_dt, "../CCISS_ShinyApp/app/cciss_spatial/flp_bounds.csv")
# edatopes <- c("B2","C4","D6")
# eda_table <- copy(E1) ##Edatopic table
# eda_table <- eda_table[is.na(SpecialCode),]
# eda_table <- eda_table[Edatopic %in% edatopes,]
# eda_table <- unique(eda_table[,.(BGC,SS_NoSpace,Edatopic)])
# dbWriteTable(con, "edatopic", eda_table, row.names = FALSE, overwrite = TRUE)

# suit <- copy(S1) ##Suitability table
# suit <- na.omit(suit, cols = "spp")
# suit <- suit[,.(bgc,spp,ss_nospace,newfeas)]
# dbWriteTable(con, "suitability", suit, row.names = FALSE, overwrite = TRUE)
# bgcs <- vect("Common_Files/BGC_BC_v13_Jan21.gdb/")
# bgcs <- bgcs["MAP_LABEL"]
# names(bgcs) <- "BGC"

# ##load dem
dem <- rast("Common_Files/BC_DEM_200m.tif")
dem_table <- climr::dem_to_table(dem)
# bgcs <- project(bgcs, crs(dem))
# bgc_rast <- make_bgc_template(dem, bgcs)
# dbPopulate(con, bgc_rast)
vars_needed <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
                 "MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", "SHM", "Tave_sm", "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", "Tmin_at", "Tmin_sm", 
                 "Tmin_sp", "Tmin_wt","CMI_an", "PPT_05","PPT_06","PPT_07","PPT_08","PPT_09","PPT_at","PPT_wt","CMD_07","CMD_an"
)

BGCmodel <- readRDS("Common_Files/BGCmodel_WNA_V4.2gini.rds")

gcms_cciss <- c("ACCESS-ESM1-5", "CNRM-ESM2-1", "EC-Earth3", "GFDL-ESM4",
                "GISS-E2-1-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0")

# dbExecute(con, "DROP TABLE IF EXISTS bgc_raw")
# dbExecute(con, "DROP TABLE IF EXISTS clim_raw")
predict_bgc(con, dem_table, BGCmodel, vars_needed, gcms_cciss, periods_use = list_gcm_periods(), 
            obs_2001_2020 = TRUE, refperiod = TRUE, start_tile = 1)
summarise_preds(con)
gcms_use <- c("ACCESS-ESM1-5","EC-Earth3","GISS-E2-1-G","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0")
runs_use <- c("r1i1p1f1","r4i1p1f1","r2i1p3f1","r2i1p1f1","r1i1p1f1","r1i1p1f1")

# predict_bgc_runs(con, dem_table, BGCmodel, vars_needed, gcms_use, runs_use, ssp_use = "ssp245", periods_use = list_gcm_periods(), 
#                  obs_2001_2020 = TRUE, refperiod = TRUE, start_tile = 1)

# for (g in seq_along(gcms_use)){
#   for (p in list_gcm_periods()) {
#     plot_bgc(con, p, plot_ensemble = FALSE, gcms_use[g], ssp = "ssp245", runs_use[g], raster_template = dem, save_location = "bgc_rasters")
#   }
# }

ensemble_predictions(con)
ensemble_predictions(con, by_zone = TRUE, table_name = "ensemble_preds_zone")
for(p in list_gcm_periods()) {
  plot_bgc(con, p, plot_ensemble = TRUE, raster_template = dem, save_location = "bgc_rasters")
}
for(p in list_gcm_periods()) {
  plot_bgc(con, p, plot_ensemble = TRUE, by_zone = TRUE, raster_template = dem, save_location = "bgc_rasters")
}

for(p in c("2001_2020_obs", "1961_1990")){
    plot_bgc(con, p, plot_ensemble = TRUE, raster_template = dem, save_location = "bgc_rasters")
    plot_bgc(con, p, plot_ensemble = TRUE, by_zone = TRUE, raster_template = dem, save_location = "bgc_rasters")
}
plot_bgc(con, "1961_1900", plot_ensemble = TRUE, raster_template = dem, save_location = "bgc_rasters")
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
dbExecute(con, "SET memory_limit = '60GB'")
dbExecute(con, "SET threads = 32")
# #siteseries_preds(con)

# spp_list <- c("Pl","Sx","Fd","Cw","Hw","Py", "Bl","At", "Ac", "Ep", "Yc", "Pw", "Ss", "Bg", "Lw", "Mb") # 
# for(spp in spp_list){
#   message("Processing ", spp)
#   cciss_full_species(con, spp)
# }

# for (spp in spp_list){
#   for (period in c("2001_2020_obs")){
#     for (edatope in c("B2","C4","D6")){
#       plot_predicted_suitability(con, spp, period, edatope, raster_template = dem, save_location = "suitability_rasters")
#       plot_suitability_change(con, spp, period, edatope, raster_template = dem, save_location = "suitability_rasters")
#     }
#   }
# }

##novelty
pts <- fread("Common_Files/points_WNA_simple200_v13_26.csv")
dbExecute(con, "drop table if exists novelty_raw")
cciss_novelty(con, dem_table, pts, observed = TRUE, table_name = "bgc_raw", append = TRUE)
ensemble_novelty(con)
for(p in list_gcm_periods()) {
  plot_novelty(con, period = p, raster_template = dem)
}
for (g in seq_along(gcms_use)){
  for (p in list_gcm_periods()) {
    message("Plotting novelty for ", gcms_use[g], " ", p)
    plot_novelty(con, raster_template = dem, period = p, ensemble = FALSE, gcm = gcms_use[g], ssp = "ssp245", run = runs_use[g],  save_location = "novelty_rasters")
  }
}
