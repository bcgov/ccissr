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
con <- dbCon_cciss("Common_Files/cciss_spatial_v14.duckdb")
dbExecute(con, "SET memory_limit = '38GB'")
dbExecute(con, "SET threads = 12")
# climr::cache_clear()
# dbExecute(con, "update bgc_summary set bgc_prop = 1 where bgc_prop IS NULL;") ##this is  needed for the historic and 2001_2020_obs periods.
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
bgcs <- vect("Common_Files/BC_v13_clipped.gpkg")


# ##load dem
dem <- rast("Common_Files/BC_DEM_200m.tif")
dem_table <- climr::dem_to_table(dem)
# bgcs <- project(bgcs, crs(dem))
bgc_rast <- make_bgc_template(dem, bgcs)
plot(bgc_rast$bgc_rast)

rbgc <- copy(bgc_rast$bgc_rast)
tid <- bgc_rast$ids
tid[, bgc := stringr::str_extract(bgc, "^[A-Z]+")]
col_use <- WNA_BGCs[,.(Zone, ZoneColour)] |> unique()
tid[col_use, on = c("bgc" = "Zone"), colour := i.ZoneColour]
coltab(rbgc) <- tid[,.(bgc_id,colour)]
rgbbgc <- colorize(rbgc, to = "rgb", alpha = T)
writeRaster(rgbbgc, paste0("bgc_rasters","/bgc_Ensemble_1961_1990_ref_Zone.tif"), overwrite=TRUE)  


#dbPopulate(con, bgc_rast)
# plot_bgc_mapped(bgc_rast, by_zone = FALSE, save_location = "bgc_rasters")
# plot_bgc_mapped(bgc_rast, by_zone = TRUE, save_location = "bgc_rasters")


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
# dbExecute(con, "delete from bgc_raw where period = '1961_1990'")
# predict_bgc(con, dem_table, BGCmodel, vars_needed, gcms_cciss, periods_use = list_gcm_periods(), 
#             obs_2001_2020 = FALSE, refperiod = FALSE, start_tile = 342)
#summarise_preds(con)
gcms_use <- c("ACCESS-ESM1-5","EC-Earth3","GISS-E2-1-G","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0")
runs_use <- c("r1i1p1f1","r4i1p1f1","r2i1p3f1","r2i1p1f1","r1i1p1f1","r1i1p1f1")

# for(p in c("1961_1990", "2001_2020_obs")){
#     plot_bgc(con, p, plot_ensemble = TRUE, raster_template = dem, save_location = "bgc_rasters")
#     plot_bgc(con, p, plot_ensemble = TRUE, by_zone = TRUE, raster_template = dem, save_location = "bgc_rasters")
#     plot_bgc(con, p, plot_ensemble = TRUE, raster_template = dem, save_location = "rasters_raw", raw = TRUE)
#     plot_bgc(con, p, plot_ensemble = TRUE, by_zone = TRUE, raster_template = dem, save_location = "rasters_raw", raw = TRUE)
# }

#plot_bgc(con, "1961_1990", plot_ensemble = FALSE, raster_template = dem, save_location = "bgc_rasters")

## reference period
spp_use <- c("Pl", "Sx","Fd","Cw","Hw","Py", "Bl","At", "Ac", "Ep", "Yc", "Pw", "Ss", "Bg", "Lw", "Mb")
dat <- ref_suit(con, spp_use, edatopes = c("B2","C4","D6"))
rast_bgcs <- bgc_rast$ids

breakpoints.suit <- c(1,2,3,5)
palette.suit <-   c("#006400", "#1E90FF", "#EEC900", "#FFFFFF")
suit_cols <- data.table(value = breakpoints.suit,Colour = palette.suit)

for(species in spp_use){
  for (edatope in c("B2","C4","D6")){
    message("Plotting suitability for ", species, " edatope ", edatope)
    sub <- dat[spp == species & Edatopic == edatope,]
    sub[,suit := as.integer(round(suit))]
    sub <- merge(sub, rast_bgcs, by = "bgc", all = TRUE)
    sub[suit > 3, suit := 5]
    X <- copy(bgc_rast$bgc_rast)
    X <- subst(X, sub$bgc_id, sub$suit)
    coltab(X) <- suit_cols
    final_rgb <- colorize(X, to = "rgb", alpha = TRUE)
    writeRaster(final_rgb, paste0("suitability_rasters","/NewFeas_1961_1990_ref_",edatope,"_",species,".tif"), overwrite = T)
  }
}

# dbGetQuery(con, "SELECT cellnum, COUNT(*) AS occurrences
# FROM bgc_raw
# WHERE period = '1961_1990'
# GROUP BY cellnum, period
# HAVING COUNT(*) > 1;")

# dbGetQuery(con, "select * from bgc_raw where period = '1961_1990' order by cellnum limit 20;")

# dbGetQuery(con, "SELECT
#                     cellnum,
#                     COUNT(*) AS n_rows,
#                     COUNT(DISTINCT bgc_pred) AS n_bgc
#                 FROM bgc_raw
#                 WHERE period = '1961_1990'
#                 GROUP BY cellnum
#                 HAVING COUNT(DISTINCT bgc_pred) > 1;")

# predict_bgc_runs(con, dem_table, BGCmodel, vars_needed, gcms_use, runs_use, ssp_use = "ssp245", periods_use = list_gcm_periods(), 
#                  obs_2001_2020 = FALSE, start_tile = 1)

# for (g in seq_along(gcms_use)){
#   for (p in list_gcm_periods()) {
#     plot_bgc(con, p, plot_ensemble = FALSE, gcms_use[g], ssp = "ssp245", runs_use[g], raster_template = dem, save_location = "bgc_rasters")
#     plot_bgc(con, p, plot_ensemble = FALSE, gcms_use[g], ssp = "ssp245", runs_use[g], raster_template = dem, save_location = "rasters_raw", raw = TRUE)

#   }
# }

# ensemble_predictions(con)
# ensemble_predictions(con, by_zone = TRUE, table_name = "ensemble_preds_zone")
# for(p in list_gcm_periods()) {
#   plot_bgc(con, p, plot_ensemble = TRUE, raster_template = dem, save_location = "bgc_rasters")
#   plot_bgc(con, p, plot_ensemble = TRUE, by_zone = TRUE, raster_template = dem, save_location = "bgc_rasters")
# }
# for(p in list_gcm_periods()) {
#   plot_bgc(con, p, plot_ensemble = TRUE, by_zone = TRUE, raster_template = dem, save_location = "bgc_rasters")
# }

for(p in list_gcm_periods()) {
  plot_bgc(con, p, plot_ensemble = TRUE, raster_template = dem, save_location = "rasters_raw", raw = TRUE)
}
for(p in list_gcm_periods()) {
  plot_bgc(con, p, plot_ensemble = TRUE, by_zone = TRUE, raster_template = dem, save_location = "rasters_raw", raw = TRUE)
}

# for(p in c("2001_2020_obs", "1961_1990")){
#     plot_bgc(con, p, plot_ensemble = TRUE, raster_template = dem, save_location = "bgc_rasters")
#     plot_bgc(con, p, plot_ensemble = TRUE, by_zone = TRUE, raster_template = dem, save_location = "bgc_rasters")
# }
# for(p in c("2001_2020_obs", "1961_1990")){
#     plot_bgc(con, p, plot_ensemble = TRUE, raster_template = dem, save_location = "rasters_raw", raw = TRUE)
#     plot_bgc(con, p, plot_ensemble = TRUE, by_zone = TRUE, raster_template = dem, save_location = "rasters_raw", raw = TRUE)
# }

# plot_bgc(con, "1961_1900", plot_ensemble = TRUE, raster_template = dem, save_location = "bgc_rasters")
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

siteseries_preds(con)

dbExecute(con, "delete from cciss_res where FuturePeriod IN ('2001_2020_obs')")
spp_list <- c("Pl", "Sx","Fd","Cw","Hw","Py", "Bl","At", "Ac", "Ep", "Yc", "Pw", "Ss", "Bg", "Lw", "Mb") # "Pl",
for(spp in spp_list){
  message("Processing ", spp)
  cciss_full_species(con, spp, periods = c("2001_2020_obs"))
  # dbExecute(con, "select distinct FuturePeriod from cciss_res")
}

# t3 <- dbGetQuery(con, "select * from cciss_res where Spp = 'Pl' and FuturePeriod = '1961_1990' and Newsuit < 1 order by Newsuit")
# setDT(t3)

# for (spp in spp_list){
#   for (period in c("1961_1990", "2001_2020_obs")){
#     for (edatope in c("B2","C4","D6")){
#       message("Plotting suitability for ", spp, " ", period, " edatope ", edatope)
#       plot_predicted_suitability(con, spp, period, edatope, raster_template = dem, save_location = "rasters_raw", raw = TRUE)
#       plot_predicted_suitability(con, spp, period, edatope, raster_template = dem, save_location = "suitability_rasters", raw = FALSE)
#       if(period != "1961_1990"){
#         message("Plotting suitability change for ", spp, " ", period, " edatope ", edatope)
#         plot_suitability_change(con, spp, period, edatope, raster_template = dem, save_location = "rasters_raw", raw = TRUE)
#         plot_suitability_change(con, spp, period, edatope, raster_template = dem, save_location = "suitability_rasters", raw = FALSE)
#       }
#     }
#   }
# }

# ##novelty
pts <- fread("Common_Files/points_WNA_simple200_v13_26.csv")
# dbExecute(con, "drop table if exists novelty_raw")
cciss_novelty(con, dem_table, pts, observed = FALSE, table_name = "bgc_raw_runs", append = FALSE)
cciss_novelty(con, dem_table, pts, observed = TRUE, table_name = "bgc_raw", append = TRUE)
ensemble_novelty(con)
for(p in list_gcm_periods()) {
  plot_novelty(con, period = p, raster_template = dem, raw = TRUE, save_location = "rasters_raw")
  plot_novelty(con, period = p, raster_template = dem, raw = FALSE, save_location = "novelty_rasters")

}
for (g in seq_along(gcms_use)){
  for (p in list_gcm_periods()) {
    message("Plotting novelty for ", gcms_use[g], " ", p)
    plot_novelty(con, raster_template = dem, period = p, ensemble = FALSE, gcm = gcms_use[g], ssp = "ssp245", run = runs_use[g],  save_location = "rasters_raw", raw = TRUE)
  }
}

## CCISS Spatial DB

dbGetQuery(con, "select * from cciss_res limit 10;")
dbGetQuery(con, "select count(*) from cciss_res where Curr <= 3 and Prop1 ")
readRenviron("ccissr/.Renviron")
library(pool)
dbCon <- dbPool(
  drv = RPostgres::Postgres(),
  dbname = "cciss_spatial",
  host = Sys.getenv("BCGOV_HOST"),
  port = 5432, 
  user = Sys.getenv("BCGOV_USR"),
  password = Sys.getenv("BCGOV_PWD")
)

spp_use <- c("Pl", "Sx","Fd","Cw","Hw","Py", "Bl","At", "Ac", "Ep", "Yc", "Pw", "Ss", "Bg", "Lw", "Mb")
#eda <- data.table(Edatope = c("B2","C4","D6"), eda_code = 1L:3L)

cw_spp <- data.table(spp = sort(spp_use), spp_id = seq_along(spp_use))
dbWriteTable(con, "cw_spp", cw_spp, row.names = FALSE, overwrite = TRUE)


message("Finished cciss_feas table...")

dbExecute(con, "DROP TABLE IF EXISTS cciss_feas_array")

dbExecute(con, "
CREATE TABLE cciss_feas_array (
  cellid INTEGER,
  spp_id SMALLINT,
  edatope SMALLINT,
  suit SMALLINT[],
  prop1 SMALLINT[],
  prop2 SMALLINT[],
  prop3 SMALLINT[]
)
")

for (spp in spp_use) {
  message('Processing ', spp)

  spp_id <- cw_spp$spp_id[cw_spp$spp == spp]

  dbExecute(con, glue::glue_sql("
    INSERT INTO cciss_feas_array
    WITH base AS (
  SELECT
    CAST(SiteRef AS INTEGER) AS cellid,
    {spp_id}::SMALLINT AS spp_id,
    CASE
      WHEN Edatope = 'B2' THEN 1
      WHEN Edatope = 'C4' THEN 2
      WHEN Edatope = 'D6' THEN 3
    END::SMALLINT AS edatope,

    CASE
      WHEN FuturePeriod = '2001_2020_obs' THEN 1981
      ELSE CAST(substr(FuturePeriod, 1, 4) AS INTEGER)
    END AS fp_code,

    CAST(Curr * 100 AS SMALLINT) AS curr,
    CAST(Newsuit * 100 AS SMALLINT) AS suit,
    CAST(Prop1 * 100 AS SMALLINT) AS prop1,
    CAST(Prop2 * 100 AS SMALLINT) AS prop2,
    CAST(Prop3 * 100 AS SMALLINT) AS prop3

  FROM cciss_res
  WHERE Spp = {spp}),

grouped AS (
  SELECT
    cellid,
    spp_id,
    edatope,

    [
      max(curr),
      max(CASE WHEN fp_code = 1981 THEN suit END),
      max(CASE WHEN fp_code = 2001 THEN suit END),
      max(CASE WHEN fp_code = 2021 THEN suit END),
      max(CASE WHEN fp_code = 2041 THEN suit END),
      max(CASE WHEN fp_code = 2061 THEN suit END),
      max(CASE WHEN fp_code = 2081 THEN suit END)
    ]::SMALLINT[] AS suit,

    [
      max(CASE WHEN fp_code = 1981 THEN prop1 END),
      max(CASE WHEN fp_code = 2001 THEN prop1 END),
      max(CASE WHEN fp_code = 2021 THEN prop1 END),
      max(CASE WHEN fp_code = 2041 THEN prop1 END),
      max(CASE WHEN fp_code = 2061 THEN prop1 END),
      max(CASE WHEN fp_code = 2081 THEN prop1 END)
    ]::SMALLINT[] AS prop1,

    [
      max(CASE WHEN fp_code = 1981 THEN prop2 END),
      max(CASE WHEN fp_code = 2001 THEN prop2 END),
      max(CASE WHEN fp_code = 2021 THEN prop2 END),
      max(CASE WHEN fp_code = 2041 THEN prop2 END),
      max(CASE WHEN fp_code = 2061 THEN prop2 END),
      max(CASE WHEN fp_code = 2081 THEN prop2 END)
    ]::SMALLINT[] AS prop2,

    [
      max(CASE WHEN fp_code = 1981 THEN prop3 END),
      max(CASE WHEN fp_code = 2001 THEN prop3 END),
      max(CASE WHEN fp_code = 2021 THEN prop3 END),
      max(CASE WHEN fp_code = 2041 THEN prop3 END),
      max(CASE WHEN fp_code = 2061 THEN prop3 END),
      max(CASE WHEN fp_code = 2081 THEN prop3 END)
    ]::SMALLINT[] AS prop3,

    sum(coalesce(prop1, 0) + coalesce(prop2, 0) + coalesce(prop3, 0)) AS prop_total

  FROM base
  GROUP BY cellid, spp_id, edatope
)

SELECT
  cellid,
  spp_id,
  edatope,
  suit,
  prop1,
  prop2,
  prop3
FROM grouped
WHERE suit[1] <= 300
   OR prop_total > 10;
  ", .con = con))
}

dbExecute(con, "COPY (
                SELECT
                  cellid,
                  spp_id,
                  edatope,

                  replace(replace(replace(CAST(suit AS VARCHAR),  '[', '{'), ']', '}'), ' ', '') AS suit,
                  replace(replace(replace(CAST(prop1 AS VARCHAR), '[', '{'), ']', '}'), ' ', '') AS prop1,
                  replace(replace(replace(CAST(prop2 AS VARCHAR), '[', '{'), ']', '}'), ' ', '') AS prop2,
                  replace(replace(replace(CAST(prop3 AS VARCHAR), '[', '{'), ']', '}'), ' ', '') AS prop3

                FROM cciss_feas_array
              )
              TO '/sapho/kdaust/FFEC/Common_Files/cciss_feas_array.csv'
              WITH (FORMAT CSV, HEADER TRUE);"
)

### BGC Preds database
dbExecute(dbCon, "drop table bgc_preds")
dbExecute(dbCon, "create table bgc_preds (cellid integer, fp_code smallint, bgc_pred varchar(12), bgc_prop real);")

dat <- dbGetQuery(con, "select * from bgc_summary")
dat[period == "2001_2020_obs", period := "1981"]
dat[,fp_code := as.integer(substr(period, 1, 4))]
dat[,cellnum := as.integer(cellnum)]
setnames(dat, old = "cellnum", new = "cellid")
dat[, period := NULL]
setcolorder(dat, c("cellid", "fp_code", "bgc_pred", "bgc_prop"))

dbWriteTable(dbCon, "bgc_preds", dat, row.names = FALSE, append = TRUE)


  dbExecute(dbCon, "create index on bgc_preds(cellid)")
  message("Finished bgc_preds table!")

## climate data for novelty
dbExecute(dbCon, "drop table if exists future_climate, historic_climate, historic_icv")
pred_vars <- BGCmodel[["forest"]][["independent.variable.names"]] ##required predictors
pts <- fread("Common_Files/points_WNA_simple200_v13_26.csv")

# DEM
dem <- aggregate(dem, fact=8, fun = "max")
# climate data for the biogeoclimatic projections
gcms_use <- c("ACCESS-ESM1-5","EC-Earth3","GISS-E2-1-G","MIROC6","MPI-ESM1-2-HR","MRI-ESM2-0")
runs_use <- c("r1i1p1f1","r4i1p1f1","r2i1p3f1","r2i1p1f1","r1i1p1f1","r1i1p1f1")
ssp_use <- "ssp245"
periods_use <- list_gcm_periods()

grid <- as.data.frame(dem, cells = TRUE, xy = TRUE)
colnames(grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects

res_ls <- list()
for(i in 1:length(gcms_use)){
  clim <- downscale(xyz = grid,
                        gcms = gcms_use[i],
                        ssps = ssp_use,
                        gcm_periods = list_gcm_periods(),
                        run_nm = runs_use[i],
                        vars = list_vars()
  )
  res_ls[[i]] <- clim
}

clim.grid <- rbindlist(res_ls)
addVars(clim.grid)
clim.grid <- clim.grid[is.finite(CMD.total)] #remove NA rows to have complete cases for RF model

#historical climate for training points
#colnames(pts) <- c("id", "BGC", "lon", "lat", "elev") # rename column names to what climr expects
clim.pts <- downscale(xyz = pts,
                    vars = list_vars())
addVars(clim.pts)
clim.pts <- pts[clim.pts, on = "id"]
# Calculate the centroid climate for the training points
clim.pts.mean <- clim.pts[, lapply(.SD, mean), by = BGC, .SDcols = -c("id","PERIOD")]

# historical interannual climatic variability at the geographic centroids of the training points
pts.mean <- pts[, lapply(.SD, mean), by = BGC]
pts.mean$id <- 1:dim(pts.mean)[1]
clim.icv.pts <- downscale(xyz = pts.mean,
                        obs_years = 1951:1990,
                        obs_ts_dataset = "cru.gpcc",
                        return_refperiod = FALSE,
                        vars = list_vars())
addVars(clim.icv.pts)
clim.icv.pts <- pts.mean[clim.icv.pts, on = "id"]

bgc.pred <- predict(BGCmodel, data = clim.grid, num.threads = 4)[['predictions']]
clim.grid[,bgc_pred := bgc.pred]
vars <- c("id", "GCM", "SSP", "RUN", "PERIOD", "bgc_pred", as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")))
clim_final <- clim.grid[,..vars]
clim_final <- clim_final[!is.na(GCM),]

##write to db
dbWriteTable(dbCon, "future_climate", clim_final, row.names = FALSE)
setnames(clim.pts,old = "BGC",new = "bgc")
vars <- c("id", "bgc", as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")))
hist.clim <- clim.pts[,..vars]
dbWriteTable(dbCon, "historic_climate", hist.clim, row.names = FALSE)
dbExecute(dbCon, "create index on historic_climate(bgc)")

clim_icv <- clim.icv.pts
setnames(clim_icv,old = "BGC",new = "bgc")
vars <- c("id", "PERIOD", "bgc", as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")))
clim_icv <- clim_icv[,..vars]
dbWriteTable(dbCon, "historic_icv", clim_icv, row.names = FALSE)
dbExecute(dbCon, "create index on historic_icv(bgc)")

temp <- clim_final[, .(num = .N), by = .(id, PERIOD, bgc_pred)]
ens_win <- temp[temp[, .I[which.max(num)], by = .(id, PERIOD)]$V1]
clim_mean <- clim_final[, lapply(.SD, mean), by = .(id, PERIOD), .SDcols = -c("GCM", "SSP", "RUN", "bgc_pred")]

clim_mean[ens_win, on = c("id", "PERIOD"), bgc_pred := i.bgc_pred]
clim_mean[,GCM := "Ensemble"]
clim_mean[,`:=`(SSP = NA, RUN = NA)]
dbWriteTable(dbCon, "future_climate", clim_mean, row.names = FALSE, append = TRUE)

##load data for SZ ensembles
# clim.grid <- fread(file.path(novelty_folder,"Ensemble_Climates.csv"))
# bgc.pred <- fread(file.path(bgc_folder, "SZ_Ensemble.csv"))
# clim.grid[bgc.pred, bgc_pred := i.bgc_pred, on = c(id = "cellnum", PERIOD = "period")]
# clim.grid[,GCM := "SZ_Ensemble"]
# clim.grid[,`:=`(SSP = NA, RUN = NA)]
# vars <- c("id", "GCM", "SSP", "RUN", "PERIOD", "bgc_pred", as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")))
# clim_final <- clim.grid[,..vars]

# clim_sml <- clim_final[,.SD[sample(.N, min(1000,.N))],by = bgc_pred]


dbWriteTable(dbCon, "future_climate", clim_sml, row.names = FALSE, append = TRUE)
dbExecute(dbCon, "create index on future_climate(\"GCM\", \"PERIOD\", bgc_pred)")
