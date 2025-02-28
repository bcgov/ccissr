library(RPostgres)
library(DBI)
library(data.table)
library(sf)
library(terra)
library(climr)
library(ranger)
library(ccissr)

source("./data-raw/scripts/functions.R")

conn <- DBI::dbConnect(
  drv = RPostgres::Postgres(),
  dbname = "cciss",
  host = Sys.getenv("BCGOV_HOST"),
  port = 5432, 
  user = Sys.getenv("BCGOV_USR"),
  password = Sys.getenv("BCGOV_PWD")
)


addVars <- function(dat) {
  dat[, PPT_MJ := PPT_05 + PPT_06]
  dat[, PPT_JAS := PPT_07 + PPT_08 + PPT_09]
  dat[, PPT.dormant := PPT_at + PPT_wt]
  dat[, CMD.def := pmax(0, 500 - PPT.dormant)]
  dat[, CMDMax := CMD_07]   ## TODO: THIS IS NOT NECESSARILY CMD MAX
  dat[, CMD.total := CMD.def + CMD]
}
vars_needed <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
                 "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", "SHM", "Tave_sm", 
                 "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", "Tmin_at", "Tmin_sm", 
                 "Tmin_sp", "Tmin_wt","CMI", "PPT_05","PPT_06","PPT_07","PPT_08",
                 "PPT_09","PPT_at","PPT_wt","CMD_07","CMD"
)

nov_vars <- as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_"))

vars_needed <- unique(c(vars_needed,nov_vars))

BGCmodel <- readRDS("../Common_Files/BGCmodel_WNA_V2.1.rds")

all_bgcs <- BGCmodel$predictions
bgcs <- data.table(bgc_id = 1:length(levels(all_bgcs)), bgc = levels(all_bgcs))
gcms <- dbGetQuery(conn, "select * from gcm") |> as.data.table()
ssps <- dbGetQuery(conn, "select * from scenario")|> as.data.table()
fps <- dbGetQuery(conn, "select * from futureperiod")|> as.data.table()
fps[,fp_full := gsub("-","_",fp_full)]

# dbExecute(conn, "drop table cciss_future13_array")
query <- "
  CREATE TABLE cciss_future_array (
    siteno INTEGER REFERENCES hex_points,
    -- [gcm][scenario][futureperiod][run]
    bgc_id SMALLINT[13][4][5][3]
  )
"
dbExecute(conn, query)

query <- "
  CREATE TABLE cciss_novelty_array (
    siteno INTEGER REFERENCES hex_points,
    -- [gcm][scenario][futureperiod][run]
    novelty SMALLINT[13][4][5][3]
  )
"
dbExecute(conn, query)

idx_ipt_len <- function(index, input, length) {
  x <- NA_integer_
  length(x) <- length
  x[index] <- input
  x
}


pnts <- fread("../Common_Files/Hex_Points_Elev.csv")
pnts <- pnts[!is.na(elev),]
splits <- c(seq(1,nrow(pnts), by = 40000),nrow(pnts)+1)

##novelty setup
pts <- fread("../Common_Files/points_WNA_simple200.csv")
clim.pts <- downscale(xyz = pts, which_refmap = "refmap_climr",
                      vars = list_vars())
addVars(clim.pts)

# Calculate the centroid climate for the training points
clim.pts.mean <- clim.pts[, lapply(.SD, mean), by = pts$BGC, .SDcols = -c(1,2)]

# historical interannual climatic variability at the geographic centroids of the training points
pts.mean <- pts[, lapply(.SD, mean), by = BGC]
pts.mean$id <- 1:dim(pts.mean)[1]
clim.icv.pts <- downscale(xyz = pts.mean,
                          which_refmap = "refmap_climr",
                          obs_years = 1961:1990,
                          obs_ts_dataset = "cru.gpcc",
                          return_refperiod = FALSE,
                          vars = list_vars())
addVars(clim.icv.pts)


for(i in 1:(length(splits) - 1)){
  #tic()
  message("Processing",i)
  res <- downscale(pnts[splits[i]:(splits[i+1]-1),], 
                   which_refmap = "refmap_climr", 
                   obs_periods = "2001_2020",
                   gcms = list_gcms(), 
                   ssps = list_ssps(), 
                   gcm_periods = list_gcm_periods(), 
                   max_run = 3L,
                   return_refperiod = FALSE,
                   vars = vars_needed)
  
  addVars(res)
  res <- res[!is.na(Tave_sm),]
  res[is.na(res)] <- 0
  res <- res[RUN != "ensembleMean",]
  
  ##predict
  temp <- predict(BGCmodel, data = res, num.threads = 4)
  dat <- cbind(res[,.(id,GCM,SSP,RUN,PERIOD)],temp$predictions)
  dat <- dat[RUN != "ensembleMean",]
  setnames(dat, old = "V2", new = "BGC")
  
  vars_temp <- c("id","GCM","SSP","RUN","PERIOD",nov_vars)
  ##novelty
  clim_nov <- res[,..vars_temp]
  clim_nov[,bgc_pred := dat$BGC]
  
  temp <- clim_nov[GCM == "ACCESS-ESM1-5" & SSP == "ssp245" & RUN == "r1i1p1f1" & PERIOD == "2041_2060",]
  novelty <- analog_novelty(clim.targets = temp, 
                            clim.analogs = clim.pts, 
                            label.targets = temp$bgc_pred, 
                            label.analogs = pts$BGC, 
                            vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
                            clim.icvs <- clim.icv.pts,
                            label.icvs <- pts.mean$BGC[clim.icv.pts$id],
                            weight.icv = 0.5,
                            threshold = 0.95,
                            pcs = NULL
                            
  )
  
  clim_nov[,novelty := analog_novelty(clim.targets = .SD, 
                                      clim.analogs = clim.pts, 
                                      label.targets = bgc_pred, 
                                      label.analogs = pts$BGC, 
                                      vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
                                      clim.icvs <- clim.icv.pts,
                                      label.icvs <- pts.mean$BGC[clim.icv.pts$id],
                                      weight.icv = 0.5,
                                      threshold = 0.95,
                                      pcs = NULL),
           by = .(GCM,SSP,RUN,PERIOD)]
  ####prepare BGC preds
  dat[gcms, gcm_id := i.gcm_id, on = c(GCM = "gcm")
      ][ssps, ssp_id := i.scenario_id, on = c(SSP = "scenario")
      ][fps, period_id := i.futureperiod_id, on = c(PERIOD = "fp_full")
      ][bgcs, bgc_id := i.bgc_id, on = c(BGC = "bgc")]
  dat[,run_id := as.integer(as.factor(RUN)), by = .(GCM,SSP,PERIOD)]
  
  dat <- dat[!is.na(gcm_id),]
  
  insert <- dat[,list(bgc_id = paste0(idx_ipt_len(run_id, bgc_id, 3L), collapse = ",")), ##run
                by = list(id, gcm_id, ssp_id, period_id)
  ][,
    list(bgc_id = paste0("{", idx_ipt_len(period_id, bgc_id, 5L), "}", collapse = ",")), ##period
    by = list(id, gcm_id, ssp_id)
  ][,
    list(bgc_id = paste0("{", idx_ipt_len(ssp_id, bgc_id, 4L), "}", collapse = ",")), ## scenario
    by = list(id, gcm_id)
  ][,
    list(bgc_id = paste0("'{", paste0("{", idx_ipt_len(gcm_id, bgc_id, 13L), "}", collapse = ","), "}'")), ## gcm
    by = list(id)
  ][,
    list(id, bgc_id = gsub("NA", "NULL", bgc_id, fixed = TRUE))
  ]
  
  ##prepare novelty
  clim_nov[gcms, gcm_id := i.gcm_id, on = c(GCM = "gcm")
  ][ssps, ssp_id := i.scenario_id, on = c(SSP = "scenario")
  ][fps, period_id := i.futureperiod_id, on = c(PERIOD = "fp_full")
  ]
  clim_nov[,run_id := as.integer(as.factor(RUN)), by = .(GCM,SSP,PERIOD)]
  clim_nov[,novelty_int := as.integer(novelty*10)]
  
  clim_nov <- clim_nov[!is.na(gcm_id),]
  
  insert_nov <- clim_nov[,list(novelty_int = paste0(idx_ipt_len(run_id, novelty_int, 3L), collapse = ",")), ##run
                by = list(id, gcm_id, ssp_id, period_id)
  ][,
    list(novelty_int = paste0("{", idx_ipt_len(period_id, novelty_int, 5L), "}", collapse = ",")), ##period
    by = list(id, gcm_id, ssp_id)
  ][,
    list(novelty_int = paste0("{", idx_ipt_len(ssp_id, novelty_int, 4L), "}", collapse = ",")), ## scenario
    by = list(id, gcm_id)
  ][,
    list(novelty_int = paste0("'{", paste0("{", idx_ipt_len(gcm_id, novelty_int, 13L), "}", collapse = ","), "}'")), ## gcm
    by = list(id)
  ][,
    list(id, novelty = gsub("NA", "NULL", novelty_int, fixed = TRUE))
  ]
  
  
  #toc()
  query_bgc <- paste0("
    INSERT INTO cciss_future_array (
      siteno,
      bgc_id
    ) VALUES ",
                  paste0("(", insert$id, ", ", insert$bgc_id, ")", collapse = ", ")
  )
  
  query_novelty <- paste0("
    INSERT INTO cciss_novelty_array (
      siteno,
      novelty
    ) VALUES ",
                      paste0("(", insert_nov$id, ", ", insert_nov$novelty, ")", collapse = ", ")
  )
  
  dbExecute(conn, query_bgc) 
  dbExecute(conn, query_novelty)
}

#########create indices#############################
dbSendQuery(conn, glue::glue("
    CREATE INDEX {tb}_idx ON {tb} USING BTREE ({idx});
  ", tb = "cciss_future13_array", idx = "siteno"))
dbSendQuery(conn, glue::glue("VACUUM ANALYZE {tb};", tb = "cciss_future13_array"))

################current period prob#######################
load("../Common_Files/BGC_RFprob.Rdata")

splits <- c(seq(1,nrow(pnts), by = 10000),nrow(pnts)+1)


for(i in 2:(length(splits) - 1)){
  #tic()
  message("Processing",i)
  res <- downscale(pnts[splits[i]:(splits[i+1]-1),], 
                   which_refmap = "refmap_climr", 
                   gcms = list_gcms(), 
                   ssps = list_ssps(), 
                   gcm_periods = list_gcm_periods(), 
                   max_run = 3L,
                   return_refperiod = FALSE,
                   vars = vars_needed)
  
  addVars(res)
  res <- res[!is.na(Tave_sm),]
  res[is.na(res)] <- 0
  
  ##predict
  temp <- predict(BGC_RFresp, data = res, num.threads = 4)
  dat <- cbind(res[,.(id,GCM,SSP,RUN,PERIOD)],temp$predictions)
  dat <- dat[RUN != "ensembleMean",]
  setnames(dat, old = "V2", new = "BGC")
  
  dat[gcms, gcm_id := i.gcm_id, on = c(GCM = "gcm")
  ][ssps, ssp_id := i.scenario_id, on = c(SSP = "scenario")
  ][fps, period_id := i.futureperiod_id, on = c(PERIOD = "fp_full")
  ][bgcs, bgc_id := i.bgc_id, on = c(BGC = "bgc")]
  dat[,run_id := as.integer(as.factor(RUN)), by = .(GCM,SSP,PERIOD)]
  
  dat <- dat[!is.na(gcm_id),]
  
  insert <- dat[,list(bgc_id = paste0(idx_ipt_len(run_id, bgc_id, 3L), collapse = ",")), ##run
                by = list(id, gcm_id, ssp_id, period_id)
  ][,
    list(bgc_id = paste0("{", idx_ipt_len(period_id, bgc_id, 5L), "}", collapse = ",")), ##period
    by = list(id, gcm_id, ssp_id)
  ][,
    list(bgc_id = paste0("{", idx_ipt_len(ssp_id, bgc_id, 4L), "}", collapse = ",")), ## scenario
    by = list(id, gcm_id)
  ][,
    list(bgc_id = paste0("'{", paste0("{", idx_ipt_len(gcm_id, bgc_id, 13L), "}", collapse = ","), "}'")), ## gcm
    by = list(id)
  ][,
    list(id, bgc_id = gsub("NA", "NULL", bgc_id, fixed = TRUE))
  ]
  #toc()
  query <- paste0("
    INSERT INTO cciss_future13_array (
      siteno,
      bgc_id
    ) VALUES ",
                  paste0("(", insert$id, ", ", insert$bgc_id, ")", collapse = ", ")
  )
  
  dbExecute(conn, query) 
}


load("../Common_Files/trainingpts_w_clim_Final.Rdata")

feas <- fread("../Common_Files/Feasibility_v13_2.csv")
dbWriteTable(conn, "feasorig", feas, row.names = FALSE)

#### Simplification ##########
dat <- vect("BGCv13_2_clip_dissolved.gpkg")
d2 <- simplifyGeom(dat, tolerance = 5)
writeVector(d2, "BGC_simplified.gpkg")

bgcs <- st_read("../Common_Files/BGC_simplified.gpkg")
bgcs <- as.data.table(bgcs)
crosswalk <- data.table(BGC = c("CWHxm1","CWHxm2","CWHdm","MHws","MHwsp","CWHms4"), 
                        BGC_new = c("CWHdm1","CWHdm2","CWHdm3","MHms","MHmsp","CWHws3"))

bgcs[crosswalk, BGC_new := i.BGC_new, on = "BGC"]
bgcs[!is.na(BGC_new), BGC := BGC_new]
bgcs[,c("BGC_new","MAP_LABEL") := NULL]
bgc2 <- st_as_sf(bgcs)
st_write(bgc2,"BGC_v13_Fixed.gpkg")

###BGC attribution
bgc <- vect("../Common_Files/BGC_simplified.gpkg")
bgc4326 <- project(bgc, "epsg:4326")
hex <- fread("../Common_Files/Hex_Points_Elev.csv")
hex <- hex[!is.na(elev),]
hexsf <- st_as_sf(hex, coords = c("lon","lat"), crs = 4326)
hexv <- vect(hexsf["id"])
hexv <- project(hexv, "epsg:3005")
bgc_att <- intersect(bgc, hexv)

dat <- fread("../../../Downloads/WNAv13_Subzone_colours_2.csv")
dat2 <- dat[,.(classification = BGC, colour = RGB)]
dat2 <- dat2[classification != "",]
fwrite(dat2, "WNAv13_SubzoneCols.csv")


run <- data.table(run_id = 1:3, run = 1:3)
dbWriteTable(conn, "run", run, row.names = F)

dbWriteTable(conn, "bgcv13", bgcs, row.names = F)

tq <- "SELECT ROW_NUMBER() OVER(ORDER BY gcm_id, scenario_id, futureperiod_id, run_id) row_idx,
               gcm,
               scenario,
               futureperiod,
               run
        FROM gcm 
        CROSS JOIN scenario
        CROSS JOIN futureperiod
        CROSS JOIN run"




# hex_grid <- st_read(conn, query = "select * from hex_points")
# hex_ll <- st_transform(hex_grid, 4326)
# dem <- rast("../Common_Files/WNA_DEM_SRT_30m_cropped.tif")
# crs(dem)
# hex_ll <- vect(hex_ll)
# elev <- extract(dem, hex_ll)
# 
# coords <- as.data.frame(hex_ll,geom = "XY")
# coords$elev <- elev$WNA_DEM_SRT_30m
# colnames(coords) <- c("id","lon","lat","elev")
# fwrite(coords, "../Common_Files/Hex_Points_Elev.csv")
# 
# rm(hex_grid, hex_ll, dem)
# pnts <- as.data.table(coords)
# pnts <- pnts[!is.na(elev),]
# 
# res <- downscale(pnts[1:10000,], 
#                  which_refmap = "refmap_climr", 
#                  gcms = list_gcms(), 
#                  ssps = list_ssps(), 
#                  gcm_periods = list_gcm_periods(), 
#                  max_run = 3L,
#                  vars = list_vars())
# 
# tmp <- c("id", "GCM", "SSP", "RUN", "PERIOD",vars_needed)
# res <- res[,..tmp]
# addVars(res)
# fwrite(res,"Test_Climr_Vars.csv")
# res <- fread("Test_Climr_Vars.csv")
# res <- res[!is.na(Tave_sm),]
# res[is.na(res)] <- 0
# 
# setnames(res, old = "DDsub0_sp", new = "DD_0_sp")
# BGCmodel <- readRDS("../Common_Files/WNA_BGCv12_10May24.rds")
# temp <- predict(BGCmodel, data = res, num.threads = 4)
# dat <- cbind(res[,.(id,GCM,SSP,RUN,PERIOD)],temp$predictions)
# dat <- dat[RUN != "ensembleMean",]
# setnames(dat, old = "V2", new = "BGC")
# 
# fwrite(dat,"TestPredictions.csv")
# dat <- fread("TestPredictions.csv")
# t1 <- res[id == 2043 & GCM == "IPSL-CM6A-LR",]
# 
# gcms <- dbGetQuery(conn, "select * from gcm") |> as.data.table()
# #gcms <- gcms[gcm_id != 9,]
# ssps <- dbGetQuery(conn, "select * from scenario")|> as.data.table()
# bgcs <- dbGetQuery(conn, "select * from bgc")|> as.data.table()
# fps <- dbGetQuery(conn, "select * from futureperiod")|> as.data.table()
# fps[,fp_full := gsub("-","_",fp_full)]
# 
# ##join with ids
# dat[gcms, gcm_id := i.gcm_id, on = c(GCM = "gcm")
#     ][ssps, ssp_id := i.scenario_id, on = c(SSP = "scenario")
#       ][fps, period_id := i.futureperiod_id, on = c(PERIOD = "fp_full")
#         ][bgcs, bgc_id := i.bgc_id, on = c(BGC = "bgc")]
# dat[,run_id := as.integer(as.factor(RUN)), by = .(GCM,SSP,PERIOD)]
# 
# dat <- dat[!is.na(gcm_id),]
# #dat <- dat[gcm_id != 9,]
# 
# setorder(dat, id, PERIOD, GCM,SSP, RUN)






# dt <- dat[id == 9511,]
# insert1 <- dt[,
#               list(bgc_id = paste0(idx_ipt_len(run_id, bgc_id, 3L), collapse = ",")), ##run
#               by = list(id, gcm_id, ssp_id, period_id)
# ]
# t1 <- insert1[,.(Num = .N), by = .(id)]
# t1 <- insert1[,.(Period = .N), by = .(gcm_id,ssp_id)]
# t1 <- insert1[,.(ssp = .N), by = .(gcm_id)]

if (nrow(dt) > 0) {
  insert <- dat[,list(bgc_id = paste0(idx_ipt_len(run_id, bgc_id, 3L), collapse = ",")), ##run
    by = list(id, gcm_id, ssp_id, period_id)
  ][,
    list(bgc_id = paste0("{", idx_ipt_len(period_id, bgc_id, 5L), "}", collapse = ",")), ##period
    by = list(id, gcm_id, ssp_id)
  ][,
    list(bgc_id = paste0("{", idx_ipt_len(ssp_id, bgc_id, 4L), "}", collapse = ",")), ## scenario
    by = list(id, gcm_id)
  ][,
    list(bgc_id = paste0("'{", paste0("{", idx_ipt_len(gcm_id, bgc_id, 13L), "}", collapse = ","), "}'")), ## gcm
    by = list(id)
  ][,
    list(id, bgc_id = gsub("NA", "NULL", bgc_id, fixed = TRUE))
  ]
  
  #temp <- insert[grep("{{{{49,393,308},{393,393,393},{393,393,393},{40,308,40},{308,40,308}},{{49,393,49},{393,393,393},{40,40,40},{40,40,308},{40,40,40}}",bgc_id, fixed = T),]
  query <- paste0("
    INSERT INTO cciss_future13_array (
      siteno,
      bgc_id
    ) VALUES ",
    paste0("(", insert$id, ", ", insert$bgc_id, ")", collapse = ", ")
  )
  
  dbExecute(conn, query) 
}


