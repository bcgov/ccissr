library(RPostgres)
library(DBI)
library(data.table)
library(sf)
library(terra)
library(climr)
library(ranger)
library(pool)
library(ccissr)

setwd("/sapho/kdaust/FFEC/")
readRenviron("./ccissr/.Renviron")
conn <- dbPool(
  drv = RPostgres::Postgres(),
  dbname = "cciss",
  host = Sys.getenv("BCGOV_HOST"),
  port = 5432, 
  user = Sys.getenv("BCGOV_USR"),
  password = Sys.getenv("BCGOV_PWD")
)

# qry <- "create table preselected_points13 as (select * from (
# select *, row_number() over (partition by bgc order by random()) as u 
# from bgc_attribution13) as a 
# where u <= 200);"

# qry <- "create table preselected_dist13 as (select * from (
# select *, row_number() over (partition by bgc, dist_code order by random()) as u 
# from bgc_dist_ids) as a 
# where u <= 200);"


# dbExecute(conn, qry)

# hex_grid <- st_read(conn, query = "select * from hex_points")
# bgc <- st_read("~/FFEC/BGC_v13_Fixed.gpkg")

# bgc_att <- st_join(hex_grid,bgc)
# bgc2 <- st_drop_geometry(bgc_att)
# setDT(bgc2)
# bgc <- bgc2[,.(siteno,BGC)]
# setnames(bgc, c("siteno","bgc"))
# bgc <- bgc[!is.na(bgc),]
# dbExecute(conn, "drop table bgc_attribution13")
# dbWriteTable(conn, "bgc_attribution13", bgc, row.names = FALSE)
# dbExecute(conn, "create index on bgc_attribution13(siteno)")
# fwrite(bgc, "bgc_attribution13.csv")
# addVars <- function(dat) {
#   dat[, PPT_MJ := PPT_05 + PPT_06]
#   dat[, PPT_JAS := PPT_07 + PPT_08 + PPT_09]
#   dat[, PPT.dormant := PPT_at + PPT_wt]
#   dat[, CMD.def := pmax(0, 500 - PPT.dormant)]
#   dat[, CMDMax := CMD_07]   ## TODO: THIS IS NOT NECESSARILY CMD MAX
#   dat[, CMD.total := CMD.def + CMD]
# }
# vars_needed <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
#                  "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", "SHM", "Tave_sm", 
#                  "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", "Tmin_at", "Tmin_sm", 
#                  "Tmin_sp", "Tmin_wt","CMI", "PPT_05","PPT_06","PPT_07","PPT_08","PPT_09","PPT_at","PPT_wt","CMD_07","CMD"
# )




################Future model ########################################################
BGC_RFresp <- readRDS("./Common_Files/BGCmodel_WNA_V4.2gini.rds")  
# all_bgcs <- BGC_RFresp$predictions
# bgcs <- data.table(bgc_id = 1:length(levels(all_bgcs)), bgc = levels(all_bgcs))
# # dbWriteTable(conn, "bgc14", bgcs, row.names = FALSE, overwrite = TRUE)
# gcms <- dbGetQuery(conn, "select * from gcm") |> as.data.table()
# ssps <- dbGetQuery(conn, "select * from scenario")|> as.data.table()
# fps <- dbGetQuery(conn, "select * from futureperiod")|> as.data.table()
# fps[,fp_full := gsub("-","_",fp_full)]

# dbExecute(conn, "drop table cciss_future13_array")
# query <- "
#   CREATE TABLE cciss_future14_array (
#     siteno INTEGER REFERENCES hex_points,
#     -- [gcm][scenario][futureperiod][run]
#     bgc_id SMALLINT[13][4][5][3]
#   )
# "
# dbExecute(conn, query)



idx_ipt_len <- function(index, input, length) {
  x <- NA_integer_
  length(x) <- length
  x[index] <- input
  x
}

vars_needed <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
                 "MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", "SHM", "Tave_sm", 
                 "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", "Tmin_at", "Tmin_sm", 
                 "Tmin_sp", "Tmin_wt","CMI_an", "PPT_05","PPT_06","PPT_07","PPT_08","PPT_09","PPT_at","PPT_wt","CMD_07","CMD_an"
)

vars_expert <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
"MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", "SHM", "Tave_sm", 
"Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", "Tmin_at", "Tmin_sm", 
"Tmin_sp", "Tmin_wt", "CMI_an", "PPT_MJ", "PPT_JAS", "CMD.total"
)

hex_pnts <- fread("./Common_Files/Hex_Points_Elev.csv")
pnts <- hex_pnts[!is.na(elev),]

# Create table CCISS Current

res <- downscale(pnts, 
                 which_refmap = "refmap_climr", 
                 obs_periods = "2001_2020",
                 return_refperiod = FALSE,
                 vars = vars_needed)

ccissr::addVars(res)
res <- res[!is.na(Tave_sm),]
res[is.na(res)] <- 0

##predict
message("Predicting...")
temp <- predict(BGC_RFresp, data = res, num.threads = 10)
dat <- cbind(res[,.(id,PERIOD)],temp$predictions)
setnames(dat, old = "V2", new = "BGC")
dat <- dat[,.(id, BGC)]
setnames(dat, c("siteno", "bgc_pred"))
dbWriteTable(conn, "cciss_current14", dat, row.names = F, overwrite = T)
# splits <- c(seq(1,nrow(pnts), by = 10000),nrow(pnts)+1)

# for(i in 280:(length(splits) - 1)){
#   #tic()
#   message("Processing",i)
#   res <- downscale(pnts[splits[i]:(splits[i+1]-1),], 
#                    which_refmap = "refmap_climr", 
#                    gcms = list_gcms(), 
#                    ssps = list_ssps(), 
#                    gcm_periods = list_gcm_periods(), 
#                    max_run = 3L,
#                    ensemble_mean = FALSE,
#                    return_refperiod = FALSE,
#                    vars = vars_needed,
#                    nthread = 8)
  
#   ccissr::addVars(res)
#   res <- res[!is.na(Tave_sm),]
#   res[is.na(res)] <- 0
  
#   ##predict
#   message("Predicting...")
#   temp <- predict(BGC_RFresp, data = res, num.threads = 16)
#   dat <- cbind(res[,.(id,GCM,SSP,RUN,PERIOD)],temp$predictions)
#   setnames(dat, old = "V2", new = "BGC")
  
#   dat[gcms, gcm_id := i.gcm_id, on = c(GCM = "gcm")
#       ][ssps, ssp_id := i.scenario_id, on = c(SSP = "scenario")
#       ][fps, period_id := i.futureperiod_id, on = c(PERIOD = "fp_full")
#       ][bgcs, bgc_id := i.bgc_id, on = c(BGC = "bgc")]
#   dat[,run_id := as.integer(as.factor(RUN)), by = .(GCM,SSP,PERIOD)]
  
#   dat <- dat[!is.na(gcm_id),]
  
#   insert <- dat[,list(bgc_id = paste0(idx_ipt_len(run_id, bgc_id, 3L), collapse = ",")), ##run
#                 by = list(id, gcm_id, ssp_id, period_id)
#   ][,
#     list(bgc_id = paste0("{", idx_ipt_len(period_id, bgc_id, 5L), "}", collapse = ",")), ##period
#     by = list(id, gcm_id, ssp_id)
#   ][,
#     list(bgc_id = paste0("{", idx_ipt_len(ssp_id, bgc_id, 4L), "}", collapse = ",")), ## scenario
#     by = list(id, gcm_id)
#   ][,
#     list(bgc_id = paste0("'{", paste0("{", idx_ipt_len(gcm_id, bgc_id, 13L), "}", collapse = ","), "}'")), ## gcm
#     by = list(id)
#   ][,
#     list(id, bgc_id = gsub("NA", "NULL", bgc_id, fixed = TRUE))
#   ]
#   #toc()
#   query <- paste0("
#     INSERT INTO cciss_future14_array (
#       siteno,
#       bgc_id
#     ) VALUES ",
#                   paste0("(", insert$id, ", ", insert$bgc_id, ")", collapse = ", ")
#   )
  
#   dbExecute(conn, query) 
# }

######### CCISS Novelty ###########################


##############################################################################
query <- "
  CREATE TABLE cciss_novelty14_array (
    siteno INTEGER REFERENCES hex_points,
    -- [gcm][scenario][futureperiod][run]
    novelty SMALLINT[13][4][5][3]
  )
"
#dbExecute(conn, query)
##novelty setup
library(DBI)
library(glue)
# pts <- fread("Common_Files/points_WNA_simple200_v13_26.csv")
# nov_vars <- as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_"))
# clim.pts <- downscale(xyz = pts, which_refmap = "refmap_climr", return_refperiod = TRUE,
#                       vars = nov_vars)
# clim.pts[pts, BGC := i.BGC, on = "id"]

# # historical interannual climatic variability at the geographic centroids of the training points
# pts.mean <- pts[, lapply(.SD, mean), by = BGC]
# pts.mean$id <- 1:dim(pts.mean)[1]
# clim.icv.pts <- downscale(xyz = pts.mean,
#                           which_refmap = "refmap_climr",
#                           obs_years = 1961:1990,
#                           obs_ts_dataset = "cru.gpcc",
#                           return_refperiod = FALSE,
#                           vars = nov_vars)

# hex_pnts <- fread("./Common_Files/Hex_Points_Elev.csv")
# hex_pnts <- hex_pnts[!is.na(elev),]

# # Idea: run novelty calculation by BGC
# con_local <- DBI::dbConnect(duckdb::duckdb(), dbdir = "Common_Files/novelty_temp.duckdb")

# BGC_RFresp <- readRDS("./Common_Files/BGCmodel_WNA_V4.2gini.rds")  
# all_bgcs <- BGC_RFresp$predictions
# bgcs <- data.table(bgc_id = 1:length(levels(all_bgcs)), bgc = levels(all_bgcs))
# gcms <- dbGetQuery(conn, "select * from gcm") |> as.data.table()
# ssps <- dbGetQuery(conn, "select * from scenario")|> as.data.table()
# fps <- dbGetQuery(conn, "select * from futureperiod")|> as.data.table()
# fps[,fp_full := gsub("-","_",fp_full)]
# splits <- c(seq(1,nrow(pnts), by = 5000000),nrow(pnts)+1)

# for(g in 2){ #1:length(gcms$gcm)
#   for(s in 1:length(ssps$scenario)){
#     for(p in 1:length(fps$futureperiod_id)){
#       message(glue("Processing GCM: {gcms$gcm[g]}, SSP: {ssps$scenario[s]}, Period: {fps$fp_full[p]}"))
#       res_ls <- list()
#       for(i in 1:(length(splits) - 1)) {
#         #message(glue("Processing split {i} of {length(splits) - 1}"))
#         res <- downscale(hex_pnts[splits[i]:(splits[i+1]-1),], 
#                          which_refmap = "refmap_climr", 
#                          gcms = gcms$gcm[g], 
#                          ssps = ssps$scenario[s], 
#                          gcm_periods = fps$fp_full[p], 
#                          max_run = 3L,
#                          return_refperiod = FALSE,
#                          ensemble_mean = FALSE,
#                          vars = nov_vars,
#                          nthread = 6)
#         res_ls[[i]] <- res[!is.na(Tmin_sm),]
#       }
#       res <- rbindlist(res_ls)
#       rm(res_ls)

#       runs <- data.table(run = unique(res$RUN))
#       runs[, run_id := as.integer(as.factor(run))]

#       bgc_ls <- list()
#       for (run_id in runs$run_id) {
#         tmp <- dbGetQuery(conn, glue("select siteno, bgc_id[1][2][2][{run_id}] from cciss_future14_array;"))
#         setDT(tmp)
#         tmp[,RUN := runs$run[run_id]]
#         bgc_ls[[run_id]] <- tmp
#       }
#       bgc_preds <- rbindlist(bgc_ls)
#       bgc_preds[bgcs, bgc_pred := i.bgc, on = c(bgc_id = "bgc_id")]
#       res[bgc_preds, bgc_pred := i.bgc_pred, on = c(id = "siteno", RUN = "RUN")]
#       res <- na.omit(res)
#       message("Predicting novelty...")
#       res[,novelty := analog_novelty_core(clim.targets = .SD, 
#                                               clim.analogs = clim.pts, 
#                                               label.targets = bgc_pred, 
#                                               label.analogs = clim.pts$BGC, 
#                                               vars = as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_")),
#                                               clim.icvs = clim.icv.pts,
#                                               label.icvs = pts.mean$BGC[clim.icv.pts$id],
#                                               weight.icv = 0.5,
#                                               threshold = 0.95,
#                                               pcs = NULL), by = .(GCM, SSP, PERIOD, RUN)]

#       res <- res[,.(id, GCM, SSP, PERIOD, RUN, novelty)]
#       dbWriteTable(con_local, "novelty_res", res, row.names = FALSE, append = TRUE)
#     }
#   }
# }

# load novelty from duckdb and write to postgres
# pnts <- fread("./Common_Files/Hex_Points_Elev.csv")
# splits <- c(seq(1,nrow(pnts), by = 40000),nrow(pnts)+1)
# message("Total splits:", length(splits) - 1)

# for(i in 132:(length(splits) - 1)) {
#   message("Processing",i)
#   novelty_res <- dbGetQuery(con_local, glue("select * from novelty_res where id >= {splits[i]} and id < {splits[i+1]}")) |> as.data.table()
#   novelty_res <- unique(novelty_res, by = c("id", "GCM", "SSP", "PERIOD", "RUN"))
  
#   novelty_res[gcms, gcm_id := i.gcm_id, on = c(GCM = "gcm")
#   ][ssps, ssp_id := i.scenario_id, on = c(SSP = "scenario")
#   ][fps, period_id := i.futureperiod_id, on = c(PERIOD = "fp_full")
#   ]
#   novelty_res[,run_id := as.integer(as.factor(RUN)), by = .(GCM,SSP,PERIOD)]
#   novelty_res[,novelty_int := as.integer(novelty*10)]
  
#   novelty_res <- novelty_res[!is.na(gcm_id),]
  
#   insert_nov <- novelty_res[,list(novelty_int = paste0(idx_ipt_len(run_id, novelty_int, 3L), collapse = ",")), ##run
#                 by = list(id, gcm_id, ssp_id, period_id)
#   ][,
#     list(novelty_int = paste0("{", idx_ipt_len(period_id, novelty_int, 5L), "}", collapse = ",")), ##period
#     by = list(id, gcm_id, ssp_id)
#   ][,
#     list(novelty_int = paste0("{", idx_ipt_len(ssp_id, novelty_int, 4L), "}", collapse = ",")), ## scenario
#     by = list(id, gcm_id)
#   ][,
#     list(novelty_int = paste0("'{", paste0("{", idx_ipt_len(gcm_id, novelty_int, 13L), "}", collapse = ","), "}'")), ## gcm
#     by = list(id)
#   ][,
#     list(id, novelty = gsub("NA", "NULL", novelty_int, fixed = TRUE))
#   ]

#   query_novelty <- paste0("
#     INSERT INTO cciss_novelty14_array (
#       siteno,
#       novelty
#     ) VALUES ",
#                       paste0("(", insert_nov$id, ", ", insert_nov$novelty, ")", collapse = ", ")
#   )
#   dbExecute(conn, query_novelty)
  
# }