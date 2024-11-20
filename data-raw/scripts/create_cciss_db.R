library(RPostgres)
library(DBI)
library(data.table)
library(sf)
library(terra)

conn <- DBI::dbConnect(
  drv = RPostgres::Postgres(),
  dbname = Sys.getenv("BCGOV_DB"),
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
                 "Tmin_sp", "Tmin_wt","CMI", "PPT_05","PPT_06","PPT_07","PPT_08","PPT_09","PPT_at","PPT_wt","CMD_07","CMD"
)

load("../Common_Files/BGC_RFresp.Rdata")

all_bgcs <- BGC_RFresp$predictions
bgcs <- data.table(bgc_id = 1:length(levels(all_bgcs)), bgc = levels(all_bgcs))
gcms <- dbGetQuery(conn, "select * from gcm") |> as.data.table()
ssps <- dbGetQuery(conn, "select * from scenario")|> as.data.table()
fps <- dbGetQuery(conn, "select * from futureperiod")|> as.data.table()
fps[,fp_full := gsub("-","_",fp_full)]

# dbExecute(conn, "drop table cciss_future13_array")
# query <- "
#   CREATE TABLE cciss_future13_array (
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


