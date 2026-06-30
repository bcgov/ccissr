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


qry <- "create table preselected_points14 as (select * from (
select *, row_number() over (partition by bgc order by random()) as u
from bgc_attribution14) as a
where u <= 200);"

qry <- "create table preselected_dist14 as (select * from (
select *, row_number() over (partition by bgc, dist_code order by random()) as u
from bgc_dist_ids) as a
where u <= 200);"


dbExecute(conn, qry)

hex_grid <- st_read(conn, query = "select * from hex_points")
bgc <- st_read("../Common_Files/WNA_BGCv13/BGC_BC_v13_Jan21.gdb")
bgc <- bgc["MAP_LABEL"]
bgc_att <- st_join(hex_grid,bgc)
bgc2 <- st_drop_geometry(bgc_att)
setDT(bgc2)
setnames(bgc2, c("siteno","bgc"))
bgc2 <- bgc2[!is.na(bgc),]
dbExecute(conn, "drop table bgc_attribution13")
dbWriteTable(conn, "bgc_attribution14", bgc2, row.names = FALSE)
dbExecute(conn, "create index on bgc_attribution14(siteno)")

###########Current period prob model #######################
# pnts <- fread("~/FFEC/Common_Files/Hex_Points_Elev.csv")
# pnts <- pnts[!is.na(elev),]

# splits <- c(seq(1,nrow(pnts), by = 10000),nrow(pnts)+1)

# # qry <- "CREATE TABLE cciss_prob13 (
# #  siteno INTEGER REFERENCES hex_points,
# #  bgc_pred VARCHAR(12),
# #  prob FLOAT8)"

# # dbExecute(conn, qry)

# for(i in 2:(length(splits) - 1)){
#   #tic()
#   message("Processing",i)
#   res <- downscale(pnts[splits[i]:(splits[i+1]-1),], 
#                    which_refmap = "refmap_climr", 
#                    obs_periods = "2001_2020",
#                    return_refperiod = FALSE,
#                    vars = vars_needed)
  
#   addVars(res)
#   res <- res[!is.na(Tave_sm),]
#   res[is.na(res)] <- 0
  
#   ##predict
#   temp <- predict(BGCmodel_prob, data = res, num.threads = 32)
#   #t2 <- as.data.table(temp$predictions)
  
#   dat <- cbind(res[,.(id,PERIOD)],temp$predictions)
#   dat2 <- melt(dat, id.vars = c("id","PERIOD"))
#   setnames(dat2, c("id","PERIOD","bgc_pred","bgc_prop"))
#   dat2 <- dat2[bgc_prop > 0.005,]
#   dbWriteTable(conn, "cciss_prob13", dat2[,.(siteno = id, bgc_pred, prob = bgc_prop)], row.names = F, append = T)
# }
# message("done!")
################Future model ########################################################
# BGC_RFresp <- readRDS("./Common_Files/BGCmodel_WNA_V4.2gini.rds")  
# all_bgcs <- BGC_RFresp$predictions
# bgcs <- data.table(bgc_id = 1:length(levels(all_bgcs)), bgc = levels(all_bgcs))
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

# vars_needed <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
#                  "MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", "SHM", "Tave_sm", 
#                  "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", "Tmin_at", "Tmin_sm", 
#                  "Tmin_sp", "Tmin_wt","CMI_an", "PPT_05","PPT_06","PPT_07","PPT_08","PPT_09","PPT_at","PPT_wt","CMD_07","CMD_an"
# )

# vars_expert <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
# "MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", "SHM", "Tave_sm", 
# "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", "Tmin_at", "Tmin_sm", 
# "Tmin_sp", "Tmin_wt", "CMI_an", "PPT_MJ", "PPT_JAS", "CMD.total"
# )

# hex_pnts <- fread("./Common_Files/Hex_Points_Elev.csv")
# pnts <- pnts[!is.na(elev),]

# splits <- c(seq(1,nrow(pnts), by = 10000),nrow(pnts)+1)

# for(i in 411:(length(splits) - 1)){
#   #tic()
#   message("Processing",i)
#   res <- downscale(pnts[splits[i]:(splits[i+1]-1),], 
#                    which_refmap = "refmap_climr", 
#                    gcms = list_gcms(), 
#                    ssps = list_ssps(), 
#                    gcm_periods = list_gcm_periods(), 
#                    max_run = 3L,
#                    return_refperiod = FALSE,
#                    vars = vars_needed,
#                    nthread = 8)
  
#   ccissr::addVars(res)
#   res <- res[!is.na(Tave_sm),]
#   res[is.na(res)] <- 0
  
#   ##predict
#   message("Predicting...")
#   res <- res[RUN != "ensembleMean",]
#   temp <- predict(BGC_RFresp, data = res, num.threads = 12)
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

analog_novelty_fast <- function(clim.targets, clim.analogs, label.targets, label.analogs, vars,
                                clim.icvs = NULL, label.icvs = NULL, weight.icv = 0.5,
                                sigma = TRUE, threshold = 0.95,
                                pcs = NULL) {
  
  stopifnot(weight.icv >= 0, weight.icv <= 1)
  
  use_icv <- !is.null(clim.icvs)
  if (use_icv && is.null(label.icvs)) {
    stop("If clim.icvs is supplied, label.icvs must also be supplied.")
  }
  
  analogs <- unique(label.targets)
  novelty <- rep(NA_real_, length(label.targets))
  
  # Precompute indices once
  target_idx <- split(seq_along(label.targets), as.character(label.targets))
  analog_idx <- split(seq_along(label.analogs), as.character(label.analogs))
  if (use_icv) {
    icv_idx <- split(seq_along(label.icvs), as.character(label.icvs))
  }
  
  weight.analog <- 1 - weight.icv
  
  ## log vars
  
  for (analog in analogs) {
    analog_key <- as.character(analog)
    
    idx_a <- analog_idx[[analog_key]]
    idx_t <- target_idx[[analog_key]]
    
    if (is.null(idx_a) || is.null(idx_t)) next
    if (length(idx_a) < 2L || length(idx_t) < 1L) next
    
    # Pull data as matrices
    tmpA <- logVars(clim.analogs[idx_a, ..vars], zero_adjust = TRUE)
    tmpT <- logVars(clim.targets[idx_t, ..vars], zero_adjust = TRUE)

    A0 <- as.matrix(tmpA)
    T0 <- as.matrix(tmpT)
    
    # Clean analogs
    A0 <- A0[complete.cases(A0), , drop = FALSE]
    if (nrow(A0) < 2L) next
    
    # Remove zero-variance / invalid columns based on analog climate
    mu_A <- colMeans(A0, na.rm = TRUE)
    sd_A <- apply(A0, 2, sd, na.rm = TRUE)
    
    keep <- is.finite(sd_A) & sd_A > 0
    if (sum(keep) < 3L) next
    
    A0 <- A0[, keep, drop = FALSE]
    T0 <- T0[, keep, drop = FALSE]
    mu_A <- mu_A[keep]
    sd_A <- sd_A[keep]
    
    # Scale analog and target to analog mean/sd
    A <- sweep(sweep(A0, 2, mu_A, "-"), 2, sd_A, "/")
    T1 <- sweep(sweep(T0, 2, mu_A, "-"), 2, sd_A, "/")
    
    # ICV data: scale using ICV mean but analog sd, matching original logic
    if (use_icv) {
      idx_i <- icv_idx[[analog_key]]
      if (is.null(idx_i) || length(idx_i) < 2L) next
      
      tmpICV <- logVars(clim.icvs[idx_i, ..vars], zero_adjust = TRUE)
      I0 <- as.matrix(tmpICV)
      I0 <- I0[complete.cases(I0), keep, drop = FALSE]
      if (nrow(I0) < 2L) next
      
      mu_I <- colMeans(I0, na.rm = TRUE)
      I <- sweep(sweep(I0, 2, mu_I, "-"), 2, sd_A, "/")
    }
    
    # prcomp cannot handle NA rows; use complete target rows only for PCA sampling
    T_complete <- T1[complete.cases(T1), , drop = FALSE]
    if (nrow(T_complete) < 1L) next
    
    s <- sample.int(
      n = nrow(T_complete),
      size = nrow(A),
      replace = nrow(T_complete) < nrow(A)
    )
    
    T_sample <- T_complete[s, , drop = FALSE]
    
    # PCA on pooled analog + sampled target
    pca <- prcomp(
      rbind(A, T_sample),
      scale = FALSE
    )
    
    pcs_use <- pcs
    if (is.null(pcs_use)) {
      cumvar <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
      pcs_use <- which(cumvar >= threshold)[1]
      pcs_use <- max(3L, pcs_use)
    }
    
    pcs_use <- min(pcs_use, ncol(pca$rotation))
    if (pcs_use < 1L) next
    
    rot <- pca$rotation[, seq_len(pcs_use), drop = FALSE]
    
    # PC scores via matrix multiplication rather than predict()
    PA <- A %*% rot
    PT <- T1 %*% rot
    if (use_icv) PI <- I %*% rot
    
    # Standardize PCs to analog mean and pooled analog/ICV sd
    pc_mu_A <- colMeans(PA, na.rm = TRUE)
    pc_sd_A <- apply(PA, 2, sd, na.rm = TRUE)
    
    if (use_icv) {
      pc_sd_I <- apply(PI, 2, sd, na.rm = TRUE)
      pc_sd_use <- weight.analog * pc_sd_A + weight.icv * pc_sd_I
    } else {
      pc_sd_use <- pc_sd_A
    }
    
    # valid_pc <- is.finite(pc_sd_use) & pc_sd_use > 0
    # if (sum(valid_pc) < 1L) next
    # 
    # PA <- PA[, valid_pc, drop = FALSE]
    # PT <- PT[, valid_pc, drop = FALSE]
    # pc_mu_A <- pc_mu_A[valid_pc]
    # pc_sd_use <- pc_sd_use[valid_pc]
    
    PA <- sweep(sweep(PA, 2, pc_mu_A, "-"), 2, pc_sd_use, "/")
    PT <- sweep(sweep(PT, 2, pc_mu_A, "-"), 2, pc_sd_use, "/")
    
    if (use_icv) {
      pc_mu_I <- colMeans(PI, na.rm = TRUE)
      PI <- sweep(sweep(PI, 2, pc_mu_I, "-"), 2, pc_sd_use, "/")
    }
    
    pcs_final <- ncol(PA)
    
    # Combine spatial covariance and ICV covariance
    cov_A <- var(PA)
    
    if (use_icv) {
      cov_I <- var(PI)
      cov_use <- weight.analog * cov_A + weight.icv * cov_I
    } else {
      cov_use <- cov_A
    }
    
    # Mahalanobis distance
    md2 <- tryCatch(
      mahalanobis(PT, center = rep(0, pcs_final), cov = cov_use),
      error = function(e) rep(NA_real_, nrow(PT))
    )
    
    md <- sqrt(md2)
    
    if (sigma) {
      p <- pchisq(md2, df = pcs_final)
      q <- sqrt(qchisq(p, df = 1))
      q[!is.finite(q)] <- 8 # set infinite values to 8 sigma (outside the decimal precision of pchi) 
      q[is.na(p)] <- NA # reset NA values as NA
      novelty[idx_t] <- q
    } else {
      novelty[idx_t] <- md
    }
  }
  
  novelty
}
##############################################################################
# query <- "
#   CREATE TABLE cciss_novelty14_array (
#     siteno INTEGER REFERENCES hex_points,
#     -- [gcm][scenario][futureperiod][run]
#     novelty SMALLINT[13][4][5][3]
#   )
# "
# dbExecute(conn, query)
##novelty setup
library(EnvStats)
library(DBI)
library(glue)
# pts <- fread("Common_Files/points_WNA_simple200_v13_26.csv")
# clim.pts <- downscale(xyz = pts, which_refmap = "refmap_climr",
#                       vars = list_vars())
# addVars(clim.pts)

# # Calculate the centroid climate for the training points
# clim.pts.mean <- clim.pts[, lapply(.SD, mean), by = pts$BGC, .SDcols = -c(1,2)]

# # historical interannual climatic variability at the geographic centroids of the training points
# pts.mean <- pts[, lapply(.SD, mean), by = BGC]
# pts.mean$id <- 1:dim(pts.mean)[1]
# clim.icv.pts <- downscale(xyz = pts.mean,
#                           which_refmap = "refmap_climr",
#                           obs_years = 1961:1990,
#                           obs_ts_dataset = "cru.gpcc",
#                           return_refperiod = FALSE,
#                           vars = list_vars())
# addVars(clim.icv.pts)
# nov_vars <- as.vector(outer(c("Tmin", "Tmax", "PPT"), c("wt", "sp", "sm", "at"), paste, sep = "_"))

# hex_pnts <- fread("./Common_Files/Hex_Points_Elev.csv")
# hex_pnts <- hex_pnts[!is.na(elev),]

## Idea: run novelty calculation by BGC
con_local <- DBI::dbConnect(duckdb::duckdb(), dbdir = "Common_Files/novelty_temp.duckdb")

BGC_RFresp <- readRDS("./Common_Files/BGCmodel_WNA_V4.2gini.rds")  
all_bgcs <- BGC_RFresp$predictions
bgcs <- data.table(bgc_id = 1:length(levels(all_bgcs)), bgc = levels(all_bgcs))
gcms <- dbGetQuery(conn, "select * from gcm") |> as.data.table()
ssps <- dbGetQuery(conn, "select * from scenario")|> as.data.table()
fps <- dbGetQuery(conn, "select * from futureperiod")|> as.data.table()
fps[,fp_full := gsub("-","_",fp_full)]

# temp <- dbGetQuery(con_local, "select distinct GCM, SSP, PERIOD from novelty_res where GCM == 'CanESM5'") |> as.data.table()
# setorder(temp, GCM, SSP, PERIOD)
# dbExecute(con_local, "delete from novelty_res where GCM == 'CanESM5'")

# for(g in 3){
#   for(s in 1:length(ssps$scenario)){
#     for(p in 1:length(fps$futureperiod_id)){
#       message(glue("Processing GCM: {gcms$gcm[g]}, SSP: {ssps$scenario[s]}, Period: {fps$fp_full[p]}"))
#       res <- downscale(hex_pnts, 
#                        which_refmap = "refmap_climr", 
#                        gcms = gcms$gcm[g], 
#                        ssps = ssps$scenario[s], 
#                        gcm_periods = fps$fp_full[p], 
#                        max_run = 3L,
#                        return_refperiod = FALSE,
#                        ensemble_mean = FALSE,
#                        vars = nov_vars,
#                        nthread = 8)
#       res <- res[!is.na(Tmin_sm),]

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

#       res[,novelty := analog_novelty_fast(clim.targets = .SD, 
#                                               clim.analogs = clim.pts, 
#                                               label.targets = bgc_pred, 
#                                               label.analogs = pts$BGC, 
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

## load novelty from duckdb and write to postgres
pnts <- fread("./Common_Files/Hex_Points_Elev.csv")
splits <- c(seq(1,nrow(pnts), by = 40000),nrow(pnts)+1)
message("Total splits:", length(splits) - 1)

for(i in 1:(length(splits) - 1)) {
  message("Processing",i)
  novelty_res <- dbGetQuery(con_local, glue("select * from novelty_res where id >= {splits[i]} and id < {splits[i+1]}")) |> as.data.table()
  
  novelty_res[gcms, gcm_id := i.gcm_id, on = c(GCM = "gcm")
  ][ssps, ssp_id := i.scenario_id, on = c(SSP = "scenario")
  ][fps, period_id := i.futureperiod_id, on = c(PERIOD = "fp_full")
  ]
  novelty_res[,run_id := as.integer(as.factor(RUN)), by = .(GCM,SSP,PERIOD)]
  novelty_res[,novelty_int := as.integer(novelty*10)]
  
  novelty_res <- novelty_res[!is.na(gcm_id),]
  
  insert_nov <- novelty_res[,list(novelty_int = paste0(idx_ipt_len(run_id, novelty_int, 3L), collapse = ",")), ##run
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

  query_novelty <- paste0("
    INSERT INTO cciss_novelty14_array (
      siteno,
      novelty
    ) VALUES ",
                      paste0("(", insert_nov$id, ", ", insert_nov$novelty, ")", collapse = ", ")
  )
  dbExecute(conn, query_novelty)
  
}

## testing
sites <- 6305115:6305125
gcm_weight <- data.table(gcm = c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CNRM-ESM2-1", "EC-Earth3", 
                                 "GFDL-ESM4", "GISS-E2-1-G", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC6", 
                                 "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"),
                         weight = c(1,0,0,1,1,1,1,0,0,1,1,1,0))

rcp_weight <- data.table(rcp = c("ssp126","ssp245","ssp370","ssp585"), 
                         weight = c(0.8,1,0.8,0))

all_weight <- as.data.table(expand.grid(gcm = gcm_weight$gcm,rcp = rcp_weight$rcp))
all_weight[gcm_weight,wgcm := i.weight, on = "gcm"]
all_weight[rcp_weight,wrcp := i.weight, on = "rcp"]
all_weight[,weight := wgcm*wrcp]


all_weight[, comb := glue::glue_sql(
  "({gcm}, {rcp}, {weight})",
  .con = conn
)]

test <- dbGetCCISS_v13(conn, sites, avg = FALSE, modWeights = all_weight,  
                       cciss_table = "cciss_future14_array", 
                       cciss_observed = "cciss_current14",
                       bgc_table = "bgc_attribution14", 
                       bgc_lookup = "bgc14")
test2 <- dbGetCCISS_novelty(conn, sites, avg = FALSE, modWeights = all_weight,  
                       cciss_table = "cciss_future14_array", 
                       novelty_table = "cciss_novelty14_array",
                       cciss_observed = "cciss_current14",
                       bgc_table = "bgc_attribution14", 
                       bgc_lookup = "bgc14", nov_cutoff = 5)
