
#' Create raster and id table of BGCs based on vector data
#' @param xyz SpatRaster or data.table of points (must have lat, lon, id).
#' @param bgcs SpatVector or sf object of bgc bondaries with column "BGC" specifying name.
#' @return List containing resultant SpatRaster and data.table of ids (or single data.table if input is a data.table).
#' @import data.table
#' @importFrom terra project vect rasterize
#' @importFrom sf st_transform
#' @export
make_bgc_template <- function(xyz, bgcs, res = 0.0008){
  if(inherits(bgcs,"SpatVector")){
    #bgcs <- project(bgcs, "epsg:4326")
    bgcs$bgc_id <- as.numeric(as.factor(bgcs$BGC))
  } else {
    #bgcs <- st_transform(bgcs, 4326)
    bgcs$bgc_id <- as.numeric(as.factor(bgcs$BGC))
    bgcs <- vect(bgcs)
  }
  if(inherits(xyz,"SpatRaster")){
    bgc_ids <- unique(data.table(bgc = bgcs$BGC, bgc_id = bgcs$bgc_id))
    bc_bgc <- rasterize(bgcs, xyz, field = "bgc_id")
    return(list(bgc_rast = bc_bgc, ids = bgc_ids))
  } else {
    bgcs$id <- 1:nrow(bgcs)
    temp_pts <- vect(xyz, crs = "epsg:4326")
    bgcs_small <- crop(bgcs, temp_pts)
    temp_r <- rast(bgcs_small, res = res)
    bgc_id <- rasterize(bgcs_small, temp_r, field = "id")
    pnt_id <- extract(bgc_id, temp_pts)
    
    res <- data.table(cell = xyz$id, BGC = bgcs$BGC[pnt_id$id])
    return(res)
  }
  
}

#' Creates BGC predictions from RF model
#' @param dbCon Database connection
#' @param xyz SpatRaster or data.table.
#' @param BGCmodel Ranger random forest model of BGCs
#' @param vars_needed Character. List of variables required for model
#' @param gcms_use Character. List of gcms used in summarised predictions
#' @param periods_use Character. List of time periods to create predictions for
#' @param ssp_use Character. List of ssps to use. Default `c("ssp126", "ssp245", "ssp370")`
#' @return NULL. Results are written to csv files in base_folder/bgc_data
#' @import climr data.table ranger duckdb
#' @export
predict_bgc <- function(dbCon,
                              xyz, 
                              BGCmodel, 
                              vars_needed, 
                              gcms_use, 
                              periods_use, 
                              ssp_use = c("ssp126", "ssp245", "ssp370"),
                              max_runs_use = 0L,
                              start_tile = 1) {
  if(duckdb_table_exists(dbCon, "bgc_raw")) {
    message("Use cached table bgc_raw :)")
    return(invisible(TRUE))
  }
  
  if(inherits(xyz, "SpatRaster")){
    points_dat <- as.data.frame(xyz, cells=T, xy=T)
    colnames(points_dat) <- c("id", "lon", "lat", "elev")
    #points_dat <- points_dat[,c(2,3,4,1)] #restructure for climr input
  } else if(!all(c("lon", "lat", "elev", "id") %in% names(xyz))){
    stop("xyz must have columns lon, lat, elev, and id if it is a dataframe")
  } else {
    points_dat <- copy(xyz)
  }
  
  splits <- c(seq(1, nrow(points_dat), by = 10000), nrow(points_dat) + 1)
  message("There are ", length(splits), " tiles")
  
  for (i in start_tile:(length(splits) - 1)){
    cat(i, "\n")
    clim_dat <- climr::downscale(points_dat[splits[i]:(splits[i+1]-1),], 
                          which_refmap = "refmap_climr",
                          gcms = gcms_use,
                          gcm_periods = periods_use,
                          ssps = ssp_use,
                          max_run = max_runs_use,
                          vars = c(vars_needed, "MAT"),
                          nthread = 6,
                          return_refperiod = FALSE)
    addVars(clim_dat)
    clim_dat <- na.omit(clim_dat)
    
    mat_dat <- clim_dat[,.(cellnum = id, ssp = SSP, gcm = GCM, run = RUN, period = PERIOD, MAT)]
    dbWriteTable(dbCon, "clim_raw", mat_dat, row.names = FALSE, append = TRUE)
    
    setnames(clim_dat, old = c("PAS_an","Tmin_an","CMI_an"), new = c("PAS","Tmin","CMI"))
    temp <- predict(BGCmodel, data = clim_dat, num.threads = 8)
    dat <- data.table(cellnum = clim_dat$id, ssp = clim_dat$SSP, gcm = clim_dat$GCM, run = clim_dat$RUN,
                      period = clim_dat$PERIOD, bgc_pred = temp$predictions)
    dbWriteTable(dbCon, "bgc_raw", dat, row.names = FALSE, append = TRUE)

    rm(clim_dat, dat, mat_dat)
    gc()
  }
  message("Created table bgc_raw and clim_raw")
  
  ref_clim <- climr::downscale(points_dat, which_refmap = "refmap_climr", vars = "MAT", return_refperiod = TRUE)
  ref_clim[,PERIOD := NULL]
  setnames(ref_clim, c("cellnum","MAT"))
  dbWriteTable(dbCon, "clim_refperiod", ref_clim, row.names = FALSE)
  
  qry <- "CREATE TABLE clim_summary AS
                      
                      WITH clim_diff AS (
                        select a.cellnum, ssp, gcm, run, period, (a.MAT - clim_refperiod.MAT) as MAT_diff
                        FROM clim_raw a
                        JOIN clim_refperiod USING (cellnum)
                      )
                      SELECT
                        ssp, gcm, run, period, AVG(MAT_diff) as MAT_diff
                      FROM clim_diff
                      GROUP BY ssp, gcm, run, period;"
  dbExecute(dbCon, qry)
  message("Created table clim_refperiod and clim_summary")
}



# summary_preds_obs <- function(raster_template, 
#                               BGCmodel, 
#                               vars_needed, 
#                               out_folder = "bgc_data") {
#   points_dat <- as.data.frame(raster_template, cells=T, xy=T)
#   colnames(points_dat) <- c("id", "lon", "lat", "elev")
#   points_dat <- points_dat[,c(2,3,4,1)] #restructure for climr input
#   splits <- c(seq(1, nrow(points_dat), by = 500000), nrow(points_dat) + 1)
#   message("There are ", length(splits), " tiles")
#   
#   period_curr <- list_obs_periods()
#   for (i in 1:(length(splits) - 1)){
#     cat(i, "\n")
#     clim_dat <- downscale(points_dat[splits[i]:(splits[i+1]-1),], 
#                           which_refmap = "refmap_climr",
#                           obs_periods = period_curr,
#                           vars = vars_needed,
#                           nthread = 6,
#                           return_refperiod = FALSE)
#     addVars(clim_dat)
#     clim_dat <- na.omit(clim_dat)
#     temp <- predict(BGCmodel, data = clim_dat, num.threads = 16)
#     dat <- data.table(cellnum = clim_dat$id,  period = clim_dat$PERIOD, bgc_pred = temp$predictions, bgc_prop = 1)
#     fwrite(dat, append = TRUE, paste0(out_folder,"/bgc_summary_obs", ".csv"))
#     rm(clim_dat, dat)
#     gc()
#   }
#   message("Done!")
# }

#' Create siteseries predictions from summarised BGC predictions
#' @param dbCon Database connection to duckdb
#' @import data.table duckdb
#' @export
siteseries_preds <- function(dbCon,
                             obs = FALSE) {
  if(duckdb_table_exists(dbCon, "siteseries_preds")) {
    message("✓ Using cached table siteseries_preds")
    return(invisible(TRUE))
  }
  
  bgc_all <- dbGetQuery(dbCon, "select * from bgc_summary") |> as.data.table()
  bgc_points <- dbGetQuery(dbCon, "select * from bgc_points") |> as.data.table()
  
  periods <- unique(bgc_all$period)
  
  eda_all <- dbGetQuery(dbCon, "select * from edatopic") |> as.data.table()
  edatopes <- unique(eda_all$Edatopic)
  
  for(period_curr in periods){
    bgc_sum <- bgc_all[period == period_curr,]
    bgc_sum[bgc_points, BGC := i.bgc, on = "cellnum"]
    setcolorder(bgc_sum, c("cellnum","period","BGC","bgc_pred","bgc_prop"))
    setnames(bgc_sum, c("SiteRef","FuturePeriod","BGC","BGC.pred","BGC.prop"))
    
    for(edatope in edatopes){
      message(period_curr,", ", edatope)
      eda_table <- copy(eda_all)
      eda_table[,HasPos := if(any(Edatopic %in% edatope)) T else F, by = .(SS_NoSpace)]
      eda_table <- unique(eda_table[(HasPos),])
      ###########################################################
      
      sites <- unique(bgc_sum$SiteRef)
      splits <- c(seq(1, length(sites), by = 200000), length(sites) + 1)
      for (i in 1:(length(splits) - 1)){
        srs <- sites[splits[i]:(splits[i+1]-1)]
        dat_sml <- bgc_sum[SiteRef %in% srs,]
        sspred <- edatopicOverlap_fast(dat_sml, E1 = eda_table)
        sspred[,Edatope := edatope]
        dbWriteTable(dbCon, "siteseries_preds", sspred, append = TRUE, row.names = FALSE)
        rm(sspred)
        gc()
      }
    }
  }
  message("✓ Created table siteseries_preds !")
}

#' Create projected suitability values from site series predictions.
#' @param dbCon duckdb database connection
#' @param species Character vector. Species codes to create projections for. 
#' @param tile_size Integer. Number of sites to process at once. May need to decrease if memory is limited. Default 4000
#' @return NULL. Writes table to database
#' @import data.table duckdb
#' @export
cciss_suitability <- function(dbCon,
                              species,
                              obs = FALSE,
                              tile_size = 4000) {
  if(duckdb_table_exists(dbCon, "cciss_res")) {
    spp <- dbGetQuery(dbCon, "select distinct Spp from cciss_res")$Spp
    missing_spp <- setdiff(species, spp)
    if(length(missing_spp) == 0){
      message("All requested species already cached :)")
      return(invisible(TRUE))
    } else {
      message("Calculating CCISS for ", paste(missing_spp, collapse = ", "))
      species <- missing_spp
    }
  }
  
  feas_table <- dbGetQuery(dbCon, "select * from suitability") |> as.data.table()
  setnames(feas_table, c("BGC", "Spp","SS_NoSpace", "Feasible"))
  edatopes <- dbGetQuery(dbCon, "select distinct Edatopic from edatopic")$Edatopic
  periods <- dbGetQuery(dbCon, "select distinct FuturePeriod from siteseries_preds")$FuturePeriod
 # stopifnot(all(c("BGC","SS_NoSpace","Sppsplit","FeasOrig","Spp","Feasible","Mod","OR") %in% names(feas_table)))
  for(period in periods){
    for(edatope in edatopes){
      sspreds <- dbGetQuery(dbCon, sprintf("select * from siteseries_preds where FuturePeriod = '%s' AND Edatope = '%s'", period, edatope)) |> as.data.table()
      sitenums <- unique(sspreds$SiteRef)
      splits <- c(seq(1, length(sitenums), by = tile_size), length(sitenums) + 1)
      for(spp in species){
        message(period, " ", edatope, " ", spp)
        for (i in 1:(length(splits) - 1)){
          temp <- sspreds[SiteRef %in% sitenums[splits[i]:(splits[i+1]-1)],]
          cciss_res <- cciss_full(temp, feas_table, spp)
          cciss_res <- na.omit(cciss_res, cols = "SiteRef")
          cciss_res[,Edatope := edatope]
          dbWriteTable(dbCon, "cciss_res", cciss_res, row.names = FALSE, append = TRUE)
        }
        rm(cciss_res)
        gc()
      } 
    }
    rm(sspreds)
    gc()
  }
  message("✓ Created or updated table cciss_res !")
}

#' Create geotif rasters of projected suitabilities for each species/edatope/period
#' @param raster_template Template SpatRaster. Must be the same raster used throughout process
#' @param base_folder Base folder to write results to.
#' @return NULL. Writes result geotifs to base_folder/final_rasters
#' @import data.table terra
#' @export
cciss_rasterize <- function(raster_template, base_folder = "spatial") {
  in_folder <- file.path(base_folder,"cciss_suit")
  if(!dir.exists(paste0(base_folder,"/final_rasters"))) dir.create(paste0(base_folder,"/final_rasters"))
  out_folder <- paste0(base_folder,"/final_rasters")
  
  fnms <- list.files(in_folder)
  vartab <- as.data.table(tstrsplit(fnms, split = "_"))
  vartab[,V1 := NULL]
  setnames(vartab, c("P1","P2","Eda"))
  vartab[,Period := paste0(P1,"_",P2)]
  vartab[,Eda := gsub(".csv","",Eda)]
  
  for(period in unique(vartab$Period)){
    for(eda in unique(vartab$Eda)) {
      dat <- fread(paste0(in_folder,"/CCISS_",period,"_",eda,".csv"))
      dat <- na.omit(dat, cols = "SiteRef")
      spps <- unique(dat$Spp)
      for(spp in spps){
        dat_spp <- dat[Spp == spp,]
        dat_spp <- dat_spp[Curr < 3.5 | Newsuit < 3.5,]
        rfinal <- copy(raster_template)
        values(rfinal) <- NA
        rfinal[dat_spp$SiteRef] <- dat_spp$Newsuit
        writeRaster(rfinal,file.path(out_folder,paste0("CCISS_",period,"_", eda,"_",spp,".tif")), overwrite=TRUE)
      }
    }
  }
}

#' Function to create rasters of historic (mapped) suitability by species and edatopic position
#' @param species Character. Vector of species codes to map
#' @param edatopes Character. Vector of edatopes (e.g., "C4")
#' @param bgc_template_list List containing SpatRaster of BGCs and id table. Usually created using `make_bgc_template`
#' @param base_folder Base folder to write results to.
#' @return NULL. Writes tifs to "base_folder/historic_suit"
#' @export
mapped_suit <- function(species, edatopes, bgc_template_list, base_folder) {
  if(!dir.exists(paste0(base_folder,"/historic_suit"))) dir.create(paste0(base_folder,"/historic_suit"))
  out_folder <- paste0(base_folder,"/historic_suit")
  for(eda_sel in edatopes){
    for(spp_sel in species) {
      bgc_ids <- copy(bgc_template_list$ids)
      bgc_rast <- copy(bgc_template_list$bgc_rast)
      eda_sub <- copy(ccissr::E1)[Edatopic == eda_sel & is.na(SpecialCode),]
      suit_sub <- copy(ccissr::S1)[spp == spp_sel,]
      eda_sub[suit_sub, suit := i.newfeas, on = c(SS_NoSpace = "ss_nospace")]
      eda_sub[is.na(suit), suit := 5]
      eda_sub <- eda_sub[,.(suit = mean(suit)), by = .(BGC)]
      bgc_ids[eda_sub, suit := i.suit, on = c(bgc = "BGC")]
      suit_rast <- subst(bgc_rast, from = bgc_ids$bgc_id, to = bgc_ids$suit)
      writeRaster(suit_rast, file.path(out_folder, paste0("HistoricSuit_",spp_sel,"_",eda_sel,".tif")))
    }
  }
}


#' Internal function to create cciss projections from site series predictions. Calculates projected suitability for each period, as well as other CCISS statistics.
#' @param SSPred Data.table of site series predictions. Usually comes from `siteseries_preds`
#' @param suit Data.table. Environmental suitability table
#' @param spp_select Character vector of species to use.
#' @return data.table containing cciss results by siteno, period, edatope and species
#' @import data.table
#' @export
cciss_full <- function(SSPred,suit,spp_select){
  suit <- suit[Spp %in% spp_select,.(BGC,SS_NoSpace,Spp,Feasible)]
  suit <- unique(suit)
  suit <- na.omit(suit)
  SSPred <- SSPred[!grepl("[0-9]a$|[0-9]b$|[0-9]c$",SS_NoSpace),] ##remove phases
  SSPred <- SSPred[!grepl("\\.1$|\\.2$|\\.3$",SS_NoSpace),]
  SSPred <- SSPred[,.(SiteRef,FuturePeriod,BGC,SS_NoSpace,SS.pred,SSprob)]
  Site_BGC <- unique(SSPred[,.(SiteRef,BGC)])
  SSPred <- na.omit(SSPred)
  setkey(SSPred,SS.pred)
  setkey(suit,SS_NoSpace)
  suitMerge <- suit[SSPred, allow.cartesian = T]
  #suitMerge <- na.omit(suitMerge)
  suitMerge[is.na(Feasible), `:=`(Feasible = 5,
                                  Spp = spp_select)]
  setnames(suitMerge, old = c("SS_NoSpace", "i.SS_NoSpace"), new = c("SS.pred", "SS_NoSpace"))
  suitVotes <- data.table::dcast(suitMerge, SiteRef + Spp + FuturePeriod + SS_NoSpace ~ Feasible, 
                                 value.var = "SSprob", fun.aggregate = sum)
  # Fill with 0 if columns does not exist, encountered the error at SiteRef 3104856 
  colNms <- c("1","2","3","X")
  set(suitVotes, j = as.character(1:5)[!as.character(1:5) %in% names(suitVotes)], value = 0)
  
  suitVotes[,VoteSum := `1`+`2`+`3`+`4`+`5`]
  suitVotes[,X := 1 - VoteSum]
  suitVotes[,VoteSum := NULL]
  suitVotes[,X := X + `5` + `4`]
  suitVotes[,`:=`(`5` = NULL, `4` = NULL)]
  setkey(suitVotes, SS_NoSpace, Spp)
  #setkey(suit, SS_NoSpace, Spp)
  temp <- suit[,.(SS_NoSpace,Feasible)]
  setnames(temp, c("SS_NoSpace","Curr"))
  suitVotes <- merge(suitVotes, temp, by = "SS_NoSpace", all = T, allow.cartesian = TRUE)
  suitVotes[is.na(Curr), Curr := 5]
  setorder(suitVotes,SiteRef,SS_NoSpace,Spp,FuturePeriod)
  suitVotes[Curr > 3.5, Curr := 4]
  
  suitVotes[,Improve := ModelDir(as.matrix(.SD), Curr = Curr, dir = "Improve"),.SDcols = colNms]
  suitVotes[,Decline := ModelDir(as.matrix(.SD), Curr = Curr, dir = "Decline"),.SDcols = colNms]
  datRot <- suitVotes[,lapply(.SD, mean),.SDcols = c("Improve","Decline"), by = list(SiteRef,SS_NoSpace,FuturePeriod,Spp,Curr)]
  datRot[,`:=`(Improve = round(Improve*100),Decline = round(Decline*100))]
  datRot[,Curr := NULL]
  
  suitVotes <- suitVotes[,lapply(.SD, sum),.SDcols = colNms, 
                         by = .(SiteRef,FuturePeriod, SS_NoSpace,Spp,Curr)]
  suitVotes[,Newsuit := `1`+(`2`*2)+(`3`*3)+(X*5)]
  suitVotes <- merge(suitVotes, datRot, by = c('SiteRef','FuturePeriod','SS_NoSpace','Spp'),all = T)
  suitRes <- suitVotes[,.(Curr = mean(Curr),Newsuit = mean(Newsuit), Improve = mean(Improve), Decline = mean(Decline), Prop1 = mean(`1`), Prop2 = mean(`2`), Prop3 = mean(`3`)), by = .(SiteRef,FuturePeriod,Spp)]
  return(suitRes)
}
