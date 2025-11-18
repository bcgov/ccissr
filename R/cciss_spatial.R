
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

#' Create summarised BGC predictions from RF model
#' @param xyz SpatRaster or data.table.
#' @param BGCmodel Ranger random forest model of BGCs
#' @param vars_needed Character. List of variables required for model
#' @param gcms_use Character. List of gcms used in summarised predictions
#' @param periods_use Character. List of time periods to create predictions for
#' @param ssp_use Character. List of ssps to use. Default `c("ssp126", "ssp245", "ssp370")`
#' @param ssp_w Numeric vector. Weights for each ssp in `ssp_use`
#' @param base_folder Character. Name of base folder to write results to.
#' @return NULL. Results are written to csv files in base_folder/bgc_data
#' @import climr data.table
#' @export
summary_preds_gcm <- function(xyz, 
                              BGCmodel, 
                              vars_needed, 
                              gcms_use, 
                              periods_use, 
                              ssp_use = c("ssp126", "ssp245", "ssp370"),
                              ssp_w = c(0.8,1,0.8),
                              base_folder = "spatial",
                              start_tile = 1) {
  
  if(!dir.exists(paste0(base_folder,"/bgc_data"))) dir.create(paste0(base_folder,"/bgc_data"))
  out_folder <- paste0(base_folder,"/bgc_data")
  
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
  ssp_weights <- data.table(ssp = ssp_use, weight = ssp_w)
  message("There are ", length(splits), " tiles")
  
  for (i in start_tile:(length(splits) - 1)){
    cat(i, "\n")
    clim_dat <- climr::downscale(points_dat[splits[i]:(splits[i+1]-1),], 
                          which_refmap = "refmap_climr",
                          gcms = gcms_use,
                          gcm_periods = periods_use,
                          ssps = ssp_use,
                          max_run = 0L,
                          vars = vars_needed,
                          nthread = 6,
                          return_refperiod = FALSE)
    addVars(clim_dat)
    clim_dat <- na.omit(clim_dat)
    setnames(clim_dat, old = c("PAS_an","Tmin_an","CMI_an"), new = c("PAS","Tmin","CMI"))
    temp <- predict(BGCmodel, data = clim_dat, num.threads = 8)
    dat <- data.table(cellnum = clim_dat$id, ssp = clim_dat$SSP, gcm = clim_dat$GCM, 
                      period = clim_dat$PERIOD, bgc_pred = temp$predictions)
    dat[ssp_weights, weight := i.weight, on = "ssp"]
    dat_sum <- dat[,.(bgc_prop = sum(weight)/20.8), by = .(cellnum, period, bgc_pred)]
    
    fwrite(dat_sum, paste0(out_folder, "/bgc_summary_",i, ".csv"), append = TRUE)
    rm(clim_dat, dat, dat_sum)
    gc()
  }
  cat("Done!")
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
#' @param edatopes Character. Vector of desired edatopic positions. Default is "B2", "C4", "D6" (poor, mesic, rich).
#' @param obs Logical. Do prediction for observed period? Default FALSE
#' @param bgc_mapped List containing SpatRaster of BGCs and id table, or data.table (must have columns `cell`,`BGC`). Usually created using `make_bgc_template`
#' @param base_folder Base folder to write results to.
#' @return NULL. Writes results to csvs in base_folder/ss_preds
#' @import data.table
#' @export

siteseries_preds <- function(edatopes = c("B2", "C4", "D6"), 
                             obs = F, 
                             bgc_mapped,
                             base_folder = "spatial") {
  in_folder <- file.path(base_folder,"bgc_data")
  if(obs){
    periods <- list_obs_periods()
    bgc_all <- fread(file.path(in_folder, "bgc_summary_obs.csv"))
    obs_nm <- "obs_"
  } else {
    temp_ls <- list.files(in_folder, full.names = TRUE)
    temp_ls <- temp_ls[!grepl("obs",temp_ls)]
    bgc_all_ls <- lapply(temp_ls, FUN = fread)
    bgc_all <- rbindlist(bgc_all_ls)
    rm(bgc_all_ls)
    obs_nm <- ""
  }
  
  periods <- unique(bgc_all$period)
  
  if(inherits(bgc_mapped, "list")){
    bgc_rast <- bgc_mapped$bgc_rast
    rast_ids <- bgc_mapped$ids
    bgc_points <- as.data.frame(bgc_rast, cells=T, xy=T) |> as.data.table()
    bgc_points[rast_ids, BGC := i.bgc, on = "bgc_id"]
  } else {
    bgc_points <- bgc_mapped |> as.data.table()
  }

  
  if(!dir.exists(paste0(base_folder,"/ss_preds"))) dir.create(paste0(base_folder,"/ss_preds"))
  out_folder <- paste0(base_folder,"/ss_preds")
  
  # eda_ss <- special_ss[,.(SS_NoSpace,SpecialCode)]
  # edatopic[eda_ss, SpecialCode := i.SpecialCode, on = "SS_NoSpace"]
  # eda_all <- edatopic[is.na(SpecialCode),]
  eda_all <- copy(ccissr::E1)
  
  for(period_curr in periods){
    bgc_sum <- bgc_all[period == period_curr,]
    bgc_sum[bgc_points, BGC := i.BGC, on = c(cellnum = "cell")]
    setcolorder(bgc_sum, c("cellnum","period","BGC","bgc_pred","bgc_prop"))
    setnames(bgc_sum, c("SiteRef","FuturePeriod","BGC","BGC.pred","BGC.prop"))
    
    for(edatope in edatopes){
      cat(period_curr, edatope, "\n")
      eda_table <- copy(eda_all)
      eda_table[,HasPos := if(any(Edatopic %in% edatope)) T else F, by = .(SS_NoSpace)]
      eda_table <- unique(eda_table[(HasPos),])
      ###########################################################
      
      sites <- unique(bgc_sum$SiteRef)
      splits <- c(seq(1, length(sites), by = 200000), length(sites) + 1)
      for (i in 1:(length(splits) - 1)){
        cat(i, "\n")
        srs <- sites[splits[i]:(splits[i+1]-1)]
        dat_sml <- bgc_sum[SiteRef %in% srs,]
        sspred <- edatopicOverlap_fast(dat_sml, E1 = eda_table)
        fwrite(sspred, append = TRUE, paste0(out_folder, "/siteseries_",obs_nm,period_curr, "_", edatope, ".csv"))
        rm(sspred)
        gc()
      }
    }
  }
  message("Done!")
}

#' Create projected suitability values from site series predictions.
#' @param species Character vector. Species codes to create projections for. 
#' @param edatopes Character vector. List of edatopes for projections. Must be a subset of the edatopes run in `siteseries_preds`
#' @param periods Character. Periods to run projections for. Default is `list_gcm_periods()`
#' @param base_folder Base folder to write results to.
#' @param tile_size Integer. Number of sites to process at once. May need to decrease if memory is limited. Default 4000
#' @return NULL. Writes results to csvs in base_folder/cciss_suit
#' @import data.table
#' @export
cciss_suitability <- function(species,
                              edatopes,
                              obs = FALSE,
                              periods = list_gcm_periods(),
                              base_folder = "spatial",
                              tile_size = 4000) {
  feas_table <- copy(ccissr::S1)
  setnames(feas_table, c("BGC","SS_NoSpace","Sppsplit","FeasOrig","Spp","Feasible","Mod","OR"))
 # stopifnot(all(c("BGC","SS_NoSpace","Sppsplit","FeasOrig","Spp","Feasible","Mod","OR") %in% names(feas_table)))
  
  in_folder <- file.path(base_folder,"ss_preds")
  if(!dir.exists(paste0(base_folder,"/cciss_suit"))) dir.create(paste0(base_folder,"/cciss_suit"))
  out_folder <- paste0(base_folder,"/cciss_suit")
  
  if(obs) {
    periods <- list_obs_periods()
    obs_nm <- "obs_"
  } else {
    obs_nm <- ""
  }
  
  for(period in periods){
    for(edatope in edatopes){
      sspreds <- fread(paste0(in_folder,"/siteseries_",obs_nm,period,"_",edatope,".csv"))
      sitenums <- unique(sspreds$SiteRef)
      splits <- c(seq(1, length(sitenums), by = tile_size), length(sitenums) + 1)
      for(spp in species){
        message(period, edatope, spp, "\n")
        for (i in 1:(length(splits) - 1)){
          temp <- sspreds[SiteRef %in% sitenums[splits[i]:(splits[i+1]-1)],]
          cciss_res <- cciss_full(temp, feas_table, spp)
          cciss_res <- na.omit(cciss_res, cols = "SiteRef")
          fwrite(cciss_res, append = TRUE, paste0(out_folder,"/CCISS_",obs_nm,period,"_",edatope,".csv"))
        }
        rm(cciss_res)
        gc()
      } 
    }
    rm(sspreds)
    gc()
  }
  cat("Done!")
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
