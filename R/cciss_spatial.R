
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
#' @param obs_2001_2020 Logical. Also predict for 2001 - 2020 observed climate? Default `FALSE`
#' @param refperiod Logical. Also predict for reference period (1961_1990)? Default `FALSE`
#' @param max_runs_use Integer. Number of individual runs for each GCM. Default `0` (only use ensembleMean)
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
                        obs_2001_2020 = FALSE,
                        refperiod = FALSE,
                        max_runs_use = 0L,
                        start_tile = 1) {
  periods_needed <- if(obs_2001_2020) c(periods_use,"2001_2020_obs") else periods_use
  periods_needed <- if(refperiod) c(periods_needed,"1961_1990") else periods_needed
  if(duckdb_table_exists(dbCon, "bgc_raw")) {
    periods_cached <- dbGetQuery(dbCon, "select distinct period from bgc_raw")$period
    if(all(periods_needed %in% periods_cached)){
      message("Use cached table bgc_raw :)")
      return(invisible(TRUE))
    } else {
      periods_needed <- setdiff(periods_needed, periods_cached)
      message("Will predict missing period ", periods_needed)
    }
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
  if("2001_2020_obs" %in% periods_needed) obs <- "2001_2020" else obs <- NULL
  if("1961_1990" %in% periods_needed) refperiod <- TRUE else FALSE
  
  periods_needed <- periods_needed[!periods_needed %in% c("2001_2020_obs", "1961_1990")]
  if(length(periods_needed) < 1) periods_needed <- gcms_use <- ssp_use <- NULL
  tmp_names <- data.table(id = numeric(), GCM = character(), SSP = character(), RUN = character())
  
  for (i in start_tile:(length(splits) - 1)){
    cat(i, "\n")
    clim_dat <- climr::downscale(points_dat[splits[i]:(splits[i+1]-1),], 
                          which_refmap = "refmap_climr",
                          gcms = gcms_use,
                          gcm_periods = periods_needed,
                          obs_periods = obs,
                          ssps = ssp_use,
                          max_run = max_runs_use,
                          vars = c(vars_needed, "MAT"),
                          nthread = 6,
                          ensemble_mean = max_runs_use < 1,
                          return_refperiod = refperiod)
    addVars(clim_dat)
    clim_dat <- na.omit(clim_dat)
    clim_dat <- rbind(clim_dat, tmp_names, use.names = TRUE, fill = TRUE)
    clim_dat[PERIOD == "2001_2020" & is.na(GCM), PERIOD := "2001_2020_obs"]
    
    mat_dat <- clim_dat[,.(cellnum = id, ssp = SSP, gcm = GCM, run = RUN, period = PERIOD, MAT)]
    dbWriteTable(dbCon, "clim_raw", mat_dat, row.names = FALSE, append = TRUE)
    
    temp <- predict(BGCmodel, data = clim_dat, num.threads = 8)
    dat <- data.table(cellnum = clim_dat$id, ssp = clim_dat$SSP, gcm = clim_dat$GCM, run = clim_dat$RUN,
                      period = clim_dat$PERIOD, bgc_pred = temp$predictions)
    dbWriteTable(dbCon, "bgc_raw", dat, row.names = FALSE, append = TRUE)

    rm(clim_dat, dat, mat_dat)
    gc()
  }
  message("Created/updated table bgc_raw and clim_raw")
  
  if(!dbExistsTable(dbCon,"clim_refperiod")) {
    ref_clim <- climr::downscale(points_dat, which_refmap = "refmap_climr", vars = "MAT", return_refperiod = TRUE)
    ref_clim[,PERIOD := NULL]
    setnames(ref_clim, c("cellnum","MAT"))
    dbWriteTable(dbCon, "clim_refperiod", ref_clim, row.names = FALSE)
  }
  
  dbExecute(dbCon, "DROP TABLE IF EXISTS clim_summary")
  
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
  
  bgc_qry <- "select * from bgc_summary"
  if(duckdb_table_exists(dbCon, "siteseries_preds")) {
    periods_raw <- dbGetQuery(dbCon, "select distinct period from bgc_summary")$period
    periods_perexp <- dbGetQuery(dbCon, "select distinct FuturePeriod from siteseries_preds")$FuturePeriod
    missing <- setdiff(periods_raw,periods_perexp)
    
    if(length(missing) == 0) {
      message("✓ Using cached table siteseries_preds")
      return(invisible(TRUE))
    }
    message("Updating database for missing periods: ", missing)
    bgc_qry <- paste0("select * from bgc_summary where period in ('", paste(missing, collapse = "','"), "')")
  }
  
  bgc_all <- dbGetQuery(dbCon, bgc_qry) |> as.data.table()
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
  # periods we *could* run (from siteseries_preds)
  periods <- dbGetQuery(
    dbCon,
    "select distinct FuturePeriod from siteseries_preds"
  )$FuturePeriod
  
  # edatopes we care about
  edatopes_raw <- dbGetQuery(
    dbCon,
    "select distinct Edatopic from edatopic"
  )$Edatopic
  
  # we'll call the combo column Edatope to match cciss_res
  edatopes <- edatopes_raw
  
  # Build the full set of requested combos: period × edatope × species
  # (only those you *want* to compute)
  full_combos <- data.table::CJ(
    FuturePeriod = periods,
    Edatope      = edatopes,
    Spp          = species,
    unique       = TRUE
  )
  
  # Figure out which combos are already cached, if cciss_res exists
  if (duckdb_table_exists(dbCon, "cciss_res")) {
    existing_combos <- DBI::dbGetQuery(
      dbCon,
      "select distinct FuturePeriod, Edatope, Spp from cciss_res"
    ) |>
      data.table::as.data.table()
    
    # Anti-join: combos we still need to compute
    missing_combos <- full_combos[
      !existing_combos,
      on = c("FuturePeriod", "Edatope", "Spp")
    ]
    
    if (nrow(missing_combos) == 0L) {
      message("All requested species/period/edatope combinations already cached :)")
      return(invisible(TRUE))
    }
    
    message(
      "Found ", nrow(existing_combos), " cached combos, ",
      nrow(missing_combos), " remaining to compute."
    )
  } else {
    # Table doesn't exist yet: everything is missing
    missing_combos <- full_combos
    message(
      "Table cciss_res does not exist yet; will compute all ",
      nrow(missing_combos), " species/period/edatope combinations."
    )
  }
  
  # Suitability (feasibility) lookup
  feas_table <- DBI::dbGetQuery(dbCon, "select * from suitability") |>
    data.table::as.data.table()
  data.table::setnames(feas_table, c("BGC", "Spp", "SS_NoSpace", "Feasible"))
  
  # For efficiency, loop over unique (period, edatope) pairs that still
  # have at least one missing species; within those, loop over the
  # missing species only.
  missing_periods  <- unique(missing_combos$FuturePeriod)
  
  for (period in missing_periods) {
    missing_edatopes <- unique(
      missing_combos[FuturePeriod == period, Edatope]
    )
    
    for (edatope in missing_edatopes) {
      # Only species that are actually missing for this (period, edatope)
      spp_to_run <- missing_combos[
        FuturePeriod == period & Edatope == edatope,
        unique(Spp)
      ]
      
      # Pull all preds for this period/edatope once
      sspreds <- DBI::dbGetQuery(
        dbCon,
        sprintf(
          "select * from siteseries_preds
           where FuturePeriod = '%s' AND Edatope = '%s'",
          period, edatope
        )
      ) |>
        data.table::as.data.table()
      
      sitenums <- unique(sspreds$SiteRef)
      splits <- c(
        seq(1, length(sitenums), by = tile_size),
        length(sitenums) + 1
      )
      
      for (spp in spp_to_run) {
        message("Computing: ", period, " / ", edatope, " / ", spp)
        
        for (i in seq_len(length(splits) - 1L)) {
          idx <- splits[i]:(splits[i + 1L] - 1L)
          temp_sitenums <- sitenums[idx]
          
          if (length(temp_sitenums) == 0L) next
          
          temp <- sspreds[SiteRef %in% temp_sitenums]
          
          if (nrow(temp) == 0L) next
          
          cciss_res <- cciss_full(temp, feas_table, spp)
          cciss_res <- na.omit(cciss_res, cols = "SiteRef")
          cciss_res[, Edatope := edatope]
          
          DBI::dbWriteTable(
            dbCon,
            "cciss_res",
            cciss_res,
            row.names = FALSE,
            append   = TRUE
          )
        }
        
        rm(cciss_res)
        gc()
      }
      
      rm(sspreds)
      gc()
    }
  }
  
  message("✓ Created or updated table cciss_res !")
  invisible(TRUE)
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

#' Function to create create table of reference species suitability by BGC for given edatopic positions
#' @param dbCon duckdb database connection
#' @param species Character. Vector of species codes to map
#' @param edatopes Character. Vector of edatopes (e.g., "C4")
#' @importFrom glue glue_sql
#' @importFrom data.table as.data.table
#' @import duckdb
#' @return Data.table of species x edatope x bgc suitabilities
#' @export
ref_suit <- function(dbCon, species, edatopes) {
  qry <- glue_sql("
                  WITH
                  -- 1) compute suit per (bgc, edatope, spp) from edatopic + suitability
                  eda_suit AS (
                    SELECT
                      e.BGC AS bgc,
                      e.Edatopic,
                      s.spp,
                      CASE WHEN s.newfeas IS NULL OR s.newfeas = 4 THEN 5 ELSE s.newfeas END AS suit
                    FROM edatopic e
                    LEFT JOIN suitability s
                      ON e.SS_NoSpace = s.ss_nospace
                     AND s.spp IN ({species*})
                    WHERE e.Edatopic IN ({edatopes*})
                  ),
                  
                  -- 2) average suit per BGC (mirrors your data.table mean)
                  suit_mean AS (
                    SELECT
                      bgc,
                      Edatopic,
                      spp,
                      AVG(suit) AS suit
                    FROM eda_suit
                    GROUP BY bgc, Edatopic, spp
                  ),
                  
                  -- 3) all BGCs you want to output
                  bgc_units AS (
                    SELECT DISTINCT bgc
                    FROM bgc_points
                  ),
                  
                  -- 4) build full grid: (bgc × edatope × spp)
                  grid AS (
                    SELECT
                      b.bgc,
                      e.Edatopic,
                      s.spp
                    FROM bgc_units b
                    CROSS JOIN (SELECT UNNEST([{edatopes*}]) AS Edatopic) e
                    CROSS JOIN (SELECT UNNEST([{species*}])  AS spp) s
                  )
                  
                  -- 5) left join computed suits onto full grid; fill missing with 5
                  SELECT
                    g.bgc,
                    g.Edatopic,
                    g.spp,
                    COALESCE(sm.suit, 5) AS suit
                  FROM grid g
                  LEFT JOIN suit_mean sm
                    ON sm.bgc = g.bgc
                   AND sm.Edatopic = g.Edatopic
                   AND sm.spp = g.spp
", .con = dbCon)
  dat <- dbGetQuery(dbCon, qry) |> as.data.table()
  return(dat)
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


map_reference_suit <- function(con,
                              raster_template,
                              species,
                              edatope,
                              period) {
  final_dem <- copy(raster_template)
  breakpoints.suit <- c(1,2,3,999)
  palette.suit <-   c("#006400", "#1E90FF", "#EEC900", "#FFFFFF")
  breakpoints.change <- c(c(seq(-2.5,2.5,0.5),-10,10,20,30) + 15, 999)
  palette.change <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,5,6)], brewer.pal(11,"RdBu")[c(7,8,9,10,11)],"#000000", brewer.pal(9,"YlOrRd")[1:3],"#FFFFFF") # nolint
  breakpoints.binary <- seq(-1,1,0.2)
  palette.binary <- c(brewer.pal(11,"RdBu")[c(1:4,6,6)], brewer.pal(11,"RdBu")[c(6,8:11)])
  
  ##feas colours
  suit_cols <- data.table(value = breakpoints.suit,Colour = palette.suit)
  
  ##mean change colours
  change_cols <- data.table(value = breakpoints.change, Colour = palette.change)
  change_cols[value == 15, Colour := "#DFDFDF"]
  ##addret colours
  #addret_cols <- data.table(value = breakpoints.binary*100, Colour = palette.binary)
  
  if(obs) {
    periods <- list_obs_periods()
    obs_nm <- "obs_"
  } else {
    obs_nm <- ""
  }
  
  for(period in periods){
    for(edatope in edatopes){
      dat <- fread(paste0(in_folder,"/CCISS_",obs_nm,period,"_",edatope,".csv"))
      for(spp in species){
        cat(period, edatope, spp, "\n")
        dat_spp <- dat[Spp == spp,]
        dat_spp <- dat_spp[Curr < 3.5 | Newsuit < 3.5,]
        dat_spp[,FeasChange := Curr - Newsuit]
        dat_spp[Newsuit > 3.5 & Curr <= 3, FeasChange := -10]
        dat_spp[Curr > 3.5, FeasChange := round(FeasChange) * 10]
        dat_spp[,FeasChange := round(FeasChange/0.5)*0.5]
        dat_spp[,FeasRound := round(Newsuit)]
        dat_spp[,CurrRound := round(Curr)]
        dat_spp[CurrRound > 3, CurrRound := 999]
        dat_spp[FeasRound > 3, FeasRound := 999]
        dat_spp[,AddRet := Improve]
        dat_spp[Decline > Improve, AddRet := -Decline]
        dat_spp[,AddRet := round(AddRet/20)*20]
        
        #historic feasibility
        # if(period == "2001_2020" & !obs){
        #     values(final_dem) <- NA
        #     final_dem[!is.na(raster_template)] <- 999
        #     final_dem[dat_spp$SiteRef] <- dat_spp$CurrRound
        #     coltab(final_dem) <- suit_cols
        #     final_rgb <- colorize(final_dem, to = "rgb", alpha = TRUE)
        #     writeRaster(final_rgb, paste0(out_folder,"/HistoricFeas_",period,"_",edatope,"_",spp,".tif"), overwrite = T)
        # }
        
        ##new feasibility
        values(final_dem) <- NA
        final_dem[!is.na(raster_template)] <- 999
        final_dem[dat_spp$SiteRef] <- dat_spp$FeasRound
        coltab(final_dem) <- suit_cols
        final_rgb <- colorize(final_dem, to = "rgb", alpha = TRUE)
        writeRaster(final_rgb, paste0(out_folder,"/NewFeas_",obs_nm,period,"_",edatope,"_",spp,".tif"), overwrite = T)
        
        ## raw rasters
        trast <- copy(final_dem)
        values(trast) <- NA
        trast[dat_spp$SiteRef] <- dat_spp$FeasRound
        trast[trast == 999] <- NA
        trast <- as.int(trast * 10)
        writeRaster(trast, paste0(out_folder_raw,"/Feasibility_",obs_nm,period,"_",edatope,"_",spp,".tif"),overwrite = T, datatype = "INT2U")
        
        ##mean change
        values(final_dem) <- NA
        final_dem[!is.na(raster_template)] <- 999
        final_dem[dat_spp$SiteRef] <- dat_spp$FeasChange + 15
        final_rgb <- subst(final_dem, change_cols$value, t(col2rgb(change_cols$Colour,alpha = TRUE)),names = c("red","green", "blue","alpha"))
        writeRaster(final_rgb, paste0(out_folder,"/MeanChange_",obs_nm,period,"_",edatope,"_",spp,".tif"),overwrite = T)
        
        trast <- copy(final_dem)
        values(trast) <- NA
        trast[dat_spp$SiteRef] <- dat_spp$FeasChange
        trast[trast == 999] <- NA
        trast <- as.int(trast * 10)
        writeRaster(trast, paste0(out_folder_raw,"/MeanChange_",period,"_",edatope,"_",spp,".tif"),overwrite = T, datatype = "INT4S")
        
        gc()
      }
    }
  }
}

map_reference_suit <- function(con,
                               raster_template,
                               species,
                               edatope,
                               table_name = "cciss_res",
                               return_raw = FALSE) {
  final_dem <- copy(raster_template)
  breakpoints.suit <- c(1,2,3,999)
  palette.suit <-   c("#006400", "#1E90FF", "#EEC900", "#FFFFFF")
  breakpoints.change <- c(c(seq(-2.5,2.5,0.5),-10,10,20,30) + 15, 999)
  palette.change <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,5,6)], brewer.pal(11,"RdBu")[c(7,8,9,10,11)],"#000000", brewer.pal(9,"YlOrRd")[1:3],"#FFFFFF") # nolint
  breakpoints.binary <- seq(-1,1,0.2)
  palette.binary <- c(brewer.pal(11,"RdBu")[c(1:4,6,6)], brewer.pal(11,"RdBu")[c(6,8:11)])
  
  ##feas colours
  suit_cols <- data.table(value = breakpoints.suit,Colour = palette.suit)
  dat_spp <- dbGetQuery(con, glue_sql("select * from {table_name} where 
                                  Spp = {species} 
                                  and Edatope = {edatope} 
                                  and FuturePeriod = '2041_2060'", .con = con))
  
  setDT(dat_spp)
  dat_spp <- dat_spp[Curr < 3.5 | Newsuit < 3.5,]
  dat_spp[,FeasChange := Curr - Newsuit]
  dat_spp[Newsuit > 3.5 & Curr <= 3, FeasChange := -10]
  dat_spp[Curr > 3.5, FeasChange := round(FeasChange) * 10]
  dat_spp[,FeasChange := round(FeasChange/0.5)*0.5]
  dat_spp[,FeasRound := round(Newsuit)]
  dat_spp[,CurrRound := round(Curr)]
  dat_spp[CurrRound > 3, CurrRound := 999]
  dat_spp[FeasRound > 3, FeasRound := 999]
  
  values(final_dem) <- NA
  if(return_raw) {
    final_dem[!is.na(raster_template)] <- 5
    final_dem[dat_spp$SiteRef] <- dat_spp$Curr
    return(final_dem)
  }
  
  final_dem[!is.na(raster_template)] <- 999
  final_dem[dat_spp$SiteRef] <- dat_spp$CurrRound
  coltab(final_dem) <- suit_cols
  final_rgb <- colorize(final_dem, to = "rgb", alpha = TRUE)
  return(final_rgb)
}

map_change_suit <- function(con,
                               raster_template,
                               species,
                               edatope,
                               period, 
                           table_name = "cciss_res") {
  final_dem <- copy(raster_template)
  breakpoints.change <- c(c(seq(-2.5,2.5,0.5),-10,10,20,30) + 15, 999)
  palette.change <- c(brewer.pal(11,"RdBu")[c(1,2,3,4,5,6)], brewer.pal(11,"RdBu")[c(7,8,9,10,11)],"#000000", brewer.pal(9,"YlOrRd")[1:3],"#FFFFFF") # nolint
  
  ##mean change colours
  change_cols <- data.table(value = breakpoints.change, Colour = palette.change)
  change_cols[value == 15, Colour := "#DFDFDF"]
  
  dat_spp <- dbGetQuery(con, glue_sql("select * from {table_name} where 
                                  Spp = {species} 
                                  and Edatope = {edatope} 
                                  and FuturePeriod = {period}", .con = con))
  setDT(dat_spp)
  dat_spp <- dat_spp[Curr < 3.5 | Newsuit < 3.5,]
  dat_spp[,FeasChange := Curr - Newsuit]
  dat_spp[Newsuit > 3.5 & Curr <= 3, FeasChange := -10]
  dat_spp[Curr > 3.5, FeasChange := round(FeasChange) * 10]
  dat_spp[,FeasChange := round(FeasChange/0.5)*0.5]
  dat_spp[,FeasRound := round(Newsuit)]
  dat_spp[,CurrRound := round(Curr)]
  dat_spp[CurrRound > 3, CurrRound := 999]
  dat_spp[FeasRound > 3, FeasRound := 999]
  
  ##mean change
  values(final_dem) <- NA
  final_dem[!is.na(raster_template)] <- 999
  final_dem[dat_spp$SiteRef] <- dat_spp$FeasChange + 15
  final_rgb <- subst(final_dem, change_cols$value, t(col2rgb(change_cols$Colour,alpha = TRUE)),
                     names = c("red","green", "blue","alpha"))
  
  return(final_rgb)
}
