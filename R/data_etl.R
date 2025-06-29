###CCISS 2020 Step 1: Pull data from postgres database
###Kiri Daust

# Copyright 2021 Will Mackenzie
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Pull points info from lat long inputs
#' @param con An active postgres DBI connection.
#' @param points A data.frame with two columns representing longitude
#'  (Long) and latitude (Lat) with GPS values (EPSG:4326) 
#' @return BEC Information, Site Reference Hex Grid Number and elevation in meter
#' @importFrom RPostgres dbGetQuery
#' @export
dbPointInfo <- function(con, points) {
  
  # It is more efficient for the query optimizer to do
  # each join in a separate sql, better use of index and caching
  # geometry have been optimized (split polygons) for the
  # whole function to return under 50ms
  
  Long <- .subset2(points, "Long")
  Lat <- .subset2(points, "Lat")
  if(!"ID" %in% colnames(points)){
    ID <- seq_along(Long)
  } else {
    ID <- .subset2(points, "ID")
  }
  
  elev_info_sql <- paste0("
    WITH pts4269 AS (
      ", paste0("SELECT st_transform(st_pointfromtext('POINT(", Long, " ", Lat, ")', 4326), 4269) geom,
      '",ID,"' id ", collapse = "\n UNION ALL \n") ,"
    )
    
    SELECT pts.id, 
    MAX(ROUND(CAST(ST_Value(dem.rast, pts.geom) as NUMERIC), 2)) elevation_m
    FROM bc_elevation dem
    CROSS JOIN pts4269 pts
    WHERE ST_Intersects(dem.rast, pts.geom)
    GROUP BY pts.id
  ")
  
  bec_info_sql <- paste0("
    WITH pts3005 AS (
      ", paste0("SELECT st_transform(st_pointfromtext('POINT(", Long, " ", Lat, ")', 4326), 3005) geom,
      '",ID,"' id ", collapse = "\n UNION ALL \n") ,"
    )
    
    SELECT pts.id,
           bec.zone,
           bec.subzone,
           bec.variant,
           bec.phase,
           bec.natural_disturbance,
           bec.map_label,
           bec.bgc_label,
           bec.zone_name,
           bec.subzone_name,
           bec.variant_name,
           bec.phase_name,
           bec.natural_disturbance_name,
           bec.feature_area_sqm,
           bec.feature_length_m
    FROM pts3005 pts
    LEFT JOIN bec_info bec
    ON ST_Intersects(pts.geom, bec.geometry)
  ")
  
  site_ref_sql <- paste0("
    WITH pts3005 AS (
      ", paste0("SELECT st_transform(st_pointfromtext('POINT(", Long, " ", Lat, ")', 4326), 3005) geom,
      '",ID,"' id ", collapse = "\n UNION ALL \n") ,"
    )
    
    SELECT pts.id, 
    hgrid.siteno site_no
    FROM pts3005 pts 
    LEFT JOIN hex_grid hgrid
    ON ST_Intersects(pts.geom, hgrid.geom)
  ")
  
  bc_land_sql <- paste0("
    WITH pts3005 AS (
      ", paste0("SELECT st_transform(st_pointfromtext('POINT(", Long, " ", Lat, ")', 4326), 3005) geom,
      '",ID,"' id ", collapse = "\n UNION ALL \n") ,"
    )
    
    SELECT pts.id,
    bcb_hres.objectid is NOT NULL onbcland
    FROM pts3005 pts
    LEFT JOIN bcb_hres
    ON ST_Within(pts.geom, bcb_hres.geometry)
  ")
  
  bc_forest_region <- paste0("
    WITH pts3005 AS (
      ", paste0("SELECT st_transform(st_pointfromtext('POINT(", Long, " ", Lat, ")', 4326), 3005) geom,
      '",ID,"' id ", collapse = "\n UNION ALL \n") ,"
    )
    
    SELECT pts.id,
    region_ref forest_region
    FROM pts3005 pts
    LEFT JOIN bc_forest_regions
    ON ST_Within(pts.geom, bc_forest_regions.geom)                       
  ")
  
  res1 <- RPostgres::dbGetQuery(con, bec_info_sql) |> as.data.table()
  setkey(res1,id)
  res2 <- RPostgres::dbGetQuery(con, elev_info_sql)|> as.data.table()
  setkey(res2,id)
  res3 <- RPostgres::dbGetQuery(con, site_ref_sql)|> as.data.table()
  setkey(res3,id)
  res4 <- RPostgres::dbGetQuery(con, bc_land_sql)|> as.data.table()
  setkey(res4,id)
  res5 <- RPostgres::dbGetQuery(con, bc_forest_region)|> as.data.table()
  setkey(res5,id)
  res <- res1[res2,][res3,][res4,][res5,]
  return(res)  
}


# Use the database, it is faster than R native sf functions
# The buffer is to make sure points are not at the edge of the map
#' Buffered bounding box
#' @return a bounding box 1km around points
#' @inheritParams dbPointInfo
#' @param buffer A numeric. Buffer distance in km.
#' @export
dbBbox <- function(con, points, buffer) {
  Long <- .subset2(points, "Long")
  Lat <- .subset2(points, "Lat")
  bbox_sql <- paste0("
    WITH pts3005 AS (
      ", paste0("SELECT st_transform(st_pointfromtext('POINT(", Long, " ", Lat, ")', 4326), 3005) geom", collapse = "\n UNION ALL \n") ,"
    )
    
    SELECT st_xmin(box),
           st_ymin(box),
           st_xmax(box),
           st_ymax(box)
    FROM (
      SELECT ST_Extent(st_transform(st_buffer(pts.geom, ", buffer, ", 3), 4326)) box
      FROM pts3005 pts
    ) box
  ")
  unname(RPostgres::dbGetQuery(con, bbox_sql))
}

# drv <- dbDriver("Postgres")
# con <- dbConnect(drv, user = "postgres", password = "jujL6cB3Wm9y", host = "138.197.168.220", 
#                  port = 5432, dbname = "cciss")
# 
# sites <- 6305115:6305125
# test <- dbGetCCISS(pool, 6518114, FALSE, scn = "ssp370")

#' Get sitenos for preselected points in BGC/district
#' @param con postgres DBI connection
#' @param bgc A character specifying BGC
#' @param district A character specifying district or NULL
#' @param maxPoints Integer - max number of points
#' @details Get siteno for BGC and District.
#' @return a vector of sitenumbers
#' @importFrom RPostgres dbGetQuery
#' @export
dbGetBGC <- function(con,bgc,district = NULL,maxPoints){
  if(is.null(district)){
    query <- paste0("select siteno from preselected_points13 where bgc IN ('",paste(bgc,collapse = "','"),"')")
  }else{
    query <- paste0("select siteno from preselected_dist13 where bgc IN ('",paste(bgc,collapse = "','"),"') and dist_code = '",district,"'")
  }
  dat <- RPostgres::dbGetQuery(con, query)$siteno
  return(dat)
}

#' Pull CCISS from a vector of SiteNo
#' @param con An active postgres DBI connection.
#' @param siteno A character vector of siteno.
#' @param avg A boolean. 
#' @param modWeights A data table of gcm and rcp weights.
#' @details Get CCISS for provided SiteNo.
#' @return A data.table containing CCISS information for each provided SiteNo.
#' @importFrom RPostgres dbGetQuery
#' @export
dbGetCCISS <- function(con, siteno, avg, modWeights){

  # Declare binding for checks
  if (FALSE) {
    comb <- gcm <- rcp <- weight <- NULL
  }
  
  groupby = "siteno"
  if (isTRUE(avg)) {
    groupby = "bgc"
  }
  modWeights[,comb := paste0("('",gcm,"','",rcp,"',",weight,")")]
  weights <- paste(modWeights$comb,collapse = ",")
  
  ##cciss_future is now test_future  
  cciss_sql <- paste0("
  WITH cciss AS (
    SELECT cciss_future12_array.siteno,
         labels.gcm,
         labels.scenario,
         labels.futureperiod,
         bgc_attribution.bgc,
         bgc.bgc bgc_pred,
         w.weight
  FROM cciss_future12_array
  JOIN bgc_attribution
    ON (cciss_future12_array.siteno = bgc_attribution.siteno),
       unnest(bgc_pred_id) WITH ordinality as source(bgc_pred_id, row_idx)
  JOIN (SELECT ROW_NUMBER() OVER(ORDER BY gcm_id, scenario_id, futureperiod_id) row_idx,
               gcm,
               scenario,
               futureperiod
        FROM gcm 
        CROSS JOIN scenario
        CROSS JOIN futureperiod) labels
    ON labels.row_idx = source.row_idx
    JOIN (values ",weights,") 
    AS w(gcm,scenario,weight)
    ON labels.gcm = w.gcm AND labels.scenario = w.scenario
  JOIN bgc
    ON bgc.bgc_id = source.bgc_pred_id
  WHERE cciss_future12_array.siteno IN (", paste(unique(siteno), collapse = ","), ")
  AND futureperiod IN ('2001', '2021','2041','2061','2081')
  
  ), cciss_count_den AS (
  
    SELECT ", groupby, " siteref,
           futureperiod,
           SUM(weight) w
    FROM cciss
    GROUP BY ", groupby, ", futureperiod
  
  ), cciss_count_num AS (
  
    SELECT ", groupby, " siteref,
           futureperiod,
           bgc,
           bgc_pred,
           SUM(weight) w
    FROM cciss
    GROUP BY ", groupby, ", futureperiod, bgc, bgc_pred
  
  ), cciss_curr AS (
      SELECT cciss_prob12.siteno,
      period,
      bgc_attribution.bgc,
      bgc_pred,
      prob
      FROM cciss_prob12
      JOIN bgc_attribution
      ON (cciss_prob12.siteno = bgc_attribution.siteno)
      WHERE cciss_prob12.siteno IN (", paste(unique(siteno), collapse = ","), ")
      
  ), curr_temp AS (
    SELECT ", groupby, " siteref,
           COUNT(distinct siteno) n
    FROM cciss_curr
    GROUP BY ", groupby, "
  )
  
  SELECT cast(a.siteref as text) siteref,
         a.futureperiod,
         a.bgc,
         a.bgc_pred,
         a.w/cast(b.w as float) bgc_prop
  FROM cciss_count_num a
  JOIN cciss_count_den b
    ON a.siteref = b.siteref
   AND a.futureperiod = b.futureperiod
   WHERE a.w <> 0
  
  UNION ALL

  SELECT cast(", groupby, " as text) siteref,
          period as futureperiod,
          bgc,
          bgc_pred,
          SUM(prob)/b.n bgc_prop
  FROM cciss_curr a
  JOIN curr_temp b
    ON a.",groupby," = b.siteref
  WHERE siteno in (", paste(unique(siteno), collapse = ","), ")
  GROUP BY ", groupby, ",period,b.n, bgc, bgc_pred
  
  UNION ALL

  SELECT DISTINCT 
            cast(", groupby, " as text) siteref,
            '1961' as futureperiod,
            bgc,
            bgc as bgc_pred,
            cast(1 as numeric) bgc_prop
    FROM cciss_curr
    WHERE siteno IN (", paste(unique(siteno), collapse = ","), ")
  ")
  
  dat <- setDT(RPostgres::dbGetQuery(con, cciss_sql))

  setnames(dat, c("SiteRef","FuturePeriod","BGC","BGC.pred","BGC.prop"))
  dat <- unique(dat) ##should fix database so not necessary
  #print(dat)
  return(dat)
}

#' Pull CCISS version 13 from a vector of SiteNo
#' @param con An active postgres DBI connection.
#' @param siteno A character vector of siteno.
#' @param avg A boolean. 
#' @param modWeights A data table of gcm and rcp weights.
#' @details Get CCISS for provided SiteNo.
#' @return A data.table containing CCISS information for each provided SiteNo.
#' @importFrom RPostgres dbGetQuery
#' @export
dbGetCCISS_v13 <- function(con, siteno, avg, modWeights){
  
  # Declare binding for checks
  if (FALSE) {
    comb <- gcm <- rcp <- weight <- NULL
  }
  
  groupby = "siteno"
  if (isTRUE(avg)) {
    groupby = "bgc"
  }
  modWeights[,comb := paste0("('",gcm,"','",rcp,"',",weight,")")]
  weights <- paste(modWeights$comb,collapse = ",")
  
  cciss_sql <- paste0("
  WITH cciss AS (
    SELECT cciss_future_array.siteno,
         labels.gcm,
         labels.scenario,
         labels.futureperiod,
         labels.run,
         bgc_attribution13_1.bgc,
         bgcv13_1.bgc bgc_pred,
         w.weight
  FROM cciss_future_array
  JOIN bgc_attribution13_1
    ON (cciss_future_array.siteno = bgc_attribution13_1.siteno),
       unnest(bgc_id) WITH ordinality as source(bgc_id, row_idx)
  JOIN (SELECT ROW_NUMBER() OVER(ORDER BY gcm_id, scenario_id, futureperiod_id, run_id) row_idx,
               gcm,
               scenario,
               futureperiod,
               run
        FROM gcm 
        CROSS JOIN scenario
        CROSS JOIN futureperiod
        CROSS JOIN run) labels
    ON labels.row_idx = source.row_idx
    JOIN (values ",weights,") 
    AS w(gcm,scenario,weight)
    ON labels.gcm = w.gcm AND labels.scenario = w.scenario
  JOIN bgcv13_1
    ON bgcv13_1.bgc_id = source.bgc_id
  WHERE cciss_future_array.siteno IN (", paste(unique(siteno), collapse = ","), ")
  
  ), cciss_count_den AS (
  
    SELECT ", groupby, " siteref,
           futureperiod,
           SUM(weight) w
    FROM cciss
    GROUP BY ", groupby, ", futureperiod
  
  ), cciss_count_num AS (
  
    SELECT ", groupby, " siteref,
           futureperiod,
           bgc,
           bgc_pred,
           SUM(weight) w
    FROM cciss
    GROUP BY ", groupby, ", futureperiod, bgc, bgc_pred
  
  ), cciss_curr AS (
      SELECT cciss_current.siteno,
      '1991' as period,
      bgc_attribution13_1.bgc,
      bgc_pred,
      cast (1 as numeric) prob
      FROM cciss_current
      JOIN bgc_attribution13_1
      ON (cciss_current.siteno = bgc_attribution13_1.siteno)
      WHERE cciss_current.siteno IN (", paste(unique(siteno), collapse = ","), ")
      
  ), curr_temp AS (
    SELECT ", groupby, " siteref,
           COUNT(distinct siteno) n
    FROM cciss_curr
    GROUP BY ", groupby, "
  )
  
  SELECT cast(a.siteref as text) siteref,
         a.futureperiod,
         a.bgc,
         a.bgc_pred,
         a.w/cast(b.w as float) bgc_prop
  FROM cciss_count_num a
  JOIN cciss_count_den b
    ON a.siteref = b.siteref
   AND a.futureperiod = b.futureperiod
   WHERE a.w <> 0
  
  UNION ALL

  SELECT cast(", groupby, " as text) siteref,
          period as futureperiod,
          bgc,
          bgc_pred,
          SUM(prob)/b.n bgc_prop
  FROM cciss_curr a
  JOIN curr_temp b
    ON a.",groupby," = b.siteref
  WHERE siteno in (", paste(unique(siteno), collapse = ","), ")
  GROUP BY ", groupby, ",period,b.n, bgc, bgc_pred
  
  UNION ALL

  SELECT DISTINCT 
            cast(", groupby, " as text) siteref,
            '1961' as futureperiod,
            bgc,
            bgc as bgc_pred,
            cast(1 as numeric) bgc_prop
    FROM cciss_curr
    WHERE siteno IN (", paste(unique(siteno), collapse = ","), ")
  ")
  
  dat <- setDT(RPostgres::dbGetQuery(con, cciss_sql))
  
  setnames(dat, c("SiteRef","FuturePeriod","BGC","BGC.pred","BGC.prop"))
  #dat <- unique(dat) ##should fix database so not necessary
  #print(dat)
  return(dat)
}


#' Pull CCISS version 13 from a vector of SiteNo with novelty 
#' @param con An active postgres DBI connection.
#' @param siteno A character vector of siteno.
#' @param avg A boolean. 
#' @param modWeights A data table of gcm and rcp weights.
#' @param nov_cutoff Numeric between 1 and 8. Sigma value specifying novelty cutoff. 
#' @details Get CCISS for provided SiteNo.
#' @return A data.table containing CCISS information for each provided SiteNo.
#' @importFrom RPostgres dbGetQuery
#' @export
#' 
dbGetCCISS_novelty <- function(con, siteno, avg, modWeights, nov_cutoff = 5){
  
  # Declare binding for checks
  if (FALSE) {
    comb <- gcm <- rcp <- weight <- NULL
  }
  
  groupby = "siteno"
  if (isTRUE(avg)) {
    groupby = "bgc"
  }
  modWeights[,comb := paste0("('",gcm,"','",rcp,"',",weight,")")]
  weights <- paste(modWeights$comb,collapse = ",")
  
  nc <- as.character(round(nov_cutoff * 10,2))
  
  cciss_sql <- paste0("
  
  WITH 
  cciss_nov AS (
    SELECT cciss_novelty_array.siteno,
    source.novelty,
    source.row_idx
    FROM cciss_novelty_array,
    unnest(novelty) WITH ordinality as source(novelty, row_idx)
    WHERE cciss_novelty_array.siteno IN (", paste(unique(siteno), collapse = ","), ")
  ),
  
    cciss_bgc AS (
    SELECT cciss_future_array.siteno,
    source.row_idx,
    labels.gcm,
    labels.scenario,
    labels.futureperiod,
    labels.run,
    bgc_attribution13_1.bgc,
    bgcv13_1.bgc bgc_pred,
    w.weight
    FROM cciss_future_array
    JOIN bgc_attribution13_1
    ON (cciss_future_array.siteno = bgc_attribution13_1.siteno),
    unnest(bgc_id) WITH ordinality as source(bgc_id, row_idx)
    JOIN (SELECT ROW_NUMBER() OVER(ORDER BY gcm_id, scenario_id, futureperiod_id, run_id) row_idx,
          gcm,
          scenario,
          futureperiod,
          run
          FROM gcm 
          CROSS JOIN scenario
          CROSS JOIN futureperiod
          CROSS JOIN run) labels
    ON labels.row_idx = source.row_idx
    JOIN (values ",weights,") 
    AS w(gcm,scenario,weight)
    ON labels.gcm = w.gcm AND labels.scenario = w.scenario
    JOIN bgcv13_1
    ON bgcv13_1.bgc_id = source.bgc_id
    WHERE cciss_future_array.siteno IN (", paste(unique(siteno), collapse = ","), ")
    
  ),
  
  cciss AS (
  SELECT cciss_bgc.siteno,
  gcm,
  scenario, 
  futureperiod,
  run,
  bgc,
  CASE WHEN cciss_nov.novelty > ",nc," THEN 'novel' ELSE bgc_pred END AS bgc_pred,
  cciss_nov.novelty,
  weight
  FROM cciss_bgc
  JOIN cciss_nov USING (siteno, row_idx)
  ),
  
  cciss_count_den AS (
    
    SELECT ", groupby, " siteref,
    futureperiod,
    SUM(weight) w
    FROM cciss
    GROUP BY ", groupby, ", futureperiod
    
  ), cciss_count_num AS (
    
    SELECT ", groupby, " siteref,
    futureperiod,
    bgc,
    bgc_pred,
    AVG(novelty) nov,
    SUM(weight) w
    FROM cciss
    GROUP BY ", groupby, ", futureperiod, bgc, bgc_pred
    
  )  ,
  
  cciss_curr AS (
      SELECT cciss_current_nov.siteno,
      '1991' as period,
      bgc_attribution13_1.bgc,
      CASE WHEN novelty > ",nc," THEN 'novel' ELSE bgc_pred END AS bgc_pred,
      cast (1 as numeric) prob,
      novelty
      FROM cciss_current_nov
      JOIN bgc_attribution13_1
      ON (cciss_current_nov.siteno = bgc_attribution13_1.siteno)
      WHERE cciss_current_nov.siteno IN (", paste(unique(siteno), collapse = ","), ")
      
  ), curr_temp AS (
    SELECT ", groupby, " siteref,
           COUNT(distinct siteno) n
    FROM cciss_curr
    GROUP BY ", groupby, "
  )
  
  SELECT cast(a.siteref as text) siteref,
         a.futureperiod,
         a.bgc,
         a.bgc_pred,
         a.w/cast(b.w as float) bgc_prop,
         a.nov novelty
  FROM cciss_count_num a
  JOIN cciss_count_den b
    ON a.siteref = b.siteref
   AND a.futureperiod = b.futureperiod
   WHERE a.w <> 0
  
  UNION ALL

  SELECT cast(", groupby, " as text) siteref,
          period as futureperiod,
          bgc,
          bgc_pred,
          SUM(prob)/b.n bgc_prop,
          AVG(novelty) novelty
  FROM cciss_curr a
  JOIN curr_temp b
    ON a.",groupby," = b.siteref
  WHERE siteno in (", paste(unique(siteno), collapse = ","), ")
  GROUP BY ", groupby, ",period,b.n, bgc, bgc_pred
  
  UNION ALL

  SELECT DISTINCT 
            cast(", groupby, " as text) siteref,
            '1961' as futureperiod,
            bgc,
            bgc as bgc_pred,
            cast(1 as numeric) bgc_prop,
            cast(0 as numeric) novelty
    FROM cciss_curr
    WHERE siteno IN (", paste(unique(siteno), collapse = ","), ")")
  
  dat <- setDT(RPostgres::dbGetQuery(con, cciss_sql))
  
  setnames(dat, c("SiteRef","FuturePeriod","BGC","BGC.pred","BGC.prop","Novelty"))
  return(dat)
}

#' Pull individual predictions by district
#' @param con An active postgres DBI connection.
#' @param siteno A character vector of siteno.
#' @param avg A boolean. 
#' @param modWeights A data table of gcm and rcp weights.
#' @details Get CCISS for provided SiteNo.
#' @return A data.table containing CCISS information for each provided SiteNo.
#' @importFrom RPostgres dbGetQuery
#' @export
dbGetCCISSRaw <- function(con, district, gcm, scenario, period){
  
  cciss_sql <- paste0("
    SELECT cciss_future12_array.siteno,
         labels.gcm,
         labels.scenario,
         labels.futureperiod,
         bgc.bgc bgc_pred
  FROM cciss_future12_array
  JOIN (select dist_code, siteno from bgc_dist_ids where dist_code = '",district, "') dists
    ON (dists.siteno = cciss_future12_array.siteno),
  unnest(bgc_pred_id) WITH ordinality as source(bgc_pred_id, row_idx)
  JOIN (SELECT ROW_NUMBER() OVER(ORDER BY gcm_id, scenario_id, futureperiod_id) row_idx,
               gcm,
               scenario,
               futureperiod
        FROM gcm 
        CROSS JOIN scenario
        CROSS JOIN futureperiod) labels
    ON labels.row_idx = source.row_idx
  JOIN bgc
    ON bgc.bgc_id = source.bgc_pred_id
  WHERE futureperiod = '",period,"'
  AND scenario = '",scenario,"'
  AND gcm = '",gcm,"'
  ")
  
  dat <- setDT(RPostgres::dbGetQuery(con, cciss_sql))
  dat <- unique(dat) ##should fix database so not necessary
  #print(dat)
  return(dat)
}

#' Return raw predictions for given sitenos
#' @param con An active postgres DBI connection.
#' @param siteno A vector of siteno.
#' @param scenario Vector of requested SSPs
#' @param period Character vector of requested future periods.
#' @details Return data from CCISS database for provided SiteNo.
#' @return A data.table containing BGC predictions.
#' @importFrom RPostgres dbGetQuery
#' @export
dbGetBGCPred <- function(con, siteno, scenario = c("ssp126","ssp245","ssp370"), period = c('2001', '2021','2041','2061','2081')){
  
  cciss_sql <- paste0("
    SELECT cciss_future12_array.siteno,
         labels.gcm,
         labels.scenario,
         labels.futureperiod,
         bgc_attribution.bgc bgc,
         bgc.bgc bgc_pred
  FROM cciss_future12_array
  JOIN bgc_attribution
    ON (cciss_future12_array.siteno = bgc_attribution.siteno),
    unnest(bgc_pred_id) WITH ordinality as source(bgc_pred_id, row_idx)

  JOIN (SELECT ROW_NUMBER() OVER(ORDER BY gcm_id, scenario_id, futureperiod_id) row_idx,
               gcm,
               scenario,
               futureperiod
        FROM gcm 
        CROSS JOIN scenario
        CROSS JOIN futureperiod) labels
    ON labels.row_idx = source.row_idx
  JOIN bgc
    ON bgc.bgc_id = source.bgc_pred_id
    WHERE cciss_future12_array.siteno IN(",paste(unique(siteno), collapse = ","),")
  AND futureperiod IN('",paste(period, collapse = "','"),"')
  AND scenario IN('",paste(scenario, collapse = "','"),"')
  ")
  
  dat <- setDT(RPostgres::dbGetQuery(con, cciss_sql))
  dat <- unique(dat) ##should fix database so not necessary
  #print(dat)
  return(dat)
}

