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
dbGetBGC <- function(con,bgc,district = NULL,maxPoints = 150){
  if(is.null(district)){
    query <- paste0("select siteno from preselected_points14 where bgc IN ('",paste(bgc,collapse = "','"),"') limit ", maxPoints)
  }else{
    query <- paste0("select siteno from preselected_dist13 where bgc IN ('",paste(bgc,collapse = "','"),"') and dist_code = '",district,"'")
  }
  dat <- RPostgres::dbGetQuery(con, query)
  if(nrow(dat) == 0){ #temporary check until preselected_dist table is updated
    query <- paste0("select siteno from preselected_points14 where bgc IN ('",paste(bgc,collapse = "','"),"') limit ", maxPoints)
    dat <- RPostgres::dbGetQuery(con, query)
  } 
  return(dat$siteno)
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
#' @importFrom glue glue glue_sql
#' @export
dbGetCCISS_v13 <- function(con, siteno, avg, modWeights, 
                           cciss_table = "cciss_novelty_array", 
                           cciss_observed = "cciss_current",
                           bgc_table = "bgc_attribution", 
                           bgc_lookup = "bgc_v13") {
  
  if (FALSE) {
    comb <- gcm <- rcp <- weight <- NULL
  }
  
  groupby <- if (isTRUE(avg)) "bgc" else "siteno"
  
  modWeights[,comb := paste0("('",gcm,"','",rcp,"',",weight,")")]
  weights <- DBI::SQL(paste(modWeights$comb, collapse = ","))
  siteno_sql <- unique(siteno)
  groupby_sql <- DBI::SQL(groupby)
  cciss_table <- DBI::SQL(cciss_table)
  cciss_observed <- DBI::SQL(cciss_observed)
  bgc_table <- DBI::SQL(bgc_table)
  bgc_lookup <- DBI::SQL(bgc_lookup)
  
  cciss_sql <- glue::glue_sql("
    WITH cciss AS (
      SELECT {cciss_table}.siteno,
             labels.gcm,
             labels.scenario,
             labels.futureperiod,
             labels.run,
             {bgc_table}.bgc,
             {bgc_lookup}.bgc bgc_pred,
             w.weight
      FROM {cciss_table}
      JOIN {bgc_table}
        ON {cciss_table}.siteno = {bgc_table}.siteno,
           unnest(bgc_id) WITH ordinality AS source(bgc_id, row_idx)
      JOIN (
        SELECT ROW_NUMBER() OVER(ORDER BY gcm_id, scenario_id, futureperiod_id, run_id) row_idx,
               gcm,
               scenario,
               futureperiod,
               run
        FROM gcm 
        CROSS JOIN scenario
        CROSS JOIN futureperiod
        CROSS JOIN run
      ) labels
        ON labels.row_idx = source.row_idx
      JOIN (VALUES {weights}) AS w(gcm, scenario, weight)
        ON labels.gcm = w.gcm 
       AND labels.scenario = w.scenario
      JOIN {bgc_lookup}
        ON {bgc_lookup}.bgc_id = source.bgc_id
      WHERE {cciss_table}.siteno IN ({siteno_sql*})
    
    ), cciss_count_den AS (
      SELECT {groupby_sql} siteref,
             futureperiod,
             SUM(weight) w
      FROM cciss
      GROUP BY {groupby_sql}, futureperiod
    
    ), cciss_count_num AS (
      SELECT {groupby_sql} siteref,
             futureperiod,
             bgc,
             bgc_pred,
             SUM(weight) w
      FROM cciss
      GROUP BY {groupby_sql}, futureperiod, bgc, bgc_pred
    
    ), cciss_curr AS (
      SELECT {cciss_observed}.siteno,
             '1991' AS period,
             {bgc_table}.bgc,
             bgc_pred,
             CAST(1 AS numeric) prob
      FROM {cciss_observed}
      JOIN {bgc_table}
        ON {cciss_observed}.siteno = {bgc_table}.siteno
      WHERE {cciss_observed}.siteno IN ({siteno_sql*})
      
    ), curr_temp AS (
      SELECT {groupby_sql} siteref,
             COUNT(DISTINCT siteno) n
      FROM cciss_curr
      GROUP BY {groupby_sql}
    )
    
    SELECT CAST(a.siteref AS text) siteref,
           a.futureperiod,
           a.bgc,
           a.bgc_pred,
           a.w / CAST(b.w AS float) bgc_prop
    FROM cciss_count_num a
    JOIN cciss_count_den b
      ON a.siteref = b.siteref
     AND a.futureperiod = b.futureperiod
    WHERE a.w <> 0
    
    UNION ALL

    SELECT CAST({groupby_sql} AS text) siteref,
           period AS futureperiod,
           bgc,
           bgc_pred,
           SUM(prob) / b.n bgc_prop
    FROM cciss_curr a
    JOIN curr_temp b
      ON a.{groupby_sql} = b.siteref
    WHERE siteno IN ({siteno_sql*})
    GROUP BY {groupby_sql}, period, b.n, bgc, bgc_pred
    
    UNION ALL

    SELECT DISTINCT 
           CAST({groupby_sql} AS text) siteref,
           '1961' AS futureperiod,
           bgc,
           bgc AS bgc_pred,
           CAST(1 AS numeric) bgc_prop
    FROM cciss_curr
  ", .con = con)
  
  dat <- data.table::setDT(RPostgres::dbGetQuery(con, cciss_sql))
  
  data.table::setnames(dat, c("SiteRef", "FuturePeriod", "BGC", "BGC.pred", "BGC.prop"))
  
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
#' @importFrom glue glue glue_sql
#' @export
#' 
dbGetCCISS_novelty <- function(con, siteno, avg, modWeights, nov_cutoff = 5,
                               cciss_table = "cciss_future_array",
                               novelty_table = "cciss_novelty_array",
                               cciss_observed = "cciss_current_nov",
                               bgc_table = "bgc_attribution",
                               bgc_lookup = "bgc_v13") {
  
  if (FALSE) {
    comb <- gcm <- rcp <- weight <- NULL
  }
  
  groupby <- if (isTRUE(avg)) "bgc" else "siteno"
  
  modWeights[,comb := paste0("('",gcm,"','",rcp,"',",weight,")")]
  weights <- DBI::SQL(paste(modWeights$comb, collapse = ","))
  siteno_sql <- unique(siteno)
  groupby_sql <- DBI::SQL(groupby)
  
  cciss_table <- DBI::SQL(cciss_table)
  novelty_table <- DBI::SQL(novelty_table)
  cciss_observed <- DBI::SQL(cciss_observed)
  bgc_table <- DBI::SQL(bgc_table)
  bgc_lookup <- DBI::SQL(bgc_lookup)
  
  nc <- round(nov_cutoff * 10, 2)
  
  cciss_sql <- glue::glue_sql("
    WITH labels AS MATERIALIZED (
  SELECT
    ROW_NUMBER() OVER (
      ORDER BY gcm_id, scenario_id, futureperiod_id, run_id
    ) AS row_idx,
    gcm,
    scenario,
    futureperiod,
    run
  FROM gcm
  CROSS JOIN scenario
  CROSS JOIN futureperiod
  CROSS JOIN run
),

weights AS MATERIALIZED (
  SELECT *
  FROM (VALUES {weights}) AS v(gcm, scenario, weight)
),

weighted_labels AS MATERIALIZED (
  SELECT
    l.row_idx,
    l.gcm,
    l.scenario,
    l.futureperiod,
    l.run,
    w.weight
  FROM labels AS l
  JOIN weights AS w
    ON w.gcm = l.gcm
   AND w.scenario = l.scenario
),

selected_bgc AS MATERIALIZED (
  SELECT
    c.siteno,
    x.bgc_id,
    x.row_idx
  FROM {cciss_table} AS c
  CROSS JOIN LATERAL unnest(c.bgc_id)
    WITH ORDINALITY AS x(bgc_id, row_idx)
  JOIN weighted_labels AS wl
    ON wl.row_idx = x.row_idx
  WHERE c.siteno IN ({siteno_sql*})
),

selected_novelty AS MATERIALIZED (
  SELECT
    n.siteno,
    x.novelty,
    x.row_idx
  FROM {novelty_table} AS n
  CROSS JOIN LATERAL unnest(n.novelty)
    WITH ORDINALITY AS x(novelty, row_idx)
  JOIN weighted_labels AS wl
    ON wl.row_idx = x.row_idx
  WHERE n.siteno IN ({siteno_sql*})
),

cciss AS MATERIALIZED (
  SELECT
    b.siteno,
    wl.gcm,
    wl.scenario,
    wl.futureperiod,
    wl.run,
    ba.bgc,
    CASE
      WHEN n.novelty > {nc} THEN 'novel'
      ELSE bl.bgc
    END AS bgc_pred,
    n.novelty,
    wl.weight
  FROM selected_bgc AS b

  JOIN selected_novelty AS n
    USING (siteno, row_idx)

  JOIN weighted_labels AS wl
    USING (row_idx)

  JOIN {bgc_table} AS ba
    USING (siteno)

  JOIN {bgc_lookup} AS bl
    ON bl.bgc_id = b.bgc_id
),
    
    cciss_count_den AS (
      SELECT {groupby_sql} siteref,
             futureperiod,
             SUM(weight) w
      FROM cciss
      GROUP BY {groupby_sql}, futureperiod
    ),
    
    cciss_count_num AS (
      SELECT {groupby_sql} siteref,
             futureperiod,
             bgc,
             bgc_pred,
             AVG(novelty) nov,
             SUM(weight) w
      FROM cciss
      GROUP BY {groupby_sql}, futureperiod, bgc, bgc_pred
    ),
    
    cciss_curr AS (
      SELECT {cciss_observed}.siteno,
             '1991' AS period,
             {bgc_table}.bgc,
             bgc_pred,
             CAST(1 AS numeric) prob,
             CAST(0 as numeric) novelty -- need to make a novelty table
      FROM {cciss_observed}
      JOIN {bgc_table}
        ON {cciss_observed}.siteno = {bgc_table}.siteno
      WHERE {cciss_observed}.siteno IN ({siteno_sql*})
    ),
    
    curr_temp AS (
      SELECT {groupby_sql} siteref,
             COUNT(DISTINCT siteno) n
      FROM cciss_curr
      GROUP BY {groupby_sql}
    )
    
    SELECT CAST(a.siteref AS text) siteref,
           a.futureperiod,
           a.bgc,
           a.bgc_pred,
           a.w / CAST(b.w AS float) bgc_prop,
           a.nov novelty
    FROM cciss_count_num a
    JOIN cciss_count_den b
      ON a.siteref = b.siteref
     AND a.futureperiod = b.futureperiod
    WHERE a.w <> 0
    
    UNION ALL

    SELECT CAST({groupby_sql} AS text) siteref,
           period AS futureperiod,
           bgc,
           bgc_pred,
           SUM(prob) / b.n bgc_prop,
           AVG(novelty) novelty
    FROM cciss_curr a
    JOIN curr_temp b
      ON a.{groupby_sql} = b.siteref
    WHERE siteno IN ({siteno_sql*})
    GROUP BY {groupby_sql}, period, b.n, bgc, bgc_pred
    
    UNION ALL

    SELECT DISTINCT 
           CAST({groupby_sql} AS text) siteref,
           '1961' AS futureperiod,
           bgc,
           bgc AS bgc_pred,
           CAST(1 AS numeric) bgc_prop,
           CAST(0 AS numeric) novelty
    FROM cciss_curr
  ", .con = con)
  
  
  dat <- data.table::setDT(RPostgres::dbGetQuery(con, cciss_sql))
  
  data.table::setnames(
    dat,
    c("SiteRef", "FuturePeriod", "BGC", "BGC.pred", "BGC.prop", "Novelty")
  )
  
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

