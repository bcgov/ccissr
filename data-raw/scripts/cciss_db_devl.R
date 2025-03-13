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


t1 <- dbGetQuery(conn, "select * from cciss_future_array limit 5")
t2 <- dbGetQuery(conn, "select * from cciss_novelty_array where siteno = 49")


### create parameter input charts
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
all_weight[,comb := paste0("('",gcm,"','",rcp,"',",weight,")")]
weights <- paste(all_weight$comb,collapse = ",")

siteno <- c(590351)
groupby = "siteno"

cciss_sql <- paste0("
  
  WITH 
   cciss_nov AS (
    SELECT cciss_novelty_array.siteno,
         labels.gcm,
         labels.scenario,
         labels.futureperiod,
         labels.run,
         source.novelty,
         source.row_idx
  FROM cciss_novelty_array,
       unnest(novelty) WITH ordinality as source(novelty, row_idx)
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
    WHERE cciss_novelty_array.siteno IN (", paste(unique(siteno), collapse = ","), ")
  ),
  
  cciss AS (
    SELECT cciss_future_array.siteno,
         labels.gcm,
         labels.scenario,
         labels.futureperiod,
         labels.run,
         bgc_attribution13_1.bgc,
         CASE WHEN cciss_nov.novelty > 50 THEN 'novel' ELSE bgcv13_1.bgc END AS bgc_pred,
         cciss_nov.novelty,
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
  JOIN cciss_nov
    ON cciss_nov.row_idx = source.row_idx
  WHERE cciss_future_array.siteno IN (", paste(unique(siteno), collapse = ","), ")
  
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
  
  ) ,
  
  cciss_curr AS (
      SELECT cciss_current_nov.siteno,
      '1991' as period,
      bgc_attribution13_1.bgc,
      CASE WHEN novelty > 50 THEN 'novel' ELSE bgc_pred END AS bgc_pred,
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
    WHERE siteno IN (", paste(unique(siteno), collapse = ","), ")
  
                    ")

test <- dbGetQuery(conn, cciss_sql)
setDT(test)
t2 <- test[,.(bgc_prop = sum(weight), nov_mean = mean(novelty), nov_sd = sd(novelty)), by = .(siteno, futureperiod, bgc, bgc_pred)]

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
    SELECT cciss_future13_array.siteno,
         labels.gcm,
         labels.scenario,
         labels.futureperiod,
         labels.run,
         bgc_attribution13.bgc,
         bgcv13.bgc bgc_pred,
         w.weight
  FROM cciss_future13_array
  JOIN bgc_attribution13
    ON (cciss_future13_array.siteno = bgc_attribution13.siteno),
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
  JOIN bgcv13
    ON bgcv13.bgc_id = source.bgc_id
  WHERE cciss_future13_array.siteno IN (", paste(unique(siteno), collapse = ","), ")
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
      SELECT cciss_prob13.siteno,
      '1991' as period,
      bgc_attribution13.bgc,
      bgc_pred,
      prob
      FROM cciss_prob13
      JOIN bgc_attribution13
      ON (cciss_prob13.siteno = bgc_attribution13.siteno)
      WHERE cciss_prob13.siteno IN (", paste(unique(siteno), collapse = ","), ")
      
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
