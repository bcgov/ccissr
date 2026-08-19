library(RPostgres)
library(DBI)
library(data.table)
library(sf)
library(terra)
library(climr)
library(ranger)
library(ccissr)
devtools::load_all()

con <- DBI::dbConnect(
  drv = RPostgres::Postgres(),
  dbname = "cciss",
  host = Sys.getenv("BCGOV_HOST"),
  port = 5432, 
  user = Sys.getenv("BCGOV_USR"),
  password = Sys.getenv("BCGOV_PWD")
)

siteno <- dbGetBGC(con, bgc = "SBSdk", maxPoints = 150)

cciss_table = "cciss_future14_array"
novelty_table = "cciss_novelty14_array"
cciss_observed = "cciss_current14"
bgc_table = "bgc_attribution14"
bgc_lookup = "bgc14"
avg <- TRUE
modWeights <- all_weight
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

nc <- round(5 * 10, 2)

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

dat <- dbGetQuery(con, cciss_sql)

explain_sql <- paste(
  "EXPLAIN (ANALYZE, BUFFERS, VERBOSE, SETTINGS, FORMAT TEXT)",
  as.character(cciss_sql)
)

plan <- DBI::dbGetQuery(con, explain_sql)
cat(paste(plan[[1]], collapse = "\n"))

res1 <- dbGetCCISS_novelty(con, siteno, avg = TRUE, modWeights = all_weight, cciss_table = "cciss_future14_array", 
                                  novelty_table = "cciss_novelty14_array",
                                  cciss_observed = "cciss_current14",
                                  bgc_table = "bgc_attribution14", 
                                  bgc_lookup = "bgc14", nov_cutoff = 5)

res2 <- dbGetCCISS_old(con, siteno, avg = TRUE, modWeights = all_weight)

dbGetCCISS_old <- function (con, siteno, avg, modWeights, nov_cutoff = 5) 
{
  if (FALSE) {
    comb <- gcm <- rcp <- weight <- NULL
  }
  groupby = "siteno"
  if (isTRUE(avg)) {
    groupby = "bgc"
  }
  modWeights[, `:=`(comb, paste0("('", gcm, "','", rcp, "',", 
                                 weight, ")"))]
  weights <- paste(modWeights$comb, collapse = ",")
  nc <- as.character(round(nov_cutoff * 10, 2))
  cciss_sql <- paste0("\n  \n  WITH \n  cciss_nov AS (\n    SELECT cciss_novelty_array.siteno,\n    source.novelty,\n    source.row_idx\n    FROM cciss_novelty_array,\n    unnest(novelty) WITH ordinality as source(novelty, row_idx)\n    WHERE cciss_novelty_array.siteno IN (", 
                      paste(unique(siteno), collapse = ","), ")\n  ),\n  \n    cciss_bgc AS (\n    SELECT cciss_future_array.siteno,\n    source.row_idx,\n    labels.gcm,\n    labels.scenario,\n    labels.futureperiod,\n    labels.run,\n    bgc_attribution13_1.bgc,\n    bgcv13_1.bgc bgc_pred,\n    w.weight\n    FROM cciss_future_array\n    JOIN bgc_attribution13_1\n    ON (cciss_future_array.siteno = bgc_attribution13_1.siteno),\n    unnest(bgc_id) WITH ordinality as source(bgc_id, row_idx)\n    JOIN (SELECT ROW_NUMBER() OVER(ORDER BY gcm_id, scenario_id, futureperiod_id, run_id) row_idx,\n          gcm,\n          scenario,\n          futureperiod,\n          run\n          FROM gcm \n          CROSS JOIN scenario\n          CROSS JOIN futureperiod\n          CROSS JOIN run) labels\n    ON labels.row_idx = source.row_idx\n    JOIN (values ", 
                      weights, ") \n    AS w(gcm,scenario,weight)\n    ON labels.gcm = w.gcm AND labels.scenario = w.scenario\n    JOIN bgcv13_1\n    ON bgcv13_1.bgc_id = source.bgc_id\n    WHERE cciss_future_array.siteno IN (", 
                      paste(unique(siteno), collapse = ","), ")\n    \n  ),\n  \n  cciss AS (\n  SELECT cciss_bgc.siteno,\n  gcm,\n  scenario, \n  futureperiod,\n  run,\n  bgc,\n  CASE WHEN cciss_nov.novelty > ", 
                      nc, " THEN 'novel' ELSE bgc_pred END AS bgc_pred,\n  cciss_nov.novelty,\n  weight\n  FROM cciss_bgc\n  JOIN cciss_nov USING (siteno, row_idx)\n  ),\n  \n  cciss_count_den AS (\n    \n    SELECT ", 
                      groupby, " siteref,\n    futureperiod,\n    SUM(weight) w\n    FROM cciss\n    GROUP BY ", 
                      groupby, ", futureperiod\n    \n  ), cciss_count_num AS (\n    \n    SELECT ", 
                      groupby, " siteref,\n    futureperiod,\n    bgc,\n    bgc_pred,\n    AVG(novelty) nov,\n    SUM(weight) w\n    FROM cciss\n    GROUP BY ", 
                      groupby, ", futureperiod, bgc, bgc_pred\n    \n  )  ,\n  \n  cciss_curr AS (\n      SELECT cciss_current_nov.siteno,\n      '1991' as period,\n      bgc_attribution13_1.bgc,\n      CASE WHEN novelty > ", 
                      nc, " THEN 'novel' ELSE bgc_pred END AS bgc_pred,\n      cast (1 as numeric) prob,\n      novelty\n      FROM cciss_current_nov\n      JOIN bgc_attribution13_1\n      ON (cciss_current_nov.siteno = bgc_attribution13_1.siteno)\n      WHERE cciss_current_nov.siteno IN (", 
                      paste(unique(siteno), collapse = ","), ")\n      \n  ), curr_temp AS (\n    SELECT ", 
                      groupby, " siteref,\n           COUNT(distinct siteno) n\n    FROM cciss_curr\n    GROUP BY ", 
                      groupby, "\n  )\n  \n  SELECT cast(a.siteref as text) siteref,\n         a.futureperiod,\n         a.bgc,\n         a.bgc_pred,\n         a.w/cast(b.w as float) bgc_prop,\n         a.nov novelty\n  FROM cciss_count_num a\n  JOIN cciss_count_den b\n    ON a.siteref = b.siteref\n   AND a.futureperiod = b.futureperiod\n   WHERE a.w <> 0\n  \n  UNION ALL\n\n  SELECT cast(", 
                      groupby, " as text) siteref,\n          period as futureperiod,\n          bgc,\n          bgc_pred,\n          SUM(prob)/b.n bgc_prop,\n          AVG(novelty) novelty\n  FROM cciss_curr a\n  JOIN curr_temp b\n    ON a.", 
                      groupby, " = b.siteref\n  WHERE siteno in (", paste(unique(siteno), 
                                                                          collapse = ","), ")\n  GROUP BY ", groupby, ",period,b.n, bgc, bgc_pred\n  \n  UNION ALL\n\n  SELECT DISTINCT \n            cast(", 
                      groupby, " as text) siteref,\n            '1961' as futureperiod,\n            bgc,\n            bgc as bgc_pred,\n            cast(1 as numeric) bgc_prop,\n            cast(0 as numeric) novelty\n    FROM cciss_curr\n    WHERE siteno IN (", 
                      paste(unique(siteno), collapse = ","), ")")
  dat <- setDT(RPostgres::dbGetQuery(con, cciss_sql))
  setnames(dat, c("SiteRef", "FuturePeriod", "BGC", "BGC.pred", 
                  "BGC.prop", "Novelty"))
  return(dat)
}

qry <- "CREATE TABLE preselected_dist13 AS
SELECT * FROM (
  SELECT siteno, bgc, dist_code,
         ROW_NUMBER() OVER (PARTITION BY bgc, dist_code ORDER BY RANDOM()) AS u
  FROM (
    SELECT dist_code, bgc, bgc_attribution13_1.siteno
    FROM bgc_attribution13_1
    JOIN district_ids USING (siteno)
  ) AS temp
) AS a
WHERE u <= 150;
"
dbExecute(conn, qry)


create_preselected <- 
"create table preselected_points13 as (select * from (
  select siteno, bgc, row_number() over (partition by bgc order by random()) as u
  from bgc_attribution13_1
) as a
where u <= 150);"

dat <- parse_qml("../../../Downloads/WNAv13_v6_Subzones.qml")

source("./data-raw/scripts/functions.R")



bc_bgc <- st_read("../Common_Files/BEC13Draft_4326.gpkg")
bc_bgc <- bc_bgc["BGC"]
plot(bc_bgc[1,])
st_write(bc_bgc, conn, "bgc_map_v13")

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

siteno <- c(1963369,4310326)
groupby = "siteno"
nc <- "50"

cciss_sql <- paste0(
  "WITH 
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
    
  ) 
  
  select * from cciss_count_num
  
  "
  )

test <- dbGetQuery(conn, cciss_sql)
setDT(test)
t2 <- test[,.(bgc_prop = sum(weight), nov_mean = mean(novelty), nov_sd = sd(novelty)), by = .(siteno, futureperiod, bgc, bgc_pred)]



"cciss_curr AS (
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
  WHERE siteno IN (", paste(unique(siteno), collapse = ","), ")"



qry <- "SELECT cciss_future_array.siteno,
         labels.gcm,
         labels.scenario,
         labels.futureperiod,
         labels.run,
         bgc_attribution13_1.bgc,
         bgcv13_1.bgc bgc_pred
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
  JOIN bgcv13_1
    ON bgcv13_1.bgc_id = source.bgc_id
  WHERE cciss_future_array.siteno IN (1963369)"

dat <- dbGetQuery(con, qry)

nruns <- dat[,.(Num = .N), by = .(futureperiod)]

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
