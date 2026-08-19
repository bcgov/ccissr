#' Create table of species relative suitabile area
#' @param con duckdb database connection
#' @param edatope Character. Edatopic position to use. Default `C4` (zonal)
#' @param fut_wt Numeric vector of length 5, corresponding to weight of each period for future period calculation. Default `c(0,0,0,0.5,0.5,0)`, which is average of 2041-2060 and 2061-2080 periods.
#' @param curr_wt Numeric vector of length 5, corresponding to weight of each period for current period calculation. Default `c(0.5,0.5,0,0,0,0)`, which is average of refperiod and 2001-2020 periods.
#' @param periods Character vector of length 5 with names of periods. Must match names in table `cciss_res`. Default `c("Curr","2001_2020","2021_2040","2041_2060","2061_2080","2081_2100")`
#' @param BGCxDistrict Logical. Should summary be by BGC (default) or by BGC and natural resource district?
#' @return data.table containing loss and gain proportion by Spp and BGC
#' @import data.table duckdb
#' @importFrom glue glue_sql
#' @export
spp_loss_gain <- function(
    con,
    edatope = "C4",
    fut_wt  = c(0,0,0,0.5,0.5,0),
    curr_wt = c(0.5,0.5,0,0,0,0),
    periods = c("Curr","2001_2020","2021_2040","2041_2060","2061_2080","2081_2100"),
    BGCxDistrict = FALSE
) {
  stopifnot(length(fut_wt) == length(periods), length(curr_wt) == length(periods))
  
  if(BGCxDistrict) {
    if(!duckdb_table_exists(con, "dist_points")){
      stop("District table does not exist in database! Please add it using the dbAddDistricts function.")
    }
  }
  
  # create a TEMP weights table in DuckDB
  wts <- data.frame(
    period  = periods,
    fut_wt  = as.numeric(fut_wt),
    curr_wt = as.numeric(curr_wt)
  )
  dbExecute(con, "DROP TABLE IF EXISTS period_weights;")
  dbWriteTable(con, "period_weights", wts, temporary = TRUE)
  
  if(BGCxDistrict){
    bgc_qry <- glue_sql("SELECT
                f.*,
                b.bgc AS BGC,
                d.district AS District
              FROM flags f
              LEFT JOIN bgc_points b
                ON b.cellnum = f.SiteRef
              LEFT JOIN dist_points d
                ON d.cellnum = f.SiteRef",.con = con)
    group <- glue_sql("BGC, District, Spp",.con = con)
  } else {
    bgc_qry <- glue_sql("SELECT
                f.*,
                b.bgc AS BGC
              FROM flags f
              LEFT JOIN bgc_points b
                ON b.cellnum = f.SiteRef",.con = con)
    group <- glue_sql("BGC, Spp",.con = con)
  }
  
  sql <- glue_sql("
    WITH base AS (
      SELECT
        SiteRef,
        FuturePeriod,
        Spp,
        Curr,
        Newsuit
      FROM cciss_res
      WHERE Edatope = {edatope}
    ),

    curr_rows AS (
      SELECT
        SiteRef,
        Spp,
        'Curr'::VARCHAR AS period,
        MAX(Curr) AS suit
      FROM base
      GROUP BY SiteRef, Spp
    ),
    
    -- future rows 
    future_rows AS (
      SELECT
        SiteRef,
        Spp,
        FuturePeriod AS period,
        Newsuit      AS suit
      FROM base
    ),
    
    long AS (
      SELECT * FROM curr_rows
      UNION ALL
      SELECT * FROM future_rows
    ),

    scored AS (
      SELECT
        l.SiteRef,
        l.Spp,
        SUM(l.suit * w.fut_wt)  AS FutSuit,
        SUM(l.suit * w.curr_wt) AS CurrSuit
      FROM long l
      JOIN period_weights w
        ON w.period = l.period
      GROUP BY l.SiteRef, l.Spp
    ),

    flags AS (
      SELECT
        s.*,
        (s.CurrSuit < 3 AND s.FutSuit > 3.5) AS Loss,
        (s.CurrSuit > 3.5 AND s.FutSuit <= 3) AS Gain
      FROM scored s
    ),

    with_bgc AS (
        {bgc_qry}
    )

    SELECT
      {group},
      SUM(CASE WHEN Loss THEN 1 ELSE 0 END) AS LossArea,
      SUM(CASE WHEN Gain THEN 1 ELSE 0 END) AS GainArea,
      COUNT(*) as TotalArea
    FROM with_bgc
    GROUP BY {group}
    ORDER BY {group};
  ", .con = con)
  
  res <- dbGetQuery(con, sql)
  data.table::as.data.table(res)
}

#' Create table of species consistency 
#' @param con duckdb database connection
#' @param edatope Character. Edatopic position to use. Default `C4` (zonal)
#' @param fut_wt Numeric vector of length 5, corresponding to weight of each period for future period calculation. Default `c(0,0,0,0.5,0.5,0)`, which is average of 2041-2060 and 2061-2080 periods.
#' @param curr_wt Numeric vector of length 5, corresponding to weight of each period for current period calculation. Default `c(0.5,0.5,0,0,0,0)`, which is average of refperiod and 2001-2020 periods.
#' @param periods Character vector of length 5 with names of periods. Must match names in table `cciss_res`. Default `c("Curr","2001_2020","2021_2040","2041_2060","2061_2080","2081_2100")`
#' @return data.table containing loss and gain proportion by Spp and BGC
#' @import data.table duckdb
#' @importFrom glue glue_sql
#' @export
spp_persistance <- function(
    con,
    edatope = "C4",
    fut_wt  = c(0,0,0,0.5,0.5,0),
    curr_wt = c(0.5,0.5,0,0,0,0),
    periods = c("Curr","2001_2020","2021_2040","2041_2060","2061_2080","2081_2100")
) {
  stopifnot(length(fut_wt) == length(periods), length(curr_wt) == length(periods))
  
  # create a TEMP weights table in DuckDB
  wts <- data.frame(
    period  = periods,
    fut_wt  = as.numeric(fut_wt),
    curr_wt = as.numeric(curr_wt)
  )
  dbExecute(con, "DROP TABLE IF EXISTS period_weights;")
  dbWriteTable(con, "period_weights", wts, temporary = TRUE)
  
  sql <- glue_sql("
    WITH base AS (
      SELECT
        SiteRef,
        FuturePeriod,
        Spp,
        Curr,
        Newsuit
      FROM cciss_res
      WHERE Edatope = {edatope}
    ),

    curr_rows AS (
      SELECT
        SiteRef,
        Spp,
        'Curr'::VARCHAR AS period,
        MAX(Curr) AS suit
      FROM base
      GROUP BY SiteRef, Spp
    ),
    
    -- future rows 
    future_rows AS (
      SELECT
        SiteRef,
        Spp,
        FuturePeriod AS period,
        Newsuit      AS suit
      FROM base
    ),
    
    long AS (
      SELECT * FROM curr_rows
      UNION ALL
      SELECT * FROM future_rows
    ),

    scored AS (
      SELECT
        l.SiteRef,
        l.Spp,
        SUM(l.suit * w.fut_wt)  AS FutSuit,
        SUM(l.suit * w.curr_wt) AS CurrSuit
      FROM long l
      JOIN period_weights w
        ON w.period = l.period
      GROUP BY l.SiteRef, l.Spp
    ),

    flags AS (
      SELECT
        s.*,
        (s.CurrSuit <= 3 AND s.FutSuit > 3.5) AS Loss,
        (s.CurrSuit > 3.5 AND s.FutSuit <= 3) AS Gain,
        (s.CurrSuit <= 3 AND s.FutSuit <= 3) AS Stable
      FROM scored s
    )

    SELECT
      SiteRef,
      SUM(CASE WHEN Loss THEN 1 ELSE 0 END) AS LossSpp,
      SUM(CASE WHEN Gain THEN 1 ELSE 0 END) AS GainSpp,
      SUM(CASE WHEN Stable THEN 1 ELSE 0 END) AS StableSpp,
      SUM(CASE WHEN CurrSuit <= 3 THEN 1 ELSE 0 END) as TotalSpp,
      SUM(CASE WHEN FutSuit <= 3 THEN 1 ELSE 0 END) as FutureSpp
    FROM flags
    GROUP BY SiteRef;
  ", .con = con)
  
  res <- dbGetQuery(con, sql)
  data.table::as.data.table(res)
}


spp_risk <- function(
    con,
    edatope = "C4",
    fut_wt  = c(0,0,0,0.5,0.5,0),
    curr_wt = c(0.5,0.5,0,0,0,0),
    periods = c("Curr","2001_2020","2021_2040","2041_2060","2061_2080","2081_2100")
) {
  stopifnot(length(fut_wt) == length(periods), length(curr_wt) == length(periods))
  
  # create a TEMP weights table in DuckDB
  wts <- data.frame(
    period  = periods,
    fut_wt  = as.numeric(fut_wt),
    curr_wt = as.numeric(curr_wt)
  )
  dbExecute(con, "DROP TABLE IF EXISTS period_weights;")
  dbWriteTable(con, "period_weights", wts, temporary = TRUE)
  
  sql <- glue_sql("
    WITH base AS (
      SELECT
        SiteRef,
        Edatope,
        FuturePeriod,
        Spp,
        Curr,
        Newsuit
      FROM cciss_res
      WHERE Edatope IN ({edatope*})
    ),

    curr_rows AS (
      SELECT
        SiteRef,
        Edatope,
        Spp,
        'Curr'::VARCHAR AS period,
        MAX(Curr) AS suit
      FROM base
      GROUP BY SiteRef, Edatope, Spp
    ),
    
    -- future rows 
    future_rows AS (
      SELECT
        SiteRef,
        Edatope,
        Spp,
        FuturePeriod AS period,
        Newsuit      AS suit
      FROM base
    ),
    
    long AS (
      SELECT * FROM curr_rows
      UNION ALL
      SELECT * FROM future_rows
    ),

    scored AS (
      SELECT
        l.SiteRef,
        l.Spp,
        Edatope,
        SUM(l.suit * w.fut_wt)  AS FutSuit,
        SUM(l.suit * w.curr_wt) AS CurrSuit
      FROM long l
      JOIN period_weights w
        ON w.period = l.period
      GROUP BY l.SiteRef, Edatope, l.Spp
    ),

    flags AS (
      SELECT
        s.*,
        (s.CurrSuit <= 2 AND s.FutSuit > 3.5) AS High,
        (s.CurrSuit <= 2 AND s.FutSuit >= 2.5 AND s.FutSuit < 3.5) AS Mod,
        (s.CurrSuit <= 2 AND s.FutSuit <= 3) AS Low
      FROM scored s
    )

    SELECT
      SiteRef,
      SUM(CASE WHEN High THEN 1 ELSE 0 END) AS HighRisk,
      SUM(CASE WHEN Mod THEN 1 ELSE 0 END) AS ModRisk,
      SUM(CASE WHEN Low THEN 1 ELSE 0 END) AS LowRisk,
      SUM(CASE WHEN CurrSuit <= 2 THEN 1 ELSE 0 END) as TotalP,
      bgc,
      region
    FROM flags
    JOIN bgc_points ON flags.SiteRef = bgc_points.cellnum
    JOIN dist_points ON flags.SiteRef = dist_points.cellnum
    GROUP BY bgc, region, SiteRef;
  ", .con = con)
  
  res <- dbGetQuery(con, sql)
  data.table::as.data.table(res)
}

spp_risk_v2 <- function(
    con,
    edatope = "C4",
    fut_wt  = c(0,0,0,0.5,0.5,0),
    curr_wt = c(0.5,0.5,0,0,0,0),
    periods = c("Curr","2001_2020","2021_2040","2041_2060","2061_2080","2081_2100")
) {
  stopifnot(length(fut_wt) == length(periods), length(curr_wt) == length(periods))
  
  # create a TEMP weights table in DuckDB
  wts <- data.frame(
    period  = periods,
    fut_wt  = as.numeric(fut_wt),
    curr_wt = as.numeric(curr_wt)
  )
  dbExecute(con, "DROP TABLE IF EXISTS period_weights;")
  dbWriteTable(con, "period_weights", wts, temporary = TRUE)
  
  sql <- glue_sql("
    WITH base AS (
      SELECT
        SiteRef,
        Edatope,
        FuturePeriod,
        Spp,
        Curr,
        Newsuit
      FROM cciss_res
      WHERE Edatope IN ({edatope*})
    ),

    curr_rows AS (
      SELECT
        SiteRef,
        Edatope,
        Spp,
        'Curr'::VARCHAR AS period,
        MAX(Curr) AS suit
      FROM base
      GROUP BY SiteRef, Edatope, Spp
    ),
    
    future_rows AS (
      SELECT
        SiteRef,
        Edatope,
        Spp,
        FuturePeriod AS period,
        Newsuit      AS suit
      FROM base
    ),
    
    long AS (
      SELECT * FROM curr_rows
      UNION ALL
      SELECT * FROM future_rows
    ),

    scored AS (
      SELECT
        l.SiteRef,
        l.Spp,
        Edatope,
        SUM(l.suit * w.fut_wt)  AS FutSuit,
        SUM(l.suit * w.curr_wt) AS CurrSuit
      FROM long l
      JOIN period_weights w
        ON w.period = l.period
      GROUP BY l.SiteRef, Edatope, l.Spp
    ),

    flags AS (
      SELECT
        s.*,
        (s.CurrSuit <= 2.5 AND s.FutSuit > 3.5) AS Loss
      FROM scored s
    )

    SELECT
      SiteRef,
      Spp,
      SUM(CASE WHEN Loss THEN 1 ELSE 0 END) AS Loss,
      bgc,
      region
    FROM flags
    JOIN bgc_points ON flags.SiteRef = bgc_points.cellnum
    JOIN dist_points ON flags.SiteRef = dist_points.cellnum
    GROUP BY bgc, region, SiteRef, Spp;
  ", .con = con)
  
  res <- dbGetQuery(con, sql)
  data.table::as.data.table(res)
}



no_suit <- function(
    con,
    edatope = "C4",
    fut_wt  = c(0,0,0,0.5,0.5,0),
    curr_wt = c(0.5,0.5,0,0,0,0),
    periods = c("Curr","2001_2020","2021_2040","2041_2060","2061_2080","2081_2100")
) {
  stopifnot(length(fut_wt) == length(periods), length(curr_wt) == length(periods))
  
  # create a TEMP weights table in DuckDB
  wts <- data.frame(
    period  = periods,
    fut_wt  = as.numeric(fut_wt),
    curr_wt = as.numeric(curr_wt)
  )
  dbExecute(con, "DROP TABLE IF EXISTS period_weights;")
  dbWriteTable(con, "period_weights", wts, temporary = TRUE)
  
  sql <- glue_sql("
    WITH base AS (
      SELECT
        SiteRef,
        Edatope,
        FuturePeriod,
        Spp,
        Curr,
        Newsuit
      FROM cciss_res
      WHERE Edatope IN ({edatope*})
    ),

    curr_rows AS (
      SELECT
        SiteRef,
        Edatope,
        Spp,
        'Curr'::VARCHAR AS period,
        MAX(Curr) AS suit
      FROM base
      GROUP BY SiteRef, Edatope, Spp
    ),
    
    -- future rows 
    future_rows AS (
      SELECT
        SiteRef,
        Edatope,
        Spp,
        FuturePeriod AS period,
        Newsuit      AS suit
      FROM base
    ),
    
    long AS (
      SELECT * FROM curr_rows
      UNION ALL
      SELECT * FROM future_rows
    ),

    scored AS (
      SELECT
        l.SiteRef,
        l.Spp,
        Edatope,
        SUM(l.suit * w.fut_wt)  AS FutSuit,
        SUM(l.suit * w.curr_wt) AS CurrSuit
      FROM long l
      JOIN period_weights w
        ON w.period = l.period
      GROUP BY l.SiteRef, Edatope, l.Spp
    ),

    flags AS (
      SELECT
        s.*,
        (s.FutSuit < 3.5) AS Suitable,
      FROM scored s
    )

    SELECT
      SiteRef,
      SUM(CASE WHEN Suitable THEN 1 ELSE 0 END) AS NumSuit,
      SUM(CASE WHEN CurrSuit < 3.5 THEN 1 ELSE 0 END) as CurrSuit
    FROM flags
    GROUP BY SiteRef;
  ", .con = con)
  
  res <- dbGetQuery(con, sql)
  data.table::as.data.table(res)
}


#' Create table of species consistency 
#' @param con duckdb database connection
#' @param edatope Character. Edatopic position to use. Default `C4` (zonal)
#' @param fut_wt Numeric vector of length 5, corresponding to weight of each period for future period calculation. Default `c(0,0,0,0.5,0.5,0)`, which is average of 2041-2060 and 2061-2080 periods.
#' @param curr_wt Numeric vector of length 5, corresponding to weight of each period for current period calculation. Default `c(0.5,0.5,0,0,0,0)`, which is average of refperiod and 2001-2020 periods.
#' @param periods Character vector of length 5 with names of periods. Must match names in table `cciss_res`. Default `c("Curr","2001_2020","2021_2040","2041_2060","2061_2080","2081_2100")`
#' @return data.table containing loss and gain proportion by Spp and BGC
#' @import data.table duckdb
#' @importFrom glue glue_sql
#' @export
spp_loss_gain_temporal <- function(
    con,
    edatope = "C4"
) {
  
  sql <- glue_sql("
    WITH base AS (
      SELECT
        SiteRef,
        FuturePeriod,
        Spp,
        Curr,
        Newsuit
      FROM cciss_res
      WHERE Edatope = {edatope}
    ),
    
    curr_rows AS (
      SELECT
        SiteRef,
        '1961_1990'::VARCHAR as FuturePeriod,
        Spp,
        Curr,
        Curr AS Newsuit
      FROM base
      WHERE FuturePeriod = '2001_2020'
    ),
    
    long AS (
      SELECT * FROM curr_rows
      UNION ALL
      SELECT * FROM base
    ),

    flags AS (
      SELECT
        s.*,
        (CASE WHEN s.Curr <= 3 AND s.Newsuit > 3.5 THEN 1 ELSE 0 END) AS Loss,
        (CASE WHEN s.Curr > 3.5 AND s.Newsuit <= 3 THEN 1 ELSE 0 END) AS Gain,
        (CASE WHEN s.Curr <= 3 AND s.Newsuit <= 3 THEN 1 ELSE 0 END) AS Stable
      FROM long s
    )

    SELECT
      Spp, FuturePeriod, SUM(Loss) AS Loss, SUM(Gain) AS Gain, SUM(Stable) AS Stable
    FROM flags
    GROUP BY Spp, FuturePeriod;
  ", .con = con)
  
  res <- dbGetQuery(con, sql)
  data.table::as.data.table(res)
}

