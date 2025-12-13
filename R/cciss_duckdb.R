### CCISS duckDb functions

#' Connect to (and possibly create) DuckDB for storing results
#' @param db_file Character. Name of database file, or path. Must end in .duckdb
#' @param read_only Logical. Default FALSE
#' @param threads Numeric. Number of threads to use in database.
#' @param verbose Logical. Be chatty?
#' @returns DuckDb database connection
#' @export
dbCon_cciss <- function(
    db_file,
    read_only   = FALSE,
    threads     = NULL,
    verbose     = TRUE
) {
  
  # --- sanity checks ----------------------------------------------------------
  db_file <- normalizePath(db_file, mustWork = FALSE)
  db_dir  <- dirname(db_file)
  
  if (!dir.exists(db_dir)) {
    dir.create(db_dir, recursive = TRUE)
  }
  
  db_exists <- file.exists(db_file)
  
  if (verbose) {
    if (db_exists) {
      message("✓ Opening DuckDB database: ", db_file)
    } else {
      message("▶ Creating DuckDB database: ", db_file, "; You will need to run `dbPopulate`.")
    }
  }
  
  # --- connect ---------------------------------------------------------------
  con <- DBI::dbConnect(
    duckdb::duckdb(),
    dbdir    = db_file,
    read_only = read_only
  )
  
  # --- pragmas ---------------------------------------------------------------
  if (!is.null(threads)) {
    DBI::dbExecute(con, sprintf("PRAGMA threads=%d;", threads))
  }
  
  # Optional but useful defaults
  DBI::dbExecute(con, "PRAGMA enable_progress_bar=true;")
  return(con)
}

#' Populate database with bgc point, edatopic and suitability tables
#' @param dbCon Database connection
#' @param bgc_template
#' @param edatopes Character vector of edatopes
#' @export
dbPopulate <- function(dbCon, bgc_template, edatopes = c("B2","C4","D6")) {
  
  bgc_rast <- bgc_template$bgc_rast ##bgc_points
  rast_ids <- bgc_template$ids
  bgc_points <- as.data.frame(bgc_rast, cells=T) |> as.data.table()
  bgc_points[rast_ids, bgc := i.bgc, on = "bgc_id"]
  bgc_points[,bgc_id := NULL]
  setnames(bgc_points, old = "cell", new = "cellnum")
  dbWriteTable(dbCon, "bgc_points", bgc_points, row.names = F, overwrite = TRUE)
  
  eda_table <- copy(E1) ##Edatopic table
  eda_table <- eda_table[is.na(SpecialCode),]
  eda_table <- eda_table[Edatopic %in% edatopes,]
  eda_table <- unique(eda_table[,.(BGC,SS_NoSpace,Edatopic)])
  dbWriteTable(dbCon, "edatopic", eda_table, row.names = FALSE, overwrite = TRUE)
  
  suit <- copy(S1) ##Suitability table
  suit <- na.omit(suit, cols = "spp")
  suit <- suit[,.(spp,ss_nospace,newfeas)]
  dbWriteTable(dbCon, "suitability", suit, row.names = FALSE, overwrite = TRUE)
  
  message("✓ Written bgc_points, edatopic, and suitability tables to duckdb!")
  return(invisible(TRUE))
}

#' Summarise raw BGC predictions
#' @param ssp_use Character. List of ssps to use. Default `c("ssp126", "ssp245", "ssp370")`
#' @param ssp_w Numeric vector. Weights for each ssp in `ssp_use`
#' @param base_folder Character. Name of base folder to write results to.
#' @return NULL. Results are written to csv files in base_folder/bgc_data
#' @import data.table
#' @export
summarise_preds <- function(dbCon,
                            ssp_use = c("ssp126", "ssp245", "ssp370"),
                            ssp_w = c(0.8,1,0.8)) {
  ssp_weights <- data.table(ssp = ssp_use, weight = ssp_w)
  dbWriteTable(dbCon, "ssp_weights", ssp_weights, temporary = TRUE)
  
  summarised_qry <- "CREATE TABLE bgc_summary AS
                      WITH wt AS (
                        SELECT r.*, w.weight
                        FROM bgc_raw r
                        LEFT JOIN ssp_weights w USING (ssp)
                      ),
                      totals AS (
                        SELECT cellnum, period, SUM(weight) AS tot_wt
                        FROM wt
                        GROUP BY cellnum, period
                      )
                      SELECT
                        w.cellnum,
                        w.period,
                        w.bgc_pred,
                        SUM(w.weight) / t.tot_wt AS bgc_prop
                      FROM wt w
                      JOIN totals t
                        ON w.cellnum = t.cellnum
                       AND w.period = t.period
                      GROUP BY w.cellnum, w.period, w.bgc_pred, t.tot_wt;
                      "
  dbExecute(con, summarised_qry)
  message("✓ Created table bgc_summary.")
  return(invisible(TRUE))
}


#' Create table of BGC subzone/variant persistance/expansion from raw BGC projections
#' @param bgc_mapped List containing SpatRaster of BGCs and id table, or data.table (must have columns `cell`,`BGC`). Usually created using `make_bgc_template`
#' @param periods Character. Either "auto" or vector of period names. If "auto", all periods present in bgc projections will be used.
#' @param base_folder Base folder to read results from.
#' @return data.table containing persistance and expansion for each model/period/scenario/bgc
#' @import data.table
#' @export
bgc_persist_expand <- function(dbCon, by_zone = TRUE){
  if(by_zone){
    tbl_nm <- "bgc_per_exp_zone"
  } else {
    tbl_nm <- "bgc_per_exp_sz"
  }
  
  if(!duckdb_table_exists(dbCon, tbl_nm)) {
    message("▶ Computing bgc persistence !")
    
    new_res <- calc_bgc_persist_expand(
      dbCon,
      by_zone = by_zone
    )
    
    # append to cache
    dbWriteTable(
      dbCon,
      tbl_nm,
      new_res
    )
  } 
  
  res <- dbGetQuery(dbCon, sprintf("
    SELECT *
    FROM %s
  ", tbl_nm)) |> as.data.table()
  return(res)
}



calc_bgc_persist_expand <- function(dbCon, by_zone = TRUE){
  stopifnot(DBI::dbIsValid(dbCon))
  
  ## --- Choose grouping variables depending on level ----
  
  if (by_zone) {
    pred_expr <- "regexp_extract(a.bgc_pred, '^[A-Z]+')"
    true_expr <- "regexp_extract(p.bgc,      '^[A-Z]+')"
    tot_group <- "regexp_extract(bgc, '^[A-Z]+')"
  } else {
    pred_expr <- "a.bgc_pred"
    true_expr <- "p.bgc"
    tot_group <- "bgc"
  }
  
  ## Construct SQL dynamically
  sql <- sprintf("
    WITH
    joined AS (
      SELECT
        a.cellnum,
        a.ssp,
        a.gcm,
        a.run,
        a.period,
        %s AS bgc_pred,
        %s AS bgc
      FROM bgc_raw a
      LEFT JOIN bgc_points p USING (cellnum)
      WHERE p.bgc IS NOT NULL          -- removes NA join
    ),
    flags AS (
      SELECT
        *,
        CASE WHEN bgc_pred = bgc THEN 1 ELSE 0 END AS Persist,
        CASE WHEN bgc_pred <> bgc THEN 1 ELSE 0 END AS Expand
      FROM joined
    ),
    agg AS (
      SELECT
        ssp, gcm, run, period, bgc_pred,
        SUM(Persist) AS Persist_Tot,
        SUM(Expand)  AS Expand_Tot
      FROM flags
      GROUP BY ssp, gcm, run, period, bgc_pred
    ),
    bgc_tot AS (
      SELECT
        %s AS bgc_true,
        COUNT(*) AS BGC_Tot
      FROM bgc_points
      GROUP BY bgc_true
    )
    SELECT
      a.*,
      t.BGC_Tot,
      (Persist_Tot * 1.0) / t.BGC_Tot AS Persistance,
      (Expand_Tot  * 1.0) / t.BGC_Tot AS Expansion
    FROM agg a
    LEFT JOIN bgc_tot t
      ON a.bgc_pred = t.bgc_true;
      ",
      # substitutions:
      pred_expr, true_expr,
      tot_group
      )
    
    # Run in DuckDB
  res <- dbGetQuery(dbCon, sql) |> as.data.table()
  return(res)
}


duckdb_table_exists <- function(dbCon, table_name) {
  sql <- sprintf("
    SELECT COUNT(*) > 0 AS exists
    FROM information_schema.tables
    WHERE table_name = '%s'
  ", table_name)
  
  DBI::dbGetQuery(dbCon, sql)$exists[1]
}

materialise_bgc_eda <- function(
    dbCon,
    rebuild          = FALSE
) {
  
  bgc_all_table    = "bgc_raw"
  bgc_points_table = "bgc_points"
  eda_table        = "edatopic"
  out_table        = "bgc_eda_mat"
  
  if (!rebuild && duckdb_table_exists(dbCon, out_table)) {
    message(sprintf("✓ Using existing table '%s'", out_table))
    return(invisible(FALSE))
  }
  
  message(sprintf("▶ Creating materialised table '%s' ...", out_table))
  
  dbExecute(dbCon, sprintf("DROP TABLE IF EXISTS %s;", out_table))
  
  sql <- sprintf("
    CREATE TABLE %s AS
    SELECT
      a.cellnum,
      a.ssp,
      a.gcm,
      a.run,
      a.period,
      a.bgc_pred,
      p.bgc AS bgc_true,
      ep.Edatopic,
      ep.SS_NoSpace AS SS_Pred,
      eh.SS_NoSpace AS SS_NoSpace
    FROM %s a
    JOIN %s p USING (cellnum)
    LEFT JOIN %s ep
      ON a.bgc_pred = ep.BGC
    LEFT JOIN %s eh
      ON p.bgc = eh.BGC
     AND ep.Edatopic = eh.Edatopic
    WHERE p.bgc IS NOT NULL;
  ",
                 out_table,
                 bgc_all_table,
                 bgc_points_table,
                 eda_table,
                 eda_table
  )
  
  DBI::dbExecute(dbCon, sql)
  
  message("✓ bgc_eda materialised")
  invisible(TRUE)
}


calc_spp_persist_expand <- function(
    con,
    spp_list,
    fractional = TRUE
) {
  
  stopifnot(DBI::dbIsValid(con))
  
  spp_sql <- paste(sprintf("'%s'", spp_list), collapse = ",")
  
  new_suit_expr <- if (fractional) {
    "1.0 - (NewSuit_code_clean - 1) / 4.0"
  } else {
    "CASE WHEN NewSuit_code_clean <> 5 THEN 1 ELSE 0 END"
  }
  
  hist_suit_expr <- if (fractional) {
    "1.0 - (HistSuit_code_clean - 1) / 4.0"
  } else {
    "CASE WHEN HistSuit_code_clean <> 5 THEN 1 ELSE 0 END"
  }
  
  mapped_expr <- if (fractional) {
    "SUM( (1.0 - (Suit_code - 1) / 4.0) * BGC_Tot )"
  } else {
    "SUM( (CASE WHEN Suit_code = 5 THEN 0 ELSE 1 END) * BGC_Tot )"
  }
  materialise_bgc_eda(con)
  
  sql <- sprintf("
    ---------------------------------------------------------------------
    -- 3. Compute mapped historic suitability
    ---------------------------------------------------------------------
    WITH bgc_sum AS (
      SELECT bgc, COUNT(*) AS BGC_Tot
      FROM bgc_points
      GROUP BY bgc
    ),

    mapped_raw AS (
      SELECT
        bs.bgc,
        bs.BGC_Tot,
        e.Edatopic,
        e.SS_NoSpace,
        s.spp,
        s.newfeas
      FROM bgc_sum bs
      LEFT JOIN edatopic e ON bs.bgc = e.BGC
      LEFT JOIN suitability s ON e.SS_NoSpace = s.ss_nospace
      WHERE s.spp IN (%s)
    ),

    mapped_clean AS (
      SELECT
        spp,
        bgc,
        BGC_Tot,
        Edatopic,
        CASE WHEN newfeas IS NULL OR newfeas = 4 THEN 5 ELSE newfeas END AS Suit_code
      FROM mapped_raw
    ),

    mapped_min AS (
      SELECT
        spp, bgc, BGC_Tot, Edatopic,
        MIN(Suit_code) AS Suit_code
      FROM mapped_clean
      GROUP BY spp, bgc, BGC_Tot, Edatopic
    ),

    mapped_suit AS (
      SELECT
        spp, Edatopic,
        %s AS MappedSuit --fractional or binary
      FROM mapped_min
      GROUP BY spp, Edatopic
    ),

    ---------------------------------------------------------------------
    -- 4. Assign BOTH historic and projected suitability
    ---------------------------------------------------------------------
    bgc_spp AS (
      SELECT
        st.spp,
        e.cellnum,
        e.ssp,
        e.gcm,
        e.run,
        e.period,
        e.Edatopic,
        --
        -- NewSuit raw
        --
        CASE
          WHEN st.newfeas IS NULL OR st.newfeas = 4 THEN 5
          ELSE st.newfeas
        END AS NewSuit_code_clean,
        --
        -- Historic Suit (based on bgc_true)
        --
        CASE
          WHEN sth.newfeas IS NULL OR sth.newfeas = 4 THEN 5
          ELSE sth.newfeas
        END AS HistSuit_code_clean
      FROM bgc_eda_mat e
      JOIN suitability st
        ON e.SS_Pred = st.ss_nospace
        AND st.spp IN (%s)
        
      LEFT JOIN suitability sth
        ON e.SS_NoSpace = sth.ss_nospace
        AND sth.spp = st.spp
    ),

    bgc_val AS (
      SELECT
        spp, cellnum, ssp, gcm, run, period, Edatopic,
        %s AS NewSuit_val, --fractional or binary
        %s AS HistSuit_val
      FROM bgc_spp
    ),

    ---------------------------------------------------------------------
    -- 5. Max suitability per location
    ---------------------------------------------------------------------
    bgc_best AS (
      SELECT
        spp, cellnum, ssp, gcm, run, period, Edatopic,
        MAX(NewSuit_val)  AS NewSuit,
        MAX(HistSuit_val) AS HistSuit
      FROM bgc_val
      GROUP BY spp, cellnum, ssp, gcm, run, period, Edatopic
    ),

    ---------------------------------------------------------------------
    -- 6. Persistence + Expansion
    ---------------------------------------------------------------------
    perexp_raw AS (
      SELECT
        spp, Edatopic, ssp, gcm, run, period,
        CASE WHEN HistSuit > 0 THEN NewSuit ELSE 0 END AS Persist,
        CASE WHEN HistSuit = 0 THEN NewSuit ELSE 0 END AS Expand
      FROM bgc_best
    ),

    perexp_tot AS (
      SELECT
        spp, Edatopic, ssp, gcm, run, period,
        SUM(Persist) AS Persist_Tot,
        SUM(Expand)  AS Expand_Tot
      FROM perexp_raw
      GROUP BY spp, Edatopic, ssp, gcm, run, period
    )

    ---------------------------------------------------------------------
    -- 7. Join with MappedSuit and normalise
    ---------------------------------------------------------------------
    SELECT
      p.*,
      m.MappedSuit,
      p.Persist_Tot / m.MappedSuit AS Persistance,
      p.Expand_Tot  / m.MappedSuit AS Expansion
    FROM perexp_tot p
    JOIN mapped_suit m
      ON p.spp = m.spp AND p.Edatopic = m.Edatopic;
    ",
                 # INSERTS:
                 spp_sql,
                 mapped_expr,
                 spp_sql,
                 new_suit_expr,
                 hist_suit_expr
  )
  
  res <- dbGetQuery(con, sql)
  return(as.data.table(res))
}

spp_persist_expand <- function(dbCon, spp_list, fractional = TRUE) {
  if(fractional){
    tbl_nm <- "spp_per_exp_frac"
  } else {
    tbl_nm <- "spp_per_exp_bin"
  }
  
  spp_sql <- paste(sprintf("'%s'", spp_list), collapse = ",")
  
  if(duckdb_table_exists(dbCon, tbl_nm)) {
    cached_spp <- dbGetQuery(dbCon, sprintf("
                           SELECT DISTINCT spp
                           FROM %s
                           WHERE spp IN (%s)", tbl_nm, spp_sql))
    missing_spp <- setdiff(spp_list, cached_spp$spp)
  } else {
    missing_spp <- spp_list
  }
  
  if (length(missing_spp) > 0) {
    message("▶ Computing persistence for: ", paste(missing_spp, collapse = ", "))
    
    # run your existing DuckDB function on missing spp only
    new_res <- calc_spp_persist_expand(
      dbCon,
      spp_list   = missing_spp,
      fractional = fractional
    )
    
    # append to cache
    dbWriteTable(
      dbCon,
      tbl_nm,
      new_res,
      append = TRUE
    )
  } else {
    message("✓ All requested species already cached")
  }
  
  res <- dbGetQuery(dbCon, sprintf("
    SELECT *
    FROM %s
    WHERE spp IN (%s)
  ", tbl_nm, spp_sql)) |> as.data.table()
  return(res)
}


calc_suit_area <- function(
    con,
    spp_list,
    fractional = TRUE
) {
  
  stopifnot(DBI::dbIsValid(con))
  
  # Quote species list for SQL IN()
  spp_sql <- paste(sprintf("'%s'", spp_list), collapse = ",")
  
  # Suit transformation for fractional vs binary
  new_suit_expr <- if (fractional) {
    "1.0 - (NewSuit_code - 1) / 4.0"
  } else {
    "CASE WHEN NewSuit_code <> 5 THEN 1 ELSE 0 END"
  }
  
  mapped_expr <- if (fractional) {
    "SUM( (1.0 - (Suit_code - 1) / 4.0) * BGC_Tot )"
  } else {
    "SUM( (CASE WHEN Suit_code <> 5 THEN 1 ELSE 0 END) * BGC_Tot )"
  }
  
  sql <- sprintf("
    -----------------------------------------------------------------------------
    -- 1. Add BGC_TRUE FROM points
    -----------------------------------------------------------------------------
    WITH bgc_joined AS (
      SELECT a.*, p.bgc AS bgc_true
      FROM bgc_raw a
      LEFT JOIN bgc_points p USING (cellnum)
      WHERE p.bgc IS NOT NULL
    ),

    -----------------------------------------------------------------------------
    -- 2. Merge edatopic on bgc_pred (SS_Pred)
    -----------------------------------------------------------------------------
    bgc_eda AS (
      SELECT
        b.*,
        e1.Edatopic,
        e1.SS_NoSpace AS SS_Pred
      FROM bgc_joined b
      LEFT JOIN edatopic e1
        ON b.bgc_pred = e1.BGC
    ),

    -----------------------------------------------------------------------------
    -- 3. Compute mapped historic suitability per species × edatopic
    -----------------------------------------------------------------------------
    bgc_sum AS (
      SELECT bgc, COUNT(*) AS BGC_Tot
      FROM bgc_points
      GROUP BY bgc
    ),

    mapped_raw AS (
      SELECT
        bs.bgc,
        bs.BGC_Tot,
        e.Edatopic,
        e.SS_NoSpace,
        s.spp,
        s.newfeas
      FROM bgc_sum bs
      FULL JOIN edatopic e ON bs.bgc = e.bgc
      FULL JOIN suitability s ON e.ss_nospace = s.ss_nospace
      WHERE s.spp IN (%s)
    ),

    mapped_clean AS (
      SELECT
        spp,
        bgc,
        BGC_Tot,
        Edatopic,
        CASE WHEN newfeas IS NULL OR newfeas = 4 THEN 5 ELSE newfeas END AS Suit_code
      FROM mapped_raw
    ),

    mapped_min AS (
      SELECT
        spp, bgc, BGC_Tot, Edatopic,
        MIN(Suit_code) AS Suit_code
      FROM mapped_clean
      GROUP BY spp, bgc, BGC_Tot, Edatopic
    ),

    mapped_suit AS (
      SELECT
        spp,
        Edatopic,
        %s AS MappedSuit --using fractional or binary
      FROM mapped_min
      GROUP BY spp, Edatopic
    ),

    -----------------------------------------------------------------------------
    -- 4. Compute projected suitability for each species
    -----------------------------------------------------------------------------
    bgc_spp_suit AS (
      SELECT
        st.spp,
        e.cellnum,
        e.ssp,
        e.gcm,
        e.run,
        e.period,
        e.Edatopic,
        --
        -- Assign suitability from suit table
        --
        CASE
          WHEN st.newfeas IS NULL OR st.newfeas = 4 THEN 5
          ELSE st.newfeas
        END AS NewSuit_code
      FROM bgc_eda e
      LEFT JOIN suitability st
        ON e.SS_Pred = st.ss_nospace
      WHERE st.spp IN (%s)
    ),

    bgc_spp_val AS (
      SELECT
        spp,
        cellnum,
        ssp,
        gcm,
        run,
        period,
        Edatopic,
        %s AS NewSuit_val --fractional or binary
      FROM bgc_spp_suit
    ),

    -----------------------------------------------------------------------------
    -- 5. Best suitability per location
    -----------------------------------------------------------------------------
    bgc_best AS (
      SELECT
        spp, cellnum, ssp, gcm, run, period, Edatopic,
        MAX(NewSuit_val) AS NewSuit
      FROM bgc_spp_val
      GROUP BY spp, cellnum, ssp, gcm, run, period, Edatopic
    ),

    -----------------------------------------------------------------------------
    -- 6. Area per species × edatopic × scenario
    -----------------------------------------------------------------------------
    suit_area_raw AS (
      SELECT
        spp, Edatopic, ssp, gcm, run, period,
        SUM(NewSuit) AS Proj_Area
      FROM bgc_best
      GROUP BY spp, Edatopic, ssp, gcm, run, period
    )

    -----------------------------------------------------------------------------
    -- 7. Join with mapped suit & compute Suit_Prop
    -----------------------------------------------------------------------------
    SELECT
      s.spp,
      s.Edatopic,
      s.ssp,
      s.gcm,
      s.run,
      s.period,
      s.Proj_Area,
      m.MappedSuit,
      s.Proj_Area / m.MappedSuit AS Suit_Prop
    FROM suit_area_raw s
    JOIN mapped_suit m
      ON s.spp = m.spp AND s.Edatopic = m.Edatopic;
  ",
                 # inserts:
                 spp_sql,
                 mapped_expr,
                 spp_sql,
                 new_suit_expr
  )
  
  res <- dbGetQuery(con, sql)
  return(as.data.table(res))
}

spp_suit_area <- function(dbCon, spp_list, fractional = TRUE) {
  if(fractional){
    tbl_nm <- "spp_suit_area_frac"
  } else {
    tbl_nm <- "spp_suit_area_bin"
  }
  
  spp_sql <- paste(sprintf("'%s'", spp_list), collapse = ",")
  
  if(duckdb_table_exists(dbCon, tbl_nm)) {
    cached_spp <- dbGetQuery(dbCon, sprintf("
                           SELECT DISTINCT spp
                           FROM %s
                           WHERE spp IN (%s)", tbl_nm, spp_sql))
    missing_spp <- setdiff(spp_list, cached_spp$spp)
  } else {
    missing_spp <- spp_list
  }
  
  if (length(missing_spp) > 0) {
    message("▶ Computing suitable area for: ", paste(missing_spp, collapse = ", "))
    
    # run your existing DuckDB function on missing spp only
    new_res <- calc_suit_area(
      dbCon,
      spp_list   = missing_spp,
      fractional = fractional
    )
    
    # append to cache
    dbWriteTable(
      dbCon,
      tbl_nm,
      new_res,
      append = TRUE
    )
  } else {
    message("✓ All requested species already cached")
  }
  
  res <- dbGetQuery(dbCon, sprintf("
    SELECT *
    FROM %s
    WHERE spp IN (%s)
  ", tbl_nm, spp_sql)) |> as.data.table()
  return(res)
}

