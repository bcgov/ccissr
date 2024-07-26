#' Prepares coordinates and obtains climate normals
#'  using `climr::downscale`
#'
#' @inheritParams .getClimVars 
#' 
#' @details
#' If `bgcs` is provided, the BGC field will be appended to
#' the output `data.table`.
#' 
#' @seealso [climr::downscale()]
#'
#' @return climate normals as a `data.table`
#'  
#' @importFrom terra intersect
#' @importFrom data.table data.table as.data.table setnames
#' @importFrom methods is
#' @export
getClimate <- function(coords, byCombo = FALSE, outFormat = "data.table",
                       filename = NULL, ...) {
  dots <- list(...)
  
  if (!is(coords, "data.table")) {
    ## keep all columns, but rename the ones created in the conversion
    coords <- if (is(coords, "SpatVector")) {
      if (!compareGeom(coords, vect(crs = crs("EPSG:4326", proj = TRUE)), 
                       stopOnError = FALSE)) {
        stop("coords must be in EPSG:4326 projection.")
      }
      
      try(as.data.table(coords, geom = "XY")) |>
        setnames(old = c("x", "y"), c("lon", "lat"))
    } else if (inherits(coords, "sf")) {
      if (requireNamespace("sf")) {
        if (sf::st_crs(coords) != sf::st_crs("EPSG:4326")) {
          stop("coords must be in EPSG:4326 projection")
        }
        
        coords <- try({
          cbind(sf::st_drop_geometry(coords), sf::st_coordinates(coords)) |>
            as.data.table() |>
            setnames(old = c("X", "Y"), c("lon", "lat")) 
        })
      } else {
        stop("Please install 'sf' package or provide a data.table or SpatVector to coords.")
      }
    } 
    
    if (is(coords, "simple-error")) {
      stop("coords is not coercible to 'data.table'. Provide coercible object class",
           " or a 'data.table'")
    }
  }
  
  climrCols <- c("lon", "lat", "elev", "id")
  if (any(!climrCols %in% names(coords))) {
    stop("coords must contain columns 'lon', 'lat', 'elev' and 'id'", 
         "\n  If providing a SpatVector/sf object, ensure it has 'elev' and 'id' attributes")
  }
  
  args <- append(list(coords = coords,
                      byCombo = byCombo, 
                      outFormat = outFormat,
                      filename = filename), 
                 dots)
  out <- do.call(.getClimVars, args)

  return(out)
}


#' Wrapper for `climr::downscale` that keeps extra columns in
#' table of point coordinates.
#' 
#' @param coords a `data.table`, or spatial points (`SpatVector` or `sf`) 
#'   with point coordinates ("lon" = longitude, "lat" = latitude), elevation ("elev") and point IDs ("id").
#' @param byCombo logical. If TRUE, `climr::downscale` is iterated by 
#'   combinations of any arguments passed to `downscale` via  `...` 
#'   (e.g., `vars`, `gcms`, `ssps`).
#' @param outFormat character. Should outputs be in the form of a 
#'  `data.table` ("data.table"), list of `data.tables` ("list") or 
#'  written directly to disk ("disk")? If `byCombo == FALSE`, only 
#'  "data.table" and "disk" will work.
#' @param filename character. Passed to `saveRDS(..., file)` if
#'   `outFormat == "disk"`. Defaults to `tempfile(fileext = ".csv")`.
#' @param ... further arguments passed to [climr::downscale()].
#'   
#' @details
#'   If `outFormat == "disk"` and `byCombo == TRUE`, `climr::downscale` 
#'   is iterated for combinations of GCMs, runs, periods and scenarios,
#'   and each output `data.table` is saved to an `.rds` file with
#'   `saveRDS(..., file = filename2, append = TRUE)`, where `basename(filename2)`
#'   is the `basename(filename)` with appended values of the arguments that were iterated
#'   through -- excluding the arguments "cache", "xyz", "nthread", "max_run", 
#'   "return_refperiod", "ppt_lr", "out_spatial", and "plot".
#' 
#' @return climate normals returned as a `data.table` or list of `data.tables`.
#'   If `outFormat == "disk"`, the file name(s) are returned instead.
#' 
#' @importFrom climr downscale
#' @importFrom data.table setDT
.getClimVars <- function(coords, byCombo = FALSE, outFormat = "data.table",
                         filename = tempfile(fileext = ".rds"), ...) {
  ## checks
  if (is.null(filename)) {
    filename <- tempfile(fileext = ".rds")
  }
  
  filename <- normPath(filename)
  
  if (outFormat == "list" & isFALSE(byCombo)) {
    stop("byCombo is FALSE, please set outFormat to 'data.table' or 'disk'")
  }
  
  
  dots <- list(...)
  if (byCombo) {
    dots <- as.data.table(expand.grid(dots, stringsAsFactors = FALSE))
    out <- list()
    for (i in 1:nrow(dots)) {
      x <- dots[i]
      outTemp <- do.call(downscale, append(list(xyz = coords), x))
      
      outTemp <- outTemp[!is.nan(get(x$vars)),] ## remove missing data (e.g. ocean)
      setDT(outTemp)
      
      ## join all columns back
      outTemp <- coords[outTemp, on = "id"]
      
      if (outFormat == "disk") {
        ## exclude certain downscale/downscale_core arguments from the file naming
        argsForFile <- setdiff(colnames(dots),
                               c("cache", "xyz", "nthread", "max_run", "return_refperiod",
                                 "ppt_lr", "out_spatial", "plot"))
        filename2 <- sub("\\.rds", "", basename(filename))
        filename2 <- paste0(filename2, "_", paste(x[, ..argsForFile], collapse = "_"), ".rds")
        filename2 <- file.path(dirname(filename), filename2)
        saveRDS(outTemp, file = filename2)
        
        out[[i]] <- filename2
      } else {
        out[[i]] <- outTemp
      }
    }
    
    if (outFormat == "data.table") {
      sharedCols <- lapply(out, names) |>
        Reduce(intersect, x = _)
      out <- lapply(out, setkeyv, cols = sharedCols)
      
      ## give a common name to climate var col (the one that is not shared across all tables)
      out <- lapply(out, .renameUniqueCol, otherCols = sharedCols)
      
      out <- rbindlist(out, use.names = TRUE)
      
      ## back to wide format
      out <- dcast.data.table(out, ... ~ var, value.var = "value")
    }
  } else {
    out <- do.call(downscale, append(list(xyz = coords), dots))
    
    out <- out[!is.nan(get(dots$vars[1])),] ## remove missing data (e.g. ocean)
    out[, PERIOD := NULL]
    
    setDT(out)
    
    ## join all columns back
    out <- coords[out, on = "id"]
    
    if (outFormat == "disk") {
      saveRDS(out, file = filename)
      out <- filename
    }
  }
  
  return(out)
}


#' Rename a column 
#'
#' @param x a `data.table`.
#' @param otherCols character. Columns not to be renamed
#'
#' @return x with target column renamed as "value" and a 
#'   new column "var" containing the original name of the 
#'   target column
#'
#' @importFrom data.table setnames
.renameUniqueCol <- function(x, otherCols) {
  targetCol <- setdiff(names(x), otherCols)
  
  if (length(targetCol) > 1) stop("There should only be one column to rename.")
  
  setnames(x, targetCol, "value")
  x[, var := targetCol]
  x
}
