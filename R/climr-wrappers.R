#' Prepares coordinates and obtains climate normals
#'  using `climr_downscale`
#'
#' @param coords a `data.table`, or spatial points (`SpatVector` or `sf`) 
#'   with point coordinates ("lon" = longitude, "lat" = latitude), elevation ("elev") and point IDs ("id").
#' @param ... further arguments passed to [climr_downscale()].
#' 
#' @details
#' If `bgcs` is provided, the BGC field will be appended to
#' the output `data.table`.
#' 
#' @seealso [climr_downscale()]
#'
#' @return climate normals as a `data.table`
#'  
#' @importFrom terra intersect
#' @importFrom data.table data.table as.data.table setnames
#' @importFrom methods is
#' @export
getClimate <- function(coords, ...) {
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
  
  args <- append(list(coords = coords), dots)
  out <- do.call(.getClimVars, args)
  
  return(out)
}


#' Wrapper for `climr_downscale` that keeps extra columns in
#' table of point coordinates.
#'
#' @inheritParams getClimate 
#' @param byCombo logical. If TRUE, `climr_downscale` is iterated by 
#'   combinations of GCM models, periods and scenarios.
#' @param outFormat character. Should outputs be in the form of a 
#'  `data.table` ("data.table"), list of `data.tables` ("list") or 
#'  written directly to disk ("disk")?
#' @param filename character. Passed to `write.csv(..., file)` if
#'   `outFormat == "disk"`. Defaults to `tempfile(fileext = ".csv")`.
#'   
#' @details
#'   If `outFormat == "disk"` and `byCombo == TRUE`, `climr_downscale` 
#'   is iterated for combinations of GCMs, runs, periods and scenarios,
#'   and output `data.tables` are saved to a .csv file with
#'   `write.csv(..., file = filename, append = TRUE)`.
#' 
#' @return climate normals returned as a `data.table`, list of `data.tables`.
#'   If `outFormat == "disk"`, `filename` is returned instead.
#' 
#' @importFrom climr climr_downscale
#' @importFrom data.table setDT
.getClimVars <- function(coords, coords_bgc, byCombo = FALSE, outFormat = "data.table",
                         filename = tempfile(fileext = ".csv"), ...) {
  browser()
  dots <- list(...)
  if (byCombo) {
    dots <- as.data.table(expand.grid(dots, stringsAsFactors = FALSE))
    
    for (i in 1:nrow(dots)) {
      x <- dots[i]
      out <- do.call(climr_downscale, append(list(xyz = coords), x))
    }
  } else {
    
  }
  
  if (outFormat == "data.table") {
    
  }
  
  ##test: 
  # clim_vars <- climr_downscale(coords[1003050:1003085,], ...)
  clim_vars <- clim_vars[!is.nan(PPT05),] ##lots of points in the ocean
  clim_vars[, PERIOD := NULL]
  
  setDT(clim_vars)
  
  ## join all columns back
  clim_vars <- coords[clim_vars, on = "id"]
  
  if (saveToDisk) {
    
  }
  
  return(clim_vars)
}