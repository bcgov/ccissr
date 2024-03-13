#' Create a set of training points with associated
#'   elevation and BGC values.
#'
#' @param bgc_poly an `sf` object (or one cohercible to an `sf` object) of BGC polygons.
#' @param elev a `SpatRaster` or `RasterLayer` of elevation covering the extent of `bgc_poly`.
#' @param gridSize numeric. Distance in m between points.
#' @param crs passed to [sf::st_as_sf()]
#' 
#' @details Points are sampled regularly from a grid with cell size
#'  defined by `gridSize` that covers `bgc_poly`.
#'
#' @return a `data.table` of point coordinates with associated IDs,
#'   elevation and BGCs.
#'   
#' @importFrom terra extract vect geom rast
#' @importFrom data.table setDT
#' @importFrom methods is
#'
#' @export
makePointCoords <- function(bgc_poly, elev, gridSize = 2000, crs = "EPSG:4326") {
  if (!is(bgc_poly, "SpatVector")) {
    bgc_poly <- try(vect(bgc_poly))
    if (is(bgc_poly, "simple-error")) {
      stop("Can't coherce bgc_poly to SpatVector. Please pass a SpatVector or another cohercible object class.")
    }
  }
  
  if (!is(elev, "SpatRaster")) {
    elev <- try(rast(elev))
    if (is(elev, "simple-error")) {
      stop("Can't coherce elev to a SpatRaster. Please pass SpatRaster or RasterLayer.")
    }
  }
  
  ## faster than using sf::st_make_grid
  # tmp_elev <- 
  bgc_grid <- rast(res = gridSize, extent = ext(bgc_poly), crs = crs(bgc_poly, proj = TRUE))
  bgc_grid[] <- 1L
  bgc_grid <- as.data.table(bgc_grid, xy = TRUE) |>
    vect(geom = c("x", "y"), crs = crs(bgc_grid, proj = TRUE)) ## faster than extracting with DT
  
  coords_elev <- terra::extract(elev, bgc_grid, method = "bilinear", xy = TRUE) |> 
    as.data.table()  ## faster than terra::geom
  setnames(coords_elev, c("id", "elev", "lon", "lat"))
  
  return(coords_elev)
}


#' Subset a group of points by spatial extent
#'
#' @param coords `data.table` of point coordinates with columns "x" (longitude)
#'   and "y" (latitude) (and any additional columns), a `SpatVector` or object 
#'   cohersible to `SpatVector`.
#' @param cropExt `SpatExtent` or `SpatVector` (whose extent will be used)
#'   to subset the data to. Defaults to the `SpatExtent` of an area in Southern BC.
#' @param crs passed to [terra::vect()] to coerce coords to a `SpatVector` 
#'   if it is not one already.
#'
#' @return a subset `coords` object.
#' 
#' @importFrom terra vect crop
#' @importFrom data.table as.data.table
#' @importFrom methods is
#'
#' @export
subsetByExtent <- function(coords, cropExt = ext(c(-123, -118, 49, 52)), crs = "EPSG:4326") {
  
  isSpatial <- FALSE
  
  if (!inherits(cropExt, c("SpatVector", "SpatExtent"))) {
    stop("cropExt must be a SpatVector or a SpatExtent")
  } 
  
  if (is(coords, "data.table")) {
    coords <- tryCatch(as.data.table(coords), error = function(e) e)
    if (is(coords, "error")) {
      stop("Can't coherce 'coords' to a data.table.",
           "  \nPlease pass a data.table, a SpatVector, or object to data.table or SpatVector")
    }
    coords_poly <- suppressWarnings(vect(coords, geom = c("x", "y"), crs = crs))
  } else {
    isSpatial <- TRUE
    if (is(coords, "SpatVector")) {
      coords_poly <- coords
    } else {
      coords_poly <- tryCatch(vect(coords), error = function(e) e)
      if (is(coords_poly, "error")) {
        stop("Can't coherce 'coords' to a SpatVector.",
             "  \nPlease pass a data.table, a SpatVector, or object to data.table or SpatVector")
      }
    }
  }
  
  coords_out <- crop(coords_poly, cropExt)
  
  if (!nrow(coords_out)) {
    warning("No points left after cropping to cropExt. Do projections match?")
  }
  
  if (!isSpatial) {
    coords_out <- as.data.table(coords_out, geom = "XY")
    coords_out <- coords_out[, .SD, .SDcols = names(coords)]  ## re-order
  }
  
  return(coords_out)
}


#' Make generate extents for gaps used in hold-out samples
#'   for model cross-validation.
#'
#' @param studyarea `SpatExtent` of study area where gaps should be generated.
#'   Defaults to an area in Southern BC. 
#' @param ngaps integer. Number of gaps wanted..
#'
#' @return a `list` of extents of length `ngaps`.
#' 
#' @export
#'
#' @importFrom terra ext
makeGapExtents <- function(studyarea = ext(c(-123, -118, 49, 52)), ngaps = 5L) {
  centre <- c(mean(studyarea[1:2]), mean(studyarea[3:4]))
  range <- c(diff(studyarea[1:2]), diff(studyarea[3:4]))
  gapcentre <- matrix(c(-1,1,1,1,1,-1,-1,-1), 4, byrow=T)
  gapextents <- ext(c(centre[1]+c(-1,1)/8*range[1], centre[2]+c(-1,1)/8*range[2]))
  
  ngaps <- ngaps - 1  ## we have one already
  extragaps <- lapply(seq_len(ngaps), function(i) {
    gap <- ext(c(centre[1]+sort(gapcentre[i,1]*c(1,3)/8)*range[1], centre[2]+sort(gapcentre[i,2]*c(1,3)/8)*range[2]))
    gap
  })
  
  return(append(gapextents, extragaps))
}


#' Crop on disk using `gdalUtilities::ogr2ogr`
#'
#' @param filename character. File name (only) to original spatial file.
#' @param filename2 character. File name (only) to use when saving cropped spatial file. 
#' @param destinationPath folder path to `filename` and where `filename2` will be saved.
#' @param studyArea a SpatVector, or NULL (default), from which the extent to crop
#'   `filename` to will be extracted. If NULL, no cropping/saving happens and `filename` is 
#'   simply loaded.
#'
#' @return a cropped SpatVector that corresponds to `filename2`.
#' @export
GDALcrop <- function(filename, filename2, destinationPath = NULL, studyArea = NULL) {
  if (!requireNamespace("gdalUtilities")) {
    stop("Please install 'gdalUtilities' to use this function")
  }
  
  if (!is.null(studyArea)) {
    if (is.null(destinationPath))
      destinationPath <- "./"
    
    ## get original projection
    origCRS <- crs(vect(file.path(destinationPath, filename)), proj = TRUE)
    ## get ext in original projection
    studyAreaExt <- as.vector(ext(project(studyArea, y = origCRS)))
    
    ## crop on disk
    gdalUtilities::ogr2ogr(src_datasource_name = file.path(destinationPath, filename), 
                           dst_datasource_name = file.path(destinationPath, filename2),
                           clipsrc = c(studyAreaExt["xmin"], studyAreaExt["ymin"], 
                                       studyAreaExt["xmax"], studyAreaExt["ymax"]),
                           overwrite = TRUE)
    filename <- filename2
  } else {
    message("studyArea not provided, not cropping.")
  }
  
  shpObj <- vect(file.path(destinationPath, filename))
  return(shpObj)
}


#' Crop on disk using `gdalUtilities::ogr2ogr`
#'
#' @param filename character. File name (only) to original spatial file.
#' @param filename2 character. File name (only) to use when saving cropped spatial file. 
#' @param destinationPath folder path to `filename` and where `filename2` will be saved.
#' @param studyArea a SpatVector, or NULL (default), from which the extent to crop
#'   `filename` to will be extracted. If NULL, no cropping/saving happens and `filename` is 
#'   simply loaded.
#'
#' @return a cropped SpatVector that corresponds to `filename2`.
#' @export
GDALintersect <- function(filename, filename2, destinationPath = NULL, studyArea = NULL) {
  if (!requireNamespace("gdalUtilities")) 
    stop("Please install 'gdalUtilities'")
  
  if (!is.null(studyArea)) {
    if (is.null(destinationPath))
      destinationPath <- "./"
    
    ## get original projection
    origCRS <- crs(vect(file.path(destinationPath, filename)), proj = TRUE)
    ## get ext in original projection
    studyAreaExt <- as.vector(ext(project(studyArea, y = origCRS)))
    
    ## crop on disk
    gdalUtilities::ogr2ogr(src_datasource_name = file.path(destinationPath, filename), 
                           dst_datasource_name = file.path(destinationPath, filename2),
                           clipsrc = c(studyAreaExt["xmin"], studyAreaExt["ymin"], 
                                       studyAreaExt["xmax"], studyAreaExt["ymax"]),
                           overwrite = TRUE)
    filename <- filename2
  } else {
    message("studyArea not provided, not cropping.")
  }
  
  shpObj <- vect(file.path(destinationPath, filename))
  return(shpObj)
}



