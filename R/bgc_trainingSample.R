#' Generate a Balanced Training Sample from BGC Units
#' 
#' This function extracts a training sample from a DEM and BGC raster dataset,
#' balancing sample sizes based on different schemes and optionally removing outliers.
#' 
#' @param dem A `SpatRaster` of the digital elevation model.
#' @param bgcs A `SpatRaster` of biogeoclimatic (BGC) units, with resolution and origin matching DEM.
#' @param bgcs_info Optional. A `data.table` containing additional BGC metadata. 
#' minimum fields are "BGC" and "Source", which must indicate non-BC units compatible with grep("USA_|AB_", bgcs_info$Source)
#' @param scheme Sampling scheme. One of `"asymptotic"` or `"squareRoot"`.
#' @param squareRoot.multiplier Numeric. Multiplier for square-root sampling (default: `10`).
#' @param removeOutliers Logical. Whether to remove climatic outliers (default: `FALSE`). The outlier metric is Mahalanobis distance.
#' @param alpha.BC Numeric. Significance level for outlier removal in BC BGC units (default of `0.0027` corresponds to 3-sigma outliers).
#' @param alpha.nonBC Numeric. Significance level for outlier removal in non-BC BGC units (default of `0.05` corresponds to 2-sigma outliers).
#' @param removeOutliers.vars Character vector. Climatic variables to check for outliers.
#' @param climaticVariance Logical. Whether to adjust sample size based on climatic variance.
#' @param climaticVariance.var Character. Variable used for climatic variance scaling. Default "MAT" is mean annual temperature.
#' @param bgcs.remove Character vector. BGC units to exclude from sampling.
#' @param plotDiagnostics Logical. Whether to generate diagnostic plots (default: `TRUE`).
#' @param plot.dir character. The directory to save the plot generated when `plotDiagnostics = TRUE`.
#' @param plot.name character. The filename for the plot generated when `plotDiagnostics = TRUE`. Do not include the file extension (.png) as this is hardcoded. 
#' @param ... Additional arguments (currently unused).
#' 
#' @return A `data.table` containing the subsampled points with columns `id`, `lon`, `lat`, `elev`, and `BGC`.
#' 
#' @importFrom data.table as.data.table setnames .SD .N := rbind
#' @importFrom terra cells
#' @importFrom climr downscale
#' 
#' @details 
#' The function follows these steps:
#' 
#' 1. Extracts point data from the DEM and BGC rasters.
#' 2. Optionally removes specified BGC units (`bgcs.remove`).
#' 3. Balances sample size across BGC units using one of two subsampling schemes:
#'    - **Asymptotic**: Uses an exponential function of population size towards a specified asymptote (default 4000). See `ccissr::subsample()`.
#'    - **SquareRoot**: Square root of population size, multiplied by a specified factor. The default multiplier of 10 would yield a 1:1 subsample for N <= 100. 
#' 4. Optionally adjusts sample size based on climatic variance (`climaticVariance`).
#' 5. Optionally removes outliers from the sample based on climate variables.
#' 6. Generates diagnostic plots if `plotDiagnostics = TRUE`.
#' 
#' @examples
#' \dontrun{
#' library(terra)
#' library(data.table)
#' 
#' dem <- rast("path/to/dem.tif")
#' bgcs <- rast("path/to/bgcs.tif")
#' 
#' sample <- bgc_trainingSample(dem, bgcs, scheme = "asymptotic")
#' head(sample)
#' }
#' 
#' @export
bgc_trainingSample <- function(dem, bgcs, bgcs_info = NULL,
                               scheme = c("asymptotic", "squareRoot"), 
                               squareRoot.multiplier = 10, 
                               removeOutliers = FALSE, alpha.BC = .0027, alpha.nonBC = .05, 
                               removeOutliers.vars = c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt"), 
                               climaticVariance = FALSE, climaticVariance.var = "MAT",
                               bgcs.remove = NULL,
                               plotDiagnostics = TRUE,
                               plot.dir = NULL, 
                               plot.name = "bgc_trainingSample_plot", 
                               ...){
  
  ## -------------------------------------------------
  ## make the points table
  points <- as.data.table(dem, cells=T, xy=T)
  colnames(points) <- c("id", "lon", "lat", "elev")
  points <- points[,c(2,3,4,1)] #restructure for climr input
  # values(X)[points$id] <- points$el ; plot(X)
  
  # extract BGC values and add to the points table
  points.bgc <- as.data.table(bgcs, cells=T, xy=T)
  colnames(points.bgc) <- c("id", "lon", "lat", "BGC")
  points[, BGC := points.bgc[.SD, on = "id", BGC]]
  points <- points[!is.na(points$BGC), ]
  points[, BGC := factor(BGC)]
  
  # values(X) <- NA; values(X)[points$id] <- points$BGC ; plot(X)
  
  ## -------------------------------------------------
  ## remove specified units
  
  if(!is.null(bgcs.remove)){
    valid_bgcs_remove <- bgcs.remove[bgcs.remove %in% points$BGC]
    
    if (length(valid_bgcs_remove) > 0) {
      points <- points[!(BGC %in% valid_bgcs_remove)]
    }    
    
    points[, BGC := factor(BGC)]
  }
  
  ## -------------------------------------------------
  ## initial sample balancing scheme. 
  scheme <- match.arg(scheme, choices = c("asymptotic", "squareRoot"))
  
  if (scheme == "asymptotic") {

    ## Option 1 - get an asymptotic subsample from each BGC unit's grid points based on population size
    pointcount1 <-  table(points$BGC)
    samplesize1 <- subsample_asymptotic(pointcount1, ...)
    points_subsample <- points[, .SD[sample(.N, min(ifelse(BGC %in% names(samplesize1), samplesize1[BGC], .N), .N))], by = BGC]
    dim(points_subsample)
    
    print("Using asymptotic scheme")
    
  } else if (scheme == "squareRoot") {

    ## Option 2 - get a square-root subsample from each BGC unit's grid points based on population size 
    pointcount1 <-  table(points$BGC)
    samplesize1 <- pointcount1^0.5*squareRoot.multiplier #square root multiplied by 10 (1:1 for N=100)
    points_subsample <- points[, .SD[sample(.N, min(ifelse(BGC %in% names(samplesize1), samplesize1[BGC], .N), .N))], by = BGC]
    dim(points_subsample)
    
    print("Using square root scheme")
  }
  
  ## ----------------------------------------------------------------------------
  ## modify sample size based on each BGC unit's climatic variance 
  
  if (climaticVariance) {
    
    # Pull climr data for the subsample:
    clim <- downscale(
      xyz = points_subsample,
      which_refmap = "refmap_climr",
      return_refperiod = TRUE, # Also return the 1961-1990 normals period.
      vars = climaticVariance.var
    )
    clim <- merge(points_subsample, clim, by = "id", all = FALSE)
    dim(clim[id %in% points_subsample$id])
    
    bgc_iqr <- aggregate(clim, by=list(clim$BGC), FUN = IQR, na.rm=T)
    
    samplesize <- table(points_subsample$BGC)*bgc_iqr[,climaticVariance.var] #BGC subsample size multiplied by climate IQR
    points_subsample <- points[, .SD[sample(.N, min(ifelse(BGC %in% names(samplesize), samplesize[BGC], .N), .N))], by = BGC]
    dim(points_subsample)

  }
  
  # ## ----------------------------------------------------------------------------
  # ## get a subsample from each BGC unit's grid points based on spatial aggregation (clumpiness) of unit
  # ## not yet implemented; still half-baked
  # 
  # # get a metric of percent of cells that have same adjacent class
  # library(landscapemetrics) # for pladj
  # pladj <- lsm_c_pladj(bgcs)  # Per-class adjacency metric
  # class_mapping <- cats(bgcs)[[1]]  # Extract numeric-class mapping
  # pladj <- merge(pladj, class_mapping, by.x = "class", by.y = "value", all.x = TRUE)
  # colnames(pladj)[colnames(pladj) == "category"] <- "BGC"
  # pladj <- pladj[which(pladj$BGC %in% points_subsample$BGC),]
  # print(pladj[order(pladj$value),])
  # hist(pladj$value)
  # 
  # # downsample in reverse proportion to pladj
  # pointcount <-  table(points_subsample$BGC)
  # pladj$BGC[-which(pladj$BGC %in% names(pointcount))]
  # names(pointcount)[-which(names(pointcount) %in% pladj$BGC)]
  # pladj_factor <- 1-((pladj$value[match(names(pointcount), pladj$BGC)]/100)^2-0.6) # this could be more elegant
  # pladj_factor[pladj_factor>1] <- 1 
  # pladj_factor[names(pointcount)=="ICHun"] <- 1
  # hist(pladj_factor)
  # samplesize <- round(pointcount*pladj_factor) #
  # points_subsample <- points_subsample[, .SD[sample(.N, min(ifelse(BGC %in% names(samplesize5), samplesize5[BGC], .N), .N))], by = BGC]
  # dim(points_subsample)
  
  ## -------------------------------------------------
  ## Remove outliers 
  
  if (removeOutliers) {
    
    # Pull climr data for the subsample:
    clim <- downscale(
      xyz = points_subsample,
      which_refmap = "refmap_climr",
      return_refperiod = TRUE, # Also return the 1961-1990 normals period.
      vars = c(removeOutliers.vars, "MAT")
    )
    clim <- merge(points_subsample, clim, by = "id", all = FALSE)
    dim(clim[id %in% points_subsample$id])
    bgc_iqr <- aggregate(clim, by=list(clim$BGC), FUN = IQR, na.rm=T)
    
    # Identify BC and non-BC BGC units
    bgcs_nonBC <- bgcs_info[grep("USA_|AB_", bgcs_info$Source), BGC]
    bgcs_BC <- setdiff(unique(points_subsample$BGC), bgcs_nonBC)  # Assuming BGC units that are not in bgcs_nonBC are BC units
    
    # Subset points_subsample for BC and non-BC BGC units
    points_BC <- points_subsample[BGC %in% bgcs_BC, ]
    points_nonBC <- points_subsample[BGC %in% bgcs_nonBC, ]
    
    # Apply outlier removal separately for BC and non-BC
    points_BC <- removeOutlier(clim[id %in% points_BC$id], alpha = alpha.BC, vars = removeOutliers.vars)
    points_nonBC <- removeOutlier(clim[id %in% points_nonBC$id], alpha = alpha.nonBC, vars = removeOutliers.vars)
    
    # Combine the BC and non-BC data back together
    points_subsample <- rbind(points_BC, points_nonBC)
    points_subsample <- points_subsample[, 1:5]  # Keep the first five columns
    
    dim(points_subsample)    
  }
  
  ## -------------------------------------------------
  ## Plots
  
  if(plotDiagnostics) {
    # visualize
    png(filename=paste0(plot.dir, "/", plot.name, ".png",sep="."), type="cairo", units="in", width=6.5, height=6.5, pointsize=10, res=300)
    par(mfrow=c(2,2), mar = c(4,3,1,1), mgp=c(2,0.25,0), tck=-0.01)
    plot(sort((pointcount1)), type = "l", ylab="Number of points", xaxt = "n")
    axis(1, at=seq(1,length(pointcount1)), labels = names(sort((pointcount1))), las=2)
    lines(sort((samplesize1)), type = "l", lty=2)
    legend("topleft", bty="n", legend = c("BGC grid point count", "sample size"), lty=c(1,2))
    box()
    
    plot(sort(log10(pointcount1)), ylim=range(log10(c(samplesize1, pointcount1))), type = "l", ylab="Number of points", yaxt = "n", xaxt = "n")
    axis(1, at=seq(1,length(pointcount1)), labels = names(sort((pointcount1))), las=2)
    axis(2, at=seq(1,9), labels = 10^seq(1,9), las=2)
    lines(sort(log10(samplesize1)), type = "l", lty=2)
    legend("bottomright", bty="n", legend = c("BGC grid point count", "sample size"), lty=c(1,2))
    box()
    
    # plot of BGC subsample size vs BGC area 
    x <- log10(as.vector(table(points$BGC[])))
    y <- log2(as.vector(table(points_subsample$BGC)))
    plot(x, y, col = "white",
         xaxt = "n", yaxt = "n",
         ylab = "sample size of BGC Unit",
         xlab = paste0("BGC unit area (number of cells)"),
    )
    axis(1, at = seq(1,20), labels = round(10^seq(1,20)))
    axis(2, at = 1:99, labels = 2^(1:99), las=2)
    text(x, y, labels = bgc_iqr$Group.1, cex = 0.5)
    
    # plot of BGC subsample size vs climatic variance 
    bgc_iqr <- aggregate(clim[id %in% points_subsample$id], by=list(clim[id %in% points_subsample$id, BGC]), FUN = IQR, na.rm=T)
    x <- log2(bgc_iqr[,climaticVariance.var])
    y <- log2(as.vector(table(points_subsample$BGC)))
    plot(x, y, col = "white",
         xaxt = "n", yaxt = "n",
         ylab = "sample size of BGC Unit",
         xlab = paste0("IQR of ", climaticVariance.var),
    )
    axis(2, at = seq(1,20), labels = round(2^seq(1,20)), las=2)
    axis(1, at = -3:3, labels = 2^(-3:3))
    text(x, y, labels = bgc_iqr$Group.1, cex = 0.5)
    
    dev.off()
    par(mfrow=c(1,1))
  }
  
  return(points_subsample)
  
} #end of function

