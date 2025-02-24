#' Remove outliers from data
#'
#' @param dat a `data.table` target data from which to remove outliers
#' @param alpha numeric. The alpha (significance) value used to determine the cutoff for outliers.
#' @param vars character. Names of columns from which outliers should be excluded.
#' 
#' @details TODO. Parallelizes computations internally with `foreach`.
#'
#' @return A `data.table` containing the non-outlying observations
#' @seealso [foreach::foreach()]
#' 
#' @importFrom foreach foreach %do%
#' @importFrom stats mahalanobis qchisq cov
#' 
#' @export
removeOutlier <- function(dat, alpha, vars){
  out <- foreach(curr = unique(as.character(dat$BGC)), .combine = rbind) %do% {
    temp <- dat[dat$BGC == curr,]
    ## log-transform ratio variables and subset for just the variables that we are analyzing for outliers
    popn <- logVars(temp[,..vars], zero_adjust = TRUE)  
    ## remove variables with non-finite values in the target population (this is an edge case that occurs when the target population has a variable (typically CMD) with only zeroes)
    popn <- popn[, lapply(.SD, function(x) if (all(is.finite(x))) x else NULL)]
    ## z-standardize
    clim.mean <- popn[, lapply(.SD, mean, na.rm = TRUE)]
    clim.sd <- popn[, lapply(.SD, sd, na.rm = TRUE)]
    popn[, (names(popn)) := lapply(names(popn), function(col) {
      (get(col) - unlist(clim.mean)[col]) / unlist(clim.sd)[col]
    })]
    ## M distance
    md <- tryCatch(mahalanobis(popn,
                               center = colMeans(popn),
                               cov = cov(popn)), error = function(e) e)
    if (!inherits(md,"error")){
      ctf <- qchisq(1-alpha, df = ncol(popn)-1)
      outl <- which(md > ctf)
      message(paste("Removing", length(outl), "outliers from", curr, "; "), sep = " ")
      if (length(outl) > 0){
        temp <- temp[-outl,]
      }
    }
    temp
  }
  return(out)
}

## rough Example (need to create sample data)
# dat <- trainData
# before <- trainData[BGC=="CWHvm1",.(BGC, Tmax_sm, PPT_sm)]
# after <- removeOutlier(before, alpha=0.0027, vars = c("Tmax_sm", "PPT_sm"))
# 
# plot(logVars(before[, .(Tmax_sm, PPT_sm)]))
# points(logVars(after[, .(Tmax_sm, PPT_sm)]), pch=16, col="red")
