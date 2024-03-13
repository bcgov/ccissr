#' Remove outliers from data
#'
#' @param dat a `data.table` target data to "clean"
#' @param alpha numeric. The alpha value used to determine the cutoff for outliers.
#' @param IDvars character. Names of columns from which outliers should be excluded.
#' 
#' @details TODO. Parallelizes computations internally with `foreach`.
#'
#' @return `data.table`
#' @seealso [foreach::foreach()]
#' 
#' @importFrom foreach foreach %do%
#' @importFrom stats mahalanobis qchisq cov
#' 
#' @export
removeOutlier <- function(dat, alpha, vars){
  out <- foreach(curr = unique(as.character(dat$BGC)), .combine = rbind) %do% {
    temp <- dat[dat$BGC == curr,]
    md <- tryCatch(mahalanobis(temp[, vars],
                               center = colMeans(temp[, vars]),
                               cov = cov(temp[, vars])), error = function(e) e)
    if (!inherits(md,"error")){
      ctf <- qchisq(1-alpha, df = ncol(temp)-1)
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

#' Remove BGCs with low sample sizes
#'
#' @param dat a `data.table` with a "BGC" column and number of rows
#'   being the sampled points.
#' @param cutoff minimum number of points necessary to retain a 
#'   given BGC.
#'
#' @return `dat` without the excluded BGCs and their samples.
#' @importFrom data.table as.data.table
#' @export
rmLowSampleBGCs <- function(dat, cutoff = 30) {
  dat <- as.data.table(dat)
  BGC_Nums <- dat[,.(Num = .N), by = BGC]
  BGC_good <- dat[!BGC %in% BGC_Nums[Num < cutoff, BGC],] 
  return(BGC_good)
}


#' `mlr3` model training wrapper
#'
#' Used to enable caching
#' 
#' @details internally runs `learner$train(task)`
#' 
#' @param learner a class `Learner` object from `mlr3`
#' @param task a class `Task` object from `mlr3`, compatible 
#'   with `learner`
#'
#' @return `learner` with the fitted model
trainModel <- function(learner, task) {
  learner$train(task)
  
  return(learner)
}