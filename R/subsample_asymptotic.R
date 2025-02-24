#' calculate an asymptotic subsample size conditional on the population size of each BGC unit. 
#'
#' This function applies an asymptotic exponential subsampling transformation to a numeric vector.
#' Values below the threshold remain unchanged, while values above it are transformed
#' using an asymptotic function to ensure a smooth transition and prevent values from exceeding
#' the original input.
#'
#' Optionally, the function can generate a plot illustrating the transformation.
#'
#' @param N A numeric vector representing the population size in each class.
#' @param asymptote A numeric value specifying the maximum possible subsampled value. Default is 2000.
#' @param threshold A numeric value specifying the point at which subsampling begins. Default is asymptote/20.
#' @param shape A numeric parameter controlling the rate of subsampling. Smaller values lead to a more gradual transition. Default is 1/asymptote.
#'
#' @return A numeric vector of the same length as `N`, where values above `threshold` are transformed
#'         and rounded to the nearest integer. 
#' 
#' @examples
#' # Generate subsampled values
#' x <- seq(0, 10000, 50)
#' y <- subsample(x)
#' plot(x, y, type = "l", col = "blue", lwd = 2, main = "Subsampling Function", xlab = "Population size", ylab = "Sample size")
#' abline(a = 0, b = 1, col = "gray", lty = 2)
#'
#' @export
subsample_asymptotic <- function(N, asymptote = 2000, threshold = NULL, shape = NULL) {
  if(is.null(threshold)) threshold <- asymptote/20
  if(is.null(shape)) shape <- 1/asymptote
  above_thresh <- N > threshold
  n <- N
  n[above_thresh] <- threshold + (asymptote - threshold) * (1 - exp(-shape * (N[above_thresh] - threshold)))
  n[n>N] <- N[n>N]
  return(round(n))
}
