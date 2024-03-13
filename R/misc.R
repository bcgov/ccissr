#' Garbage collect a few times
#' 
#' Sometimes running `gc(reset = TRUE)` once doesn't
#' release all unused memory.
#' TODO: move to ccissr package.
#' 
#' @param times integer. Number of times to repeat `gc(reset = TRUE)`
#'
#' @return NULL
#' 
#' @export
.gc <- function(times = 3L) {
  for(i in 1:times) invisible({gc(reset = TRUE)})
}


vertDist <- function(x) {
  (sideLen(x)*3)/2
}

sideLen <- function(x) {
  x/sqrt(3)
}