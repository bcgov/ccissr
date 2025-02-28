#' Log-transform climate variables
#'
#' @param dat a `data.table` with columns of climate variables corresponding to 
#'   the selected climate `elements`.
#' @param elements character. Climate elements to search for in `dat`.
#' @param base numeric. Log base.
#' @param add.fields logical. If `TRUE`, the new logged variables are added to `dat` 
#'   (TRUE). Otherwise, original column values are replaced with the logs (FALSE).
#' @param zero_adjust logical. If `TRUE` adjusts zeroes in raw data as:
#'   \eqn{base^{\log_base{x_min} - 1}}.
#'   where \eqn{x_min} is the minimum non-zero, non-NA value.
#'
#' @details
#'   All column names that partially match strings in `elements` will be 
#'   log-transformed.
#' 
#' @return `data.table` 
#' @export
logVars <- function(dat,
                    elements = c("AHM", "DD", "Eref", "FFP", "NFFD", "PAS", "PPT", "SHM", "CMD"),
                    base = exp(1),
                    add.fields = FALSE,
                    zero_adjust = FALSE) {
  
  dat <- copy(dat)
  
  # Fields to operate on (generally these should be ratio (zero-limited) variable)
  logFields <- grep(paste(elements, collapse = "|"), names(dat), value = TRUE)
  dat.log <- dat[, .SD, .SDcols = logFields]
  
  # If specified by the user, give zero values a positive value that is one order of magnitude less than the minimum positive value
  if (zero_adjust) {
    dat.log <- dat.log[, lapply(.SD, function(x) {
      x[x <= 0] <- base^(log(min(x[x > 0], na.rm = TRUE), base = base) - 1)
      return(x)
    })]
  }
  
  # Perform log transformation
  dat.log <- dat.log[, lapply(.SD, function(x) log(x, base = base))]
  
  # Add 
  if(add.fields){
    setnames(dat.log, logFields, paste0(logFields, "_log"))
    dat <- cbind(dat, dat.log)
  } else {
    dat[, (logFields) := Map(x =.SD, xname = logFields, f = function(x, xname) {
      x <- dat.log[[xname]]
      return(x)
    }), .SDcols = logFields]
  }
  return(dat)
}
