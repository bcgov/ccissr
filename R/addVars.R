
#' Add extra climate variables
#'
#' @param dat a `data.table` with columns:
#'    * PPT05, PPT06, PPT07, PPT08, PPT09 (May, ..., September precip), 
#'    PPT_at (autumn precip), PPT_wt (winter precip)
#'    * CMD07 (July climate moisture deficit), CMD (annual CMD)
#'    * DD_0_at (autumn degree-days below 0 deg), DD_0_wt (winter degree-days below 0 deg)
#'    
#' @details This function calculates more climate variables derived from those
#'   output by `[climr::downscale()]`. Presently it adds the following:
#'   -   May to June precip.: \eqn{PPT_MJ = PPT05 + PPT06}
#'   -   July to September precip \eqn{PPT_JAS = PPT07 + PPT08 + PPT09}
#'   -   Precipitation during vegetation dormant period: \eqn{PPT.dormant = PPT_at + PPT_wt}
#'   -   CMD deficit \eqn{CMD.def = 500 - PPT.dormant} (bounded to 0)
#'   -   \eqn{CMDMax = CMD07}
#'   -   \eqn{CMD.total = CMD.def + CMD}
#'   -   \eqn{DD_delayed = ((DD_0_at + DD_0_wt)*0.0238) - 1.8386} bounded to 0)
#'  
#' @return `dat` with all of the above added climate variables.
#' @export
addVars <- function(dat) {
  dat[, PPT_MJ := PPT05 + PPT06]
  dat[, PPT_JAS := PPT07 + PPT08 + PPT09]
  dat[, PPT.dormant := PPT_at + PPT_wt]
  dat[, CMD.def := pmax(0, 500 - PPT.dormant)]
  dat[, CMDMax := CMD07]   ## TODO: THIS IS NOT NECESSARILY CMD MAX
  dat[, CMD.total := CMD.def + CMD]
  dat[, DD_delayed := pmax(0, ((DD_0_at + DD_0_wt)*0.0238) - 1.8386)]
}