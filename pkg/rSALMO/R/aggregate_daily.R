#' Calculate Daily Sums or Averages
#'
#' Calculate daily sums or averages, e.g. from from hourly irradiation data.
#'
#' @param time vector of POSIXct times
#' @param y corresponding vector of data (e.g. hourly)
#' @param basedate reference date where the series starts from, e.g.,
#'   "1970-01-01 00:00.00 UTC"
#' @param FUN function that is applied to the data of each day
#' @return  Data frame with the following columns:
#' \describe{
#'   \item{day}{POSIXct time of the days}
#'   \item{y}{vector of aggregated values}
#' }
#'
#'
#' @export aggregate_daily
#'
#' @examples
#'
#' ## global hourly irradiation data
#' data(irad)
#' daily <- aggregate_daily(irad$time, irad$irad2, basedate= "1970-01-01 00:00.00 UTC")
#'

aggregate_daily <- function (time, y, basedate, FUN = sum) {

  ## derive daily sums (J/cm^2/d)
  ## derive daily sums (J/cm^2/d)
  day   <- floor(as.numeric(time)/60/60/24)
  daily <- aggregate(list(y = y), list(day = day), FUN)

  ## convert time format
  daily$day  <- as.POSIXct(trunc(daily$day) * 60*60*24, origin = "1970-01-01 00:00.00 UTC")
  ## day only; discard hours
  daily$day <- as.POSIXct(format(daily$day, "%Y-%m-%d"))

  return(daily)
}


