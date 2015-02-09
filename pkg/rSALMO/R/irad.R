#' Irradiation Data Set
#'
#' Example data set for the year 2005 with hourly global irradiation data. 
#'
#' @format Data frame with the following columns:
#' \describe{
#'   \item{time}{time in UTC, the year is 2005}
#'   \item{irad}{global irradiation, hourly sum in J/cm^2. The series contains some missing data.}
#'   \item{irad2}{global irradiation, hourly sum in J/cm^2; missing date were arbitrarily replaced with corresponding data from the next year.}
#'   \item{valid}{logical if the value in irad2 is a measured value from the original data set or not.}
#' }
#'
#' @source Data taken from Deutscher Wetterdienst \url{http://www.dwd.de}, station Dresden Klotzsche (Germany), station number1048.
#'
#' @name irad
#' @docType data
#' @keywords data
#'
#' @examples
#' ## derive daily sums (J/cm^2/d)
#' data(irad)
#' irad$day   <- floor(as.numeric(irad$time)/60/60/24)
#' daily      <- aggregate(list(irad = irad$irad2[-1]), list(day=irad$day[-1]), sum)
#' daily$time <- as.POSIXct((daily$day) * 60*60*24, origin = "1970-01-01 00:00.00 UTC")
#'
#' daily <- data.frame(
#'    time = as.POSIXct(format(daily$time, "%Y-%m-%d")),
#'    irad = daily$irad
#' )
#'
#' plot(daily, type="l", xlab="2005", ylab = "global irradiation (J/cm^2/d)")
NULL
