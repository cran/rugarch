\name{ForwardDates-methods}
\docType{methods}
\alias{ForwardDates}
\alias{ForwardDates-method}
\alias{ForwardDates,ANY-method}
\alias{ForwardDates,character-method}
\title{function: Generate Future Dates}
\description{
Given a starting date, this helper function generates a set of future dates
(excl.weekends) for use in forecasting and simulation when using external
regressors.
}
\usage{
ForwardDates(Dates, n.ahead, date.format, periodicity = "days")
}
\arguments{
\item{Dates}{A character vector of dates. The last date is used as the starting
date for the forward date creation.}
\item{n.ahead}{ The number of dates to generate forward.}
\item{date.format}{ The format of the dates e.g. \dQuote{\%Y-\%m-\%d"} .}
\item{periodicity}{ Currently only days is supported.}
}
\value{
A POSIXct vector of future dates.
}
\note{
This is a helper function particularly useful when used with the weekday dummy
variable for simulation and forecasting in light of weekday dummy external 
regressors in the mean or variance equation. For example, if fitting a GARCH
model with a "Monday" dummy variable in the mean equation, then for simulation 
or forecasting, one needs a set of forward deterministic dummy variables for the
Mondays going forward.
}
\examples{
\dontrun{
data(sp500ret)
Dates = rownames(sp500ret)
# generate the 100 forward non-weekend days
fwd100 = ForwardDates(Dates, n.ahead=100, date.format = "\%Y-\%m-\%d", 
periodicity = "days")
# create a dummy vector for those forward days which are Mondays
fwdMonday = WeekDayDummy(as.character(fwd100), date.format = "\%Y-\%m-\%d", 
weekday = "Monday")
}}
\author{Alexios Ghalanos}
\keyword{methods}
