\name{WeekDayDummy-methods}
\docType{methods}
\alias{WeekDayDummy}
\alias{WeekDayDummy-methods}
\alias{WeekDayDummy,ANY-method}
\alias{WeekDayDummy,character-method}
\title{function: Create Dummy Day-of-Week Variable}
\description{
Helper function to create a dummy, day of the week variable given a set of dates.
}
\usage{
WeekDayDummy(Dates, date.format, weekday = "Monday")
}
\arguments{
\item{Dates}{ A character vector of dates .}
\item{date.format}{ The format of the dates e.g. \dQuote{\%Y-\%m-\%d"} .}
\item{weekday}{ Character string indicating day of week.}
}
\value{
A numeric vector of 0s and 1s (date-dummy variable)/
}
\author{Alexios Ghalanos}
\examples{
data(sp500ret)
Dates=rownames(sp500ret)
# create Monday dummy
monday=WeekDayDummy(Dates, date.format="\%Y-\%m-\%d", weekday = "Monday")
}
\keyword{methods}
