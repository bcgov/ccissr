% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{rmLowSampleBGCs}
\alias{rmLowSampleBGCs}
\title{Remove BGCs with low sample sizes}
\usage{
rmLowSampleBGCs(dat, cutoff = 30)
}
\arguments{
\item{dat}{a \code{data.table} with a "BGC" column and number of rows
being the sampled points.}

\item{cutoff}{minimum number of points necessary to retain a
given BGC.}
}
\value{
\code{dat} without the excluded BGCs and their samples.
}
\description{
Remove BGCs with low sample sizes
}
