% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{.gc}
\alias{.gc}
\title{Garbage collect a few times}
\usage{
.gc(times = 3L)
}
\arguments{
\item{times}{integer. Number of times to repeat \code{gc(reset = TRUE)}}
}
\description{
Sometimes running \code{gc(reset = TRUE)} once doesn't
release all unused memory.
TODO: move to ccissr package.
}
