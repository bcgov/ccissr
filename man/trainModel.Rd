% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{trainModel}
\alias{trainModel}
\title{\code{mlr3} model training wrapper}
\usage{
trainModel(learner, task)
}
\arguments{
\item{learner}{a class \code{Learner} object from \code{mlr3}}

\item{task}{a class \code{Task} object from \code{mlr3}, compatible
with \code{learner}}
}
\value{
\code{learner} with the fitted model
}
\description{
Used to enable caching
}
\details{
internally runs \code{learner$train(task)}
}
