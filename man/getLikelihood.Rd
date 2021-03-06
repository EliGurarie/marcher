% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getLikelihood.R
\name{getLikelihood}
\alias{getLikelihood}
\alias{getLikelihood.res}
\alias{getAIC.nls}
\title{Estimate likelihoods and AICs}
\usage{
getLikelihood(
  p,
  T,
  X,
  Y,
  model = c("mouf", "mou", "mwn", "ouf", "ou", "wn")[1]
)
}
\arguments{
\item{p}{initial parameters: [tau.z, tau.v, t1, t2, x1, x2, y1, y2]}

\item{T, X, Y}{time,x and y coordinates}

\item{model}{"wn", "ou", "ouf", "mou" or "mouf",  - whether or not to estimate tau.v}
}
\description{
Estimate likelihoods and AIC for several possible migration models.
}
