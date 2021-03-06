% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeform.R
\name{makeform}
\alias{makeform}
\title{Construct a model formula for use by MCMCfit()}
\usage{
makeform(zdata, interactions = NULL, includetwofac = T)
}
\arguments{
\item{zdata}{input data matrix}

\item{interactions}{a \eqn{2 \times m} matrix of two-factor interactions}

\item{includetwofac}{If \code{T} include all two factor \emph{excluding} those in \code{interactions}.  If \code{F} then include all main effects but only those two factor interactions in \code{interactions}}
}
\description{
Given the data frame \code{zdata} in \code{Rcapture} form, 
   creates a formula with the interactions given.
}
\details{
This routine is called by the Markov Chain Monte Carlo routine \code{\link{MCMCfit}}.  The routine \code{\link{MSEfit}} uses a formula with 
 numbered variables \eqn{c1, c2, \ldots } rather than using the names of the variables.
}
