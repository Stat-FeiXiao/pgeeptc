\name{BIC}
\alias{BIC}
\title{The Bayesian information criterion for selection of the tuning parameter}
\usage{
BIC(X, Status, id, beta, fhat, Fhat, H, E)
}
\arguments{
\item{X}{a matrix of covariates that may have effect on the failure times.}
\item{Status}{the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
\item{beta}{the given values of covariates coefficients.}
\item{fhat}{the estimation of the baseline probability density function.}
\item{Fhat}{the estimation of the baseline probability density function.}
\item{H}{the first-order differential matrix of likelihood corresponding to \code{beta}.}
\item{E}{the SCAD penalty matrix used in penalized generalized estimation equations.}
}
\description{
Calculate the BIC score for selection of the tuning parameter.
}
