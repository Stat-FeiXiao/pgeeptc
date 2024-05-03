\name{mycor}
\alias{mycor}
\title{Estimation of the working correlation matrix and scale parameter}
\usage{
mycor(Status, id, X, beta_new, Fhat, corstr)
}
\arguments{
\item{Status}{the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
\item{X}{a matrix of covariates that may have effect on the failure times.}
\item{beta_new}{the given values of covariates coefficients.}
\item{Fhat}{the estimation of baseline cumulative hazard function.}
\item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.}
}
\description{
Estimate the working correlation matrix and scale parameter, where the structure of working correlation matrix can be specified as \code{independence}, \code{exchangeable} and \code{AR1}.
}
