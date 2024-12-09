% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/U_H_E_M.R
\name{U_H_E_M}
\alias{U_H_E_M}
\title{Title The element used in iterative algorithm}
\usage{
U_H_E_M(Time, Status, id, X, beta_new, Qhat, pphi, lam, nopindex)
}
\arguments{
\item{Time}{right censored data which is the follow up time.}

\item{Status}{the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.}

\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}

\item{X}{a matrix of covariates that that may have effect on the failure times.}

\item{beta_new}{the estimation of covariates coefficients.}

\item{Qhat}{the estimation of working correlation matrix.}

\item{pphi}{the estimation of scale parameter.}

\item{lam}{the tunning parameter in SCAD penalty.}

\item{nopindex}{a vector indicating which covariates are not penalized.}
}
\description{
Generate specific matrices for the iterative algorithm process.
}
