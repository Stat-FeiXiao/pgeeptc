\name{U_H_E_M}
\alias{U_H_E_M}
\title{The element used in iterative algorithm}
\usage{
U_H_E_M(Status, id, X, beta_new, Qhat, pphi, Fhat, lam, eps=1e-6)
}
\arguments{
\item{Status}{the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
\item{X}{a matrix of covariates that may have effect on the failure times.}
\item{beta_new}{the given values of covariates coefficients.}
\item{Qhat}{the estimation of working correlation matrix.}
\item{pphi}{the estimation of scale parameter.}
\item{Fhat}{the estimation of baseline cumulative hazard function.}
\item{lam}{the specific tuning parameter.}
\item{eps}{a fixed small number used in the iterative algorithm. The default is \code{eps = 1e-6}.}
}
\description{
Generate specific matrices for the iterative algorithm process.
}
