\name{boots}
\alias{boots}
\title{Variance estimation based on the bootstrap method}
\usage{
boots(formula, id, data, nboot, corstr, beta_int, lambda, nlambda, eps, maxiter, tol)
}
\arguments{
\item{formula}{a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring. The expression to the right of the "~" specifies the effect of covariates on the failure time.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
\item{data}{a data frame in which to interpret the variables named in the \code{formula}.}
\item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.}
\item{beta_int}{the initial value of the covariate coefficients.}
\item{lambda}{a user supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence based on \code{nlambda}. Supplying a value of \code{lambda} overrides this.}
\item{nlambda}{the number of \code{lambda} values. The default is \code{nlambda = 100}.}
\item{eps}{a fixed small number used in the iterative algorithm. The default is \code{eps = 1e-6}.}
\item{maxiter}{specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.}
\item{tol}{tolerance for convergence. The default is \code{eps = 1e-3}. Iteration stops once the relative change in deviance is less than \code{tol}.}
}
\description{
The variance estimation of the covariate coefficients and correlation coefficient in the working correlation matrix.
}
