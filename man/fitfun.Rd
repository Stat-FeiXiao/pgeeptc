\name{fitfun}
\alias{fitfun}
\title{Selection of the tuning parameter for the penalized generalized estimating equations}
\usage{
fitfun(Time,Status, id,Var,X,corstr ,beta_int, N,tau,nlambda,eps,lambda,maxiter, tol)
}
\arguments{
\item{Time}{right censored data which is the follow up time.}
\item{Status}{the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
\item{Var}{if it is TRUE, the program returns Std.Error. By default, \code{Var = TRUE}.}
\item{X}{a matrix of covariates that may have effect on the failure times.}
\item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.}
\item{beta_int}{the initial value of the covariate coefficients.}
\item{N}{the degree of Bernstein polynomial used to approximate the baseline distribution function in model.}
\item{tau}{the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.}
\item{nlambda}{the number of \code{lambda} values. The default is \code{nlambda = 100}.}
\item{eps}{a fixed small number used in the iterative algorithm. The default is \code{eps = 1e-6}.}
\item{lambda}{a user supplied \code{lambda} sequence. Typical usage is to have the program compute its own \code{lambda} sequence based on \code{nlambda}. Supplying a value of \code{lambda} overrides this.}
\item{maxiter}{specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.}
\item{tol}{tolerance for convergence. The default is \code{eps = 1e-3}. Iteration stops once the relative change in deviance is less than \code{tol}.}
}
\description{
We select the tuning parameters based on the BIC and fit the PTC model using SCAD-based penalized estimating equations.
}
