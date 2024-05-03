\name{pgee}
\alias{pgee}
\title{The penalized generalized estimating equations}
\usage{
pgee(Time, Status, id, X, corstr, beta_int, N, tau, lam, eps , maxiter, tol)
}
\arguments{
\item{Time}{right censored data which is the follow up time.}
\item{Status}{the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
\item{X}{a matrix of covariates that may have effect on the failure times.}
\item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.}
\item{beta_int}{the initial value of the covariate coefficients.}
\item{N}{the degree of Bernstein polynomial used to approximate the baseline distribution function in model.}
\item{tau}{the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.}
\item{lam}{the selected tuning parameter based on the BIC or a single tuning parameter provided by the user.}
\item{eps}{a fixed small number used in the iterative algorithm. The default is \code{eps = 1e-6}.}
\item{maxiter}{specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.}
\item{tol}{tolerance for convergence. The default is \code{eps = 1e-3}. Iteration stops once the relative change in deviance is less than \code{tol}.}
}
\description{
The SCAD-based penalized estimating equations with the specified tuning parameter.
}
