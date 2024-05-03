\name{pgeeptc}
\alias{pgeeptc}
\title{Variable selection in the marginal semiparametric promotion time cure model based on PGEE}
\usage{
pgeeptc(formula, id, data, corstr = "independence", Ibeta=NULL, Var=FALSE, N, tau,lambda=NULL, boots=FALSE, nboot=100, QIC =FALSE, nlambda=100, eps = 1e-6, maxiter = 100, tol = 1e-3){

}
\arguments{
\item{formula}{a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring. The expression to the right of the "~" specifies the effect of covariates on the failure time.}
\item{id}{a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.}
\item{data}{a data frame in which to interpret the variables named in the \code{formula}.}
\item{corstr}{a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.}
\item{Ibeta}{the initial value of the covariate coefficients.}
\item{Var}{if it is TRUE, the program returns Std.Error. By default, \code{Var = TRUE}.}
\item{N}{the degree of Bernstein polynomial used to approximate the baseline distribution function in model.}
\item{tau}{the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.}
\item{lambda}{a user supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence based on \code{nlambda}. Supplying a value of \code{lambda} overrides this.}
\item{boots}{if it is TRUE, the program returns Std.Error by the bootstrap method. By default, \code{boots = FALSE}.}
\item{nboot}{the number of bootstrap samples. The default is \code{nboot = 100}.#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.}
\item{QIC}{if it is TRUE, the program returns the quasilikelihood information criterion score of the fitted model. By default, \code{boots = FALSE}.}
\item{nlambda}{the number of \code{lambda} values. The default is \code{nlambda = 100}.}
\item{eps}{a fixed small number used in the iterative algorithm. The default is \code{eps = 1e-6}.}
\item{maxiter}{specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.}
\item{tol}{tolerance for convergence. The default is \code{eps = 1e-3}. Iteration stops once the relative change in deviance is less than \code{tol}.}
}
\description{
Variable selection in the marginal semiparametric promotion time cure model based on penalized general estimation equations (GEE) method. We consider three common correlation structures in this funciton. The baseline cumulative distribution function in PTC model is appromximated by the Bernstein polynomial.
}
