\name{baseF}
\alias{baseF}
\title{Estimate the baseline cumulative distribution function in the PTC model using Bernstein polynomial}
\usage{
baseF(Time, Status, X, beta, N, tau)
}
\arguments{
\item{Time}{right censored data which is the follow up time.}
\item{Status}{the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.}
\item{X}{a matrix of covariates that may have effect on the failure times.}
\item{beta}{estimation of the covariate coefficients based on the GEE or QIF approach.}
\item{N}{the degree of Bernstein polynomial used to approximate the baseline distribution function in model.}
\item{tau}{the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.}
}
\description{
The estimation of baseline survival function based on the Bernstein polynomial.
}
