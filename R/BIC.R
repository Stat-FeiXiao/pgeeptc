#' Title The Bayesian information criterion for selection of the tuning parameter
#'
#' @description Calculate the BIC score for selection of the tuning parameter.
#'
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param Status the censoring indicator, 0 = event of interest happens, and 0 = censoring.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param beta the given values of covariates coefficients.
#' @param fhat the estimation of the baseline probability density function.
#' @param Fhat the estimation of the baseline probability density function.
#' @param H the first-order differential matrix of likelihood corresponding to \code{beta}.
#' @param E the SCAD penalty matrix used in penalized generalized estimation equations.
#'
#' @export BIC

BIC<-function(X, Status, id, beta, fhat, Fhat, H, E){

  K1 <- length(unique(id) )

  beta <- as.double(beta)

 fun1<- Status*log(exp(X%*%beta)*fhat)
 fun2<- -exp(X%*%beta)*Fhat
 logl<- sum(fun1)+sum(fun2)
 MATrix <- MASS::ginv(H+K1*E)%*%(H) 
 df <- sum(diag(MATrix)) 

 BIC<-  - logl /(K1)   +  df*log(K1)/(2*K1)     

  return(BIC = BIC)
}

