#' Title  The penalized generalized estimating equations
#'
#' @description The SCAD-based penalized estimating equations with the specified tuning parameter. 
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param beta_int the initial value of the covariate coefficients.
#' @param nopindex a vector indicating which covariates are not penalized.
#' @param lam the tunning parameter in SCAD penalty.
#' @param eps tolerance for convergence. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param maxiter specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{maxiter} iterations and the estimates will be based on the last iteration.
#'
#' @export pgee

pgee <- function(Time, Status, id, X, corstr, beta_int, nopindex, lam, eps , maxiter){

  K<-length(unique(id))
  beta_new<-beta_int

  iter <- 1
  while(iter <= maxiter) {
  beta_old<-beta_new

  R.fi.hat=mycor(Status,id,X,beta_new,Time,corstr)
  Qhat=R.fi.hat$Qhat
  pphi=R.fi.hat$pphi 
  U.H.E.M.val=U_H_E_M(Time,Status,id,X,beta_new,Qhat,pphi,lam,nopindex)
  U<-U.H.E.M.val$U
  H<-U.H.E.M.val$H
  E<-U.H.E.M.val$E
  M<-U.H.E.M.val$M

  beta_new<-matrix(beta_old)+MASS::ginv(H+K*E)%*%(U-K*E%*%matrix(beta_old))

  iter<-iter+1 

  if (all ( abs(beta_old-beta_new)<= eps ) ) break
  }

  fit <- list()
  
  
  fit$beta <- as.double(beta_new)
  fit$lam <- lam
  fit$convergence <- (iter < maxiter)
  return(fit)
}