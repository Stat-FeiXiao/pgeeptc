#' Title  Selection of the tuning parameter for the penalized generalized estimating equations
#'
#' @description We select the tuning parameters based on the BIC and fit the PTC model using SCAD-based penalized estimating equations.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param Var if it is TRUE, the program returns Std.Error. By default, \code{Var = TRUE}.
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param beta_int the initial value of the covariate coefficients.
#' @param N the degree of Bernstein polynomial used to approximate the baseline distribution function in model.
#' @param tau the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.
#' @param nlambda the number of \code{lambda} values. The default is \code{nlambda = 100}.
#' @param eps a fixed small number used in the iterative algorithm. The default is \code{eps = 1e-6}.
#' @param lambda a user supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence based on \code{nlambda}. Supplying a value of \code{lambda} overrides this.
#' @param maxiter specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.
#' @param tol tolerance for convergence. The default is \code{eps = 1e-3}. Iteration stops once the relative change in deviance is less than \code{tol}.
#'
#' @export fitfun

fitfun <- function(Time,Status, id,Var,X,corstr ,beta_int, N,tau,nlambda,eps,lambda,maxiter, tol){
      
   if(is.null(lambda)){

    step_length <- .01  
    lam_search1 <- .001  
  SK <- 0
  repeat{
  lam_search1 <- lam_search1+step_length*SK
  lam_search <- lam_search1 
  SK <- SK + 1
  beta1 <- pgee(Time,Status, id,X,corstr ,beta_int, N,tau,lam_search,eps , maxiter, tol)$beta
  if(all(abs( beta1[-1] )<tol)) break
  }

  lamgrid <-  exp(seq( log(lam_search1), log(.001*lam_search1),len=nlambda) )
  L_grid <- min(lamgrid)  
  R_grid <- max(lamgrid)  
  }else{
    lam_search1 <-lamgrid <- lambda
    lam_search <-lam_search1 
    nlambda <- length(lambda)
    L_grid <- min(lambda) 
    R_grid <- max(lambda)  
  }


  p<-ncol(X)


  BIC<-rep(0,nlambda)

  for(s in 1:nlambda){
    lam1<- lamgrid[s]
    lam<- lam1
    fit_PGEE <- pgee(Time,Status, id,X,corstr ,beta_int,N,tau,lam,eps , maxiter, tol)
    beta1<-fit_PGEE$beta
    Fhat <- fit_PGEE$Fhat
    fhat <- fit_PGEE$fhat
    pphi <- fit_PGEE$pphi
    H <- fit_PGEE$H
    E <- fit_PGEE$E           
    BIC[s]<- as.double(BIC(X,Status,id,beta1,fhat,Fhat,H,E,tol) )
  }
  lamfin1<-lamgrid[which.min(BIC)]
  lamfin<-lamfin1

 
  fit <- pgee(Time,Status,id,X,corstr, beta_int,N,tau,lamfin,eps , maxiter, tol)

 if(Var){
  H<-fit_PGEE$H
  E<-fit_PGEE$E
  M<-fit_PGEE$M
  fit$V.b<- MASS::ginv(H+K*E)%*%M%*%MASS::ginv(H+K*E)    
 }
 
  

  fit$lamgrid<-lamgrid
  fit$BIC <- BIC
 

  if(lamfin1 == L_grid |  lamfin1 ==R_grid ) {
    fit$indicate <- TRUE
  }else{
    fit$indicate <- FALSE
  }
  if( lamfin1 ==R_grid ) {
    fit$indicateR <- TRUE
  }else{
    fit$indicateR <- FALSE
  }
  return(fit)
}
