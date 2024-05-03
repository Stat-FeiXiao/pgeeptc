#' Title  The penalized generalized estimating equations
#'
#' @description The SCAD-based penalized estimating equations with the specified tuning parameter. 
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param beta_int the initial value of the covariate coefficients.
#' @param N the degree of Bernstein polynomial used to approximate the baseline distribution function in model.
#' @param tau the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.
#' @param lam the selected tuning parameter based on the BIC or a single tuning parameter provided by the user.
#' @param eps a fixed small number used in the iterative algorithm. The default is \code{eps = 1e-6}.
#' @param maxiter specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.
#' @param tol tolerance for convergence. The default is \code{eps = 1e-3}. Iteration stops once the relative change in deviance is less than \code{tol}.
#'
#' @export pgee

pgee <- function(Time, Status, id, X, corstr, beta_int, N, tau, lam, eps , maxiter, tol){

      iter<-0
      K<-length(unique(id))
      beta_new<-beta_int
      Fhat<- baseF(Time, Status, X, beta_new,N,tau)$Fhat
      R.fi.hat=mycor(Status,id,X,beta_new,Fhat,corstr)
      Qhat=R.fi.hat$Qhat
      pphi=R.fi.hat$pphi

      U.H.E.M.val=U_H_E_M(Status,id,X,beta_new,Qhat,pphi,Fhat,lam,eps)
      U<-U.H.E.M.val$U
      H<-U.H.E.M.val$H
      E<-U.H.E.M.val$E
      M<-U.H.E.M.val$M

      while(iter < maxiter) {

        beta_old<-beta_new

        beta_new<-matrix(beta_old)+MASS::ginv(H+K*E)%*%(U-K*E%*%matrix(beta_old))

        Fhat <- baseF(Time, Status, X, beta_new,N,tau)$Fhat
        R.fi.hat=mycor(Status,id, X,beta_new,Fhat,corstr)
        Qhat=R.fi.hat$Qhat
        pphi=R.fi.hat$pphi

        U.H.E.M.val=U_H_E_M(Status,id, X,beta_new,Qhat,pphi,Fhat,lam,eps)
        U<-U.H.E.M.val$U
        H<-U.H.E.M.val$H
        E<-U.H.E.M.val$E
        M<-U.H.E.M.val$M

        iter<-iter+1
        if (all ( abs(beta_old-beta_new)<= tol ) ) break
      }
    fit <- list()

    beta_new <- as.matrix(beta_new)
    BaseF <- baseF(Time, Status, X, beta_new,N,tau)
    Fhat <- BaseF$Fhat
    fhat <- BaseF$fhat
    gam <- BaseF$gam
    
    R.fi.hat=mycor(Status,id, X,beta_new,Fhat,corstr)
    Qhat=R.fi.hat$Qhat
    pphi=R.fi.hat$pphi
    rho=R.fi.hat$rho

    U.H.E.M.val=U_H_E_M(Status,id, X,beta_new,Qhat,pphi,Fhat,lam,eps)
    U<-U.H.E.M.val$U
    H<-U.H.E.M.val$H
    E<-U.H.E.M.val$E
    M<-U.H.E.M.val$M
 
   
    fit$beta <- as.double(beta_new)
    fit$pphi <- pphi
    fit$rho <-  rho
    fit$lam.value<-unique(lam)
    fit$lam<-lam
    fit$Fhat <-  Fhat
    fit$fhat <-  fhat
    fit$gam <- gam
    fit$H <- H
    fit$E <- E
    fit$M <- M
    return(fit)
}
