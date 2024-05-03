#' Title Variance estimation based on the bootstrap method
#'
#' @description The variance estimation of the covariate coefficients and correlation coefficient in the working correlation matrix.
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring. The expression to the right of the "~" specifies the effect of covariates on the failure time.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param data  a data frame in which to interpret the variables named in the \code{formula}.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param beta_int the initial value of the covariate coefficients.
#' @param lambda a user supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence based on \code{nlambda}. Supplying a value of \code{lambda} overrides this.
#' @param nlambda the number of \code{lambda} values. The default is \code{nlambda = 100}.
#' @param eps a fixed small number used in the iterative algorithm. The default is \code{eps = 1e-6}.
#' @param maxiter specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.
#' @param tol tolerance for convergence. The default is \code{eps = 1e-3}. Iteration stops once the relative change in deviance is less than \code{tol}.
#' 
#' @export boots

boots<-function(formula, id, data, nboot, corstr, beta_int, lambda, nlambda, eps, maxiter, tol){
  Bootsample <- nboot

  corstr <- corstr
  K <- length(unique(id))
  BMs <- matrix(0, Bootsample, ncol(X))
  BMnus <-  rep(0, Bootsample)

  for (rrn in 1:Bootsample) {
    repeat{
    bootid<- sample((1:K), replace = TRUE)
    bootdata <- data[id == bootid[1], ]
    bootdata$id <- rep(1, sum(id == bootid[1]))
    for (ll in 2:K) {
      bootdata1 <- data[id == bootid[ll], ]
      bootdata1$id <- rep(ll, sum(id == bootid[ll]))
      bootdata <- rbind(bootdata, bootdata1)
    }
    id_boot <- bootdata$id
    Kn_boot <- length(id_boot)
    K_boot <- length(unique(id_boot))
    n_boot <- as.vector(table(id_boot))
    mf_boot <- model.frame(formula, bootdata)
    X_boot <- model.matrix(attr(mf_boot, "terms"), mf_boot)[, -1]
    X_boot <- as.matrix(X_boot)
    colnames(X_boot) <- colnames(model.matrix(attr(mf_boot, "terms"), mf_boot))[-1]
    Y_boot <- model.extract(mf_boot, "response")
    if (!inherits(Y_boot, "Surv"))
      stop("Response must be a survival object")
    Time_boot <- Y_boot[, 1]
    Status_boot <- Y_boot[, 2]
    
    tryboot <- try(fitfun(Time_boot,Status_boot, id_boot, Var=FALSE, X_boot, corstr ,beta_int, N, tau, nlambda, eps,lambda, maxiter, tol), silent = F)
    if(is(tryboot,"try-error") == FALSE)  break
    }
    esfitboot <- tryboot
    BMs[rrn, ] <- esfitboot$beta
    BMnus[rrn] <- esfitboot$rho
 

  }
  var_beta_boots <- apply(BMs, 2, var)
  sd_beta_boots <- sqrt(var_beta_boots)
  var_rho_boots <- var(BMnus)
  sd_rho_boots <- sqrt( var_rho_boots)
  
  return(sd_beta_boots=sd_beta_boots,sd_rho_boots=sd_rho_boots)

}