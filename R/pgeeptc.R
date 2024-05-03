#============== The main function for variable selection in the PTC model based on PGEE ==============#

#' Title  Variable selection in the marginal semiparametric promotion time cure model based on PGEE
#'
#' @description Variable selection in the marginal semiparametric promotion time cure model based on penalized general estimation equations (GEE) method. 
#' We consider three common correlation structures in this funciton. 
#' The baseline cumulative distribution function in PTC model is appromximated by the Bernstein polynomial.
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring. The expression to the right of the "~" specifies the effect of covariates on the failure time.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param data  a data frame in which to interpret the variables named in the \code{formula}.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param Ibeta the initial value of the covariate coefficients.
#' @param Var if it is TRUE, the program returns Std.Error. By default, \code{Var = TRUE}.
#' @param N the degree of Bernstein polynomial used to approximate the baseline distribution function in model.
#' @param tau the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.
#' @param lambda a user supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence based on \code{nlambda}. Supplying a value of \code{lambda} overrides this.
#' @param boots if it is TRUE, the program returns Std.Error by the bootstrap method. By default, \code{boots = FALSE}.
#' @param nboot the number of bootstrap samples. The default is \code{nboot = 100}.#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.
#' @param QIC if it is TRUE, the program returns the quasilikelihood information criterion score of the fitted model. By default, \code{boots = FALSE}.
#' @param nlambda the number of \code{lambda} values. The default is \code{nlambda = 100}.
#' @param eps a fixed small number used in the iterative algorithm. The default is \code{eps = 1e-6}.
#' @param maxiter specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{esmax} iterations and the estimates will be based on the last iteration. The default \code{itermax = 100}.
#' @param tol tolerance for convergence. The default is \code{eps = 1e-3}. Iteration stops once the relative change in deviance is less than \code{tol}.
#' @return An object of class \code{pgeeptc} is returned. It can be examined by \code{print.pgeeptc()}.
#'
#' @export pgeeptc

pgeeptc<-function(formula, id, data, corstr = "independence", Ibeta=NULL, Var=FALSE, N, tau,lambda=NULL, boots=FALSE, nboot=100, QIC =FALSE, nlambda=100, eps = 1e-6, maxiter = 100, tol = 1e-3){
  call <- match.call()
  uid <- sort(unique(id))
  newid <- rep(0, length(id))
  for (i in 1:length(id)) {
    j <- 1
    repeat {
      if (id[i] != uid[j])
        j <- j + 1 else {
          newid[i] <- j         
          break                 
        }
    }
  }
  data$id <- newid
  data1 <- data[data$id == 1, ]
  for (i in 2:length(uid)) {
    data1 <- rbind(data1, data[data$id == i, ])
  }
  data <- data1
  id <- data$id
  nobs <- length(id)
  K <- length(unique(id))
  nt <- as.vector(table(id))
  maxclsz <-max(nt)
  nobs<-sum(nt)
  mf <- stats::model.frame(formula, data)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  X <- as.matrix(X)
  colnames(X) <- colnames(stats::model.matrix(attr(mf, "terms"), mf))
  beta_name <- colnames(X)
  beta_length <- ncol(X)
  Y <- stats::model.extract(mf, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  Time <- Y[, 1]
  Status <- Y[, 2]

  if(length(id) != length(Time))  stop("Id and Time do not have the same length!")

  if(!(is.double(X)))  X <- as.double(X)
  if(!(is.double(Time)))  Time <- as.double(Time)
  if(!(is.double(id))) id <- as.double(id)
  if(is.null(Ibeta)) Ibeta <- rep(0,dim(X)[2])
  
  p<-ncol(X)
  if(beta_length != p) {stop("Dimension of beta != ncol(X)!")}
  
  if(!(is.double(K)))      K <- as.double(K)
  if(!(is.double(maxclsz)))  maxclsz <- as.double(maxclsz)
  if(!(is.double(nobs)))   nobs <- as.double(nobs)

  corstrs <- c("independence", "exchangeable", "AR1")

  corstrv <- as.integer(match(corstr, corstrs, -1))
  if(corstrv < 1) stop("unknown corstr!")



  fitf <- fitfun(Time,Status, id,Var,X,corstr ,Ibeta, N,tau,nlambda,eps,lambda,maxiter, tol)

  fitf$num_of_clusters <- K
  fitf$max_cluster_size <- maxclsz
  
  fitf$call <- call
  fitf$Var <- Var
  fitf$boots <- boots
  fitf$QIC <- QIC
  fitf$beta_name <- beta_name
  
  if(Var){
    
    fitf$beta_var <- diag(fitf$V.b)
    fitf$beta_sd <-  sqrt(fitf$beta_var)
    fitf$beta_zvalue <- fitf$beta/fitf$beta_sd
    fitf$beta_pvalue <- (1 - stats::pnorm(abs(fitf$beta_zvalue))) * 2
  }

  if(boots){
    
  boots_sd<-boots(formula, id, data, nboot, corstr, Ibeta, lambda, nlambda, eps, maxiter, tol)
  fitf $ boot_rho_sd <- boots_sd$boot_rho_sd
  fitf $ boot_beta_sd <- boots_sd$boot_beta_sd
  fitf$boot_beta_zvalue <- fitf$beta/fitf$boot_beta_sd
  fitf$boot_beta_pvalue <- (1 - stats::pnorm(abs(fitf$boot_beta_zvalue))) * 2

  }
  
  if(QIC){
    Kn <- length(Time)
    t2 <- Time
    t11 <- sort(Time)
    c11 <- Status[order(Time)]
    x111 <- as.matrix(X[order(Time), ])
    tt1 <- unique(t11[c11 == 1])
    
    kk <- length(table(t11[c11 == 1]))
    K <- length(unique(id))
    nt <- as.vector(table(id))
    dd <- as.matrix(table(t11[c11 == 1]))
    OM <- matrix(0,p,p)
    
    for(r in 1:p ){
      for(s in 1:p){
        OM[r,s] <- - sum( - fitf $ Fhat * (exp(X %*% fitf $ beta)*X[,r]*X[,s] )   )
      }
    }
    TrM <- sum(diag( OM%*%(fitf$V.b) ) )
    c11 <- Status[order(Time)]
    x111 <- as.matrix(X[order(Time), ])
    xx11 <- as.matrix(x111[c11 == 1,])
    fhat1 <- fitf $fhat [c11 == 1]
    fun1<- log(exp(xx11%*%fitf $ beta)*fhat1)
    fun2<- -exp(X%*%fitf $ beta)*fitf $ Fhat
    loglikelihood<- sum(fun1)+sum(fun2)
    
    fitf $QIC.value <- -2*loglikelihood+ 2*TrM
  }
  class(fitf) <- "pgeeptc"
  return(fitf)
}

#' Title Print pgeeptc object
#'
#' @param x an object of \code{pgeeptc}.
#' @param ... further arguments to be added in the \code{print.pgeeptc} function.
#'
#' @export print.pgeeptc

print.pgeeptc <- function(x, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nEstimation of Marginal Promotion Time Cure Model Based on PGEE:\n")
  if (x$Var) {
    bt <- array(c(x$beta), c(length(x$beta), 4))
    rownames(bt) <- c(x$beta_name)
    colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
    if(x$boots){
      bt[, 2] <- c(x$boots_beta_sd)
      bt[, 3] <- c(x$boots_beta_zvalue)
      bt[, 4] <- c(x$boots_beta_pvalue)
    }else{
      bt[, 2] <- c(x$beta_sd)
      bt[, 3] <- c(x$beta_zvalue)
      bt[, 4] <- c(x$beta_pvalue)
    }
    
  }
  else {
    bt <- array(c(x$beta), c(length(x$beta), 1))
    rownames(bt) <-  c(x$beta_name)
    colnames(bt) <- "Estimate"
  }
  print(bt)
  cat("\n")
  if(x$QIC){  
    cat("QIC:", x$QIC.value) 
    cat("\n") 
    }
  cat("Number of clusters:", x$num_of_clusters)
  cat("        Maximum cluster size:", x$max_cluster_size)
  cat("\n")
  invisible(x)
}