#' Title Selection of the tuning parameter for the penalized generalized estimating equations
#'
#' @description Select the tuning parameter by EBIC and fit the PTC model using SCAD-based penalized estimating equations.
#' @param Time right censored follow-up time.
#' @param Status censoring indicator, 1 = event and 0 = censoring.
#' @param id cluster identifier.
#' @param Var if TRUE, return standard errors.
#' @param X covariate matrix.
#' @param corstr correlation structure: "independence", "exchangeable", or "AR1".
#' @param stad if TRUE, standardize non-intercept covariates.
#' @param beta_int initial coefficient values.
#' @param nopindex indices of covariates not penalized.
#' @param nlambda number of lambda values.
#' @param lambda.min.ratio ratio between minimum lambda and lambdaMax.
#' @param eps convergence tolerance.
#' @param lambda user-supplied lambda sequence; NULL means construct automatically.
#' @param maxiter maximum number of iterations.
#' @param tol coefficients with absolute value below tol are treated as zero.
#' @param QIC if TRUE, calculate QIC.
#' @export fitfun

fitfun <- function(Time, Status, id, Var, X, corstr, stad, beta_int, nopindex=NULL,
                   nlambda, lambda.min.ratio, eps, lambda=NULL, maxiter, tol, QIC){
  
  N <- length(Time)
  K <- length(unique(id))
  p <- ncol(X)
  stadX <- X
  
  if(length(Status) != N) stop("Time and Status do not have the same length.")
  if(length(id) != N) stop("Time and id do not have the same length.")
  if(nrow(X) != N) stop("nrow(X) must equal the number of observations.")
  if(length(beta_int) != p) stop("length(beta_int) must equal ncol(X).")
  if(any(!Status %in% c(0,1))) stop("Status must contain only 0 and 1.")
  if(sum(Status == 1) == 0) stop("At least one observed event is required.")
  

  X.mean <- X.sd <- NULL
  if(stad && p > 1){
    X.mean <- apply(X[, -1, drop=FALSE], 2, mean)
    X.sd <- apply(X[, -1, drop=FALSE], 2, sd)
    if(any(!is.finite(X.sd)) || any(X.sd == 0))
      stop("At least one non-intercept covariate has zero or non-finite standard deviation.")
    
    Ibeta1 <- rep(0, length(beta_int))
    cumIbe <- 0
    for(lln in 2:length(beta_int)){
      cumIbe <- cumIbe +  mean(X[, lln]) * beta_int[lln]  
      Ibeta1[lln] <- sd(X[, lln]) * beta_int[lln]
    }
    Ibeta1[1] <- beta_int[1] + cumIbe
    beta_int <- Ibeta1
    
    for(i in 2:ncol(X))
      stadX[, i] <- (X[, i] - mean(X[, i])) / sd(X[, i])
  }
  

  
  if(is.null(lambda)){
    KM.fit <- survfit(Surv(Time, Status) ~ 1)
    cure.rate <- min(KM.fit$surv)
    bet.init <- c(ifelse(cure.rate == 0, 0, log(-log(cure.rate))), rep(0, p-1))
    
    optim.fit <- optim(par=bet.init, fn=Loss,
                       control=list(maxit=5e4, fnscale=-1, reltol=1e-3), 
                       Time=Time, Status=Status, X=stadX)
    Bthat <- optim.fit$par
    if(any(!is.finite(Bthat))) stop("Initial optimization produced non-finite coefficients.")
    if(optim.fit$convergence != 0)
      warning("Initial optimization used to construct lambdaMax did not fully converge.")
    
    IS <- lambda_IS(beta=Bthat, Time=Time, Status=Status, X=stadX) 
    

    if(is.null(nopindex)){
      lambda.temp <- abs((diag(IS$I) * Bthat + IS$S)[-1] * Bthat[-1])
    }else{
      lambda.temp <- abs(((diag(IS$I) * Bthat + IS$S)[-1] * Bthat[-1])[-nopindex])
    }
    
    lambdaMax <- max(lambda.temp, na.rm=TRUE)
    if(!is.finite(lambdaMax) || lambdaMax <= 0)
      stop("Unable to construct a positive finite lambdaMax.")
    
    lamgrid <- exp(seq(log(lambdaMax * 1.5),
                       log(lambdaMax * lambda.min.ratio),
                       length.out=nlambda))
  }else{
    lamgrid <- as.double(lambda)
    nlambda <- length(lamgrid)
  }
  
  L_grid <- min(lamgrid)
  R_grid <- max(lamgrid)
  

  EBIC <- rep(Inf, nlambda)
  res <- vector("list", nlambda)
  
  for(s in 1:nlambda){
    lam <- lamgrid[s]
    
    fit_PGEE <- tryCatch(
      pgee(Time, Status, id, stadX, corstr, beta_int, nopindex, lam, eps, maxiter),
      error=function(e) NULL
    )
    
    if(is.null(fit_PGEE)){
      res[[s]] <- list(lam=lam, EBIC=Inf, convergence=FALSE)
      next
    }
    
    beta1 <- fit_PGEE$beta
    lam1 <- fit_PGEE$lam
    if(length(beta1) != p || any(!is.finite(beta1))){
      fit_PGEE$EBIC <- Inf
      fit_PGEE$convergence <- FALSE
      res[[s]] <- fit_PGEE
      next
    }
    
    if(!isTRUE(fit_PGEE$convergence)){
      conv.check <- tryCatch({
        R.check <- mycor(Status, id, stadX, beta1, Time, corstr)
        tmp <- U_H_E_M(Time, Status, id, stadX, beta1,
                       R.check$Qhat, R.check$pphi, lam, nopindex)
        beta.check <- matrix(beta1) + MASS::ginv(tmp$H + K*tmp$E) %*%
          (tmp$U - K*tmp$E %*% matrix(beta1))
        all(is.finite(beta.check)) &&
          all(abs(as.double(beta.check) - beta1) <= eps)
      }, error=function(e) FALSE)
      fit_PGEE$convergence <- isTRUE(conv.check)
    }
    
    if(!isTRUE(fit_PGEE$convergence)){
      fit_PGEE$EBIC <- Inf
      res[[s]] <- fit_PGEE
      next
    }
    
    ebic <- tryCatch(
      EBIC(X, Status, Time, id, beta1, stad, tol),
      error=function(e) list(EBIC=Inf, logl=-Inf, fhat=NULL, Fhat=NULL)
    )
    
    if(!is.finite(ebic$EBIC)){
      fit_PGEE$EBIC <- Inf
      fit_PGEE$logl <- ebic$logl
      fit_PGEE$fhat <- ebic$fhat
      fit_PGEE$Fhat <- ebic$Fhat
      res[[s]] <- fit_PGEE
      next
    }
    
    EBIC[s] <- fit_PGEE$EBIC <- as.double(ebic$EBIC)
    fit_PGEE$logl <- ebic$logl
    fit_PGEE$fhat <- ebic$fhat
    fit_PGEE$Fhat <- ebic$Fhat
    res[[s]] <- fit_PGEE
  }
  
  if(!any(is.finite(EBIC)))
    stop("No lambda value produced a converged fit with finite EBIC.")
  
  index.minEBIC <- which.min(EBIC)
  fit <- res[[index.minEBIC]]
  fit$lam.minEBIC <- fit$lam
  fit$lambda.grid <- lamgrid
  fit$EBIC.grid <- EBIC
  fit$convergence.grid <- vapply(res, function(z) isTRUE(z$convergence), logical(1))
  fit$index.minEBIC <- index.minEBIC
  
  beta1 <- fit$beta
  if(stad && p > 1){
    beta2 <- c(beta1[1] - sum(beta1[-1] * X.mean / X.sd),
               beta1[-1] / X.sd)
  }else{
    beta2 <- beta1
  }
  

  fit$beta.raw <- beta2
  beta.minEBIC <- beta2
  beta2 <- as.matrix(beta2)
  

  R.fi.hat <- mycor(Status, id, X, beta2, Time, corstr)
  fit$Qhat <- R.fi.hat$Qhat
  fit$pphi <- R.fi.hat$pphi
  fit$rho <- R.fi.hat$rho
  

  active <- rep(TRUE, length(beta.minEBIC))
  if(length(beta.minEBIC) > 1){
    for(ss in 2:length(beta.minEBIC)){
      if(abs(beta.minEBIC[ss]) < tol){
        beta.minEBIC[ss] <- 0
        active[ss] <- FALSE
      }
    }
  }
  
  fit$beta <- beta.minEBIC
  fit$beta.minEBIC <- beta.minEBIC
  fit$active <- active
  
  if(Var){
    U.H.E.M.val <- U_H_E_M(Time, Status, id, X, beta2,
                           fit$Qhat, fit$pphi, fit$lam.minEBIC, nopindex)
    H <- U.H.E.M.val$H
    E <- U.H.E.M.val$E
    M <- U.H.E.M.val$M
    
    A <- H + K*E
    sv.A <- svd(A, nu=0, nv=0)$d
    fit$condition.number.full <- if(min(sv.A) > 0) max(sv.A)/min(sv.A) else Inf
    fit$min.singular.value.full <- min(sv.A)
    if(!is.finite(fit$condition.number.full) || fit$condition.number.full > 1e10)
      warning("H + K*E is singular or ill-conditioned; variance estimates may be unstable.")
    
    A.inv <- MASS::ginv(A)
    fit$V.b <- A.inv %*% M %*% t(A.inv)
    

    M1 <- M[active, active, drop=FALSE]
    H1 <- H[active, active, drop=FALSE]
    E1 <- E[active, active, drop=FALSE]
    A1 <- H1 + K*E1
    
    sv.A1 <- svd(A1, nu=0, nv=0)$d
    fit$condition.number.selected <- if(min(sv.A1) > 0) max(sv.A1)/min(sv.A1) else Inf
    fit$min.singular.value.selected <- min(sv.A1)
    if(!is.finite(fit$condition.number.selected) || fit$condition.number.selected > 1e10)
      warning("H1 + K*E1 is singular or ill-conditioned; reported SEs may be unstable.")
    
    A1.inv <- MASS::ginv(A1)
    V1 <- A1.inv %*% M1 %*% t(A1.inv)
    
    bet.var <- rep(NA_real_, length(beta.minEBIC))
    bet.var[active] <- diag(V1)
    tiny.neg <- which(!is.na(bet.var) & bet.var < 0 & bet.var > -1e-12)
    if(length(tiny.neg) > 0) bet.var[tiny.neg] <- 0
    bad.neg <- which(!is.na(bet.var) & bet.var < 0)
    if(length(bad.neg) > 0){
      warning("Negative variance estimates were set to NA.")
      bet.var[bad.neg] <- NA_real_
    }
    fit$V.b1 <- bet.var
  }
  
  if(QIC){
    OM <- matrix(0, ncol(X), ncol(X))
    eta.QIC <- as.vector(X %*% beta2)
    if(any(!is.finite(eta.QIC)) || any(eta.QIC > log(.Machine$double.xmax))){
      fit$QIC.value <- Inf
      warning("QIC could not be evaluated reliably because exp(X beta) overflowed.")
    }else{
      mu.QIC <- exp(eta.QIC)
      for(r in 1:ncol(X)){
        for(s in 1:ncol(X))
          OM[r,s] <- sum(fit$Fhat * mu.QIC * X[,r] * X[,s])
      }
      TrM <- sum(diag(OM %*% fit$V.b))
      fit$QIC.value <- -2*fit$logl + 2*TrM
    }
  }
  
  fit$lam.grad.improper <- fit$lam.minEBIC == L_grid | fit$lam.minEBIC == R_grid
  fit$indicateR <- fit$lam.minEBIC == R_grid
  return(fit)
}


#' @title Loss function for beta
#' @description Loss function for the promotion cure model.
#' @param beta coefficient vector.
#' @param Time follow-up time.
#' @param Status event indicator.
#' @param X covariate matrix.
#' @export Loss

Loss <- function(beta, Time, Status, X){
  eta <- as.vector(X %*% beta)
  if(any(!is.finite(eta)) || any(eta > log(.Machine$double.xmax))) return(-Inf)
  
  mu <- exp(eta)
  N <- length(Time)
  dRbet <- sapply(1:length(Time), function(i) sum(mu * (Time >= Time[i])))
  Rbet.min <- min(dRbet[Status == 1])
  interval <- c(Rbet.min - sum(Status), Rbet.min - 1)
  
  lambet <- tryCatch(
    uniroot(getlambda, interval, tol=.Machine$double.eps^0.75,
            dRbet=dRbet, Status=Status)$root,
    error=function(e) NA_real_
  )
  if(!is.finite(lambet)) return(-Inf)
  
  denom <- dRbet[Status == 1] - lambet
  if(any(!is.finite(denom)) || any(denom <= 0)) return(-Inf)
  

  loss <- sum(eta[Status == 1] + log(N) - log(denom)) - lambet
  if(!is.finite(loss)) return(-Inf)
  return(loss)
}


#' @title Calculate the maximum lambda value
#' @description Calculate quantities used for lambdaMax.
#' @param beta coefficient vector.
#' @param Time follow-up time.
#' @param Status event indicator.
#' @param X covariate matrix.
#' @export lambda_IS

lambda_IS <- function(beta, Time, Status, X){
  Time1 <- Time[Status == 1]
  if(length(Time1) == 0) stop("lambda_IS() requires at least one observed event.")
  
  tau <- max(Time[Status == 1])
  Delta <- 1 * (Time <= tau)
  X1 <- X[Status == 1, , drop=FALSE]
  Kn <- length(Time)
  
  eta <- as.vector(X %*% beta)
  if(any(!is.finite(eta)) || any(eta > log(.Machine$double.xmax)))
    stop("Non-finite or overflowing linear predictor in lambda_IS().")
  mu <- exp(eta)
  
  dRbet <- sapply(1:length(Time), function(i) sum(mu * (Time >= Time[i])))
  Rbet.min <- min(dRbet[Status == 1])
  interval <- c(Rbet.min - sum(Status), Rbet.min - 1)
  lambet <- uniroot(getlambda, interval, tol=.Machine$double.eps^0.75,
                    dRbet=dRbet, Status=Status)$root
  Rbetd <- dRbet[Status == 1]
  
  denom <- Rbetd * (Rbetd - lambet)
  if(any(!is.finite(denom)) || any(denom <= 0))
    stop("Invalid denominator in lambda_IS().")
  
  Dd <- t(sapply(1:length(Time1), function(j)
    apply(X * (mu * (Delta * (Time >= Time1[j]) + 1 - Delta)), 2, mean)))
  
  ch2 <- sum(1/denom)
  ch1 <- apply(Dd/denom, 2, sum)
  if(!is.finite(ch2) || ch2 == 0) stop("Invalid ch2 in lambda_IS().")
  
  ch <- ch1/ch2
  hd <- t(t(Dd) - ch)/Rbetd
  S1 <- apply(X1/Kn - hd, 2, sum)
  S2 <- ch
  S <- S1 - S2
  
  I <- 0
  for(j in 1:length(Time1)){
    I <- I +
      t(t(t(X) - Kn*hd[j,]) * (mu * (Delta*(Time >= Time1[j]) + 1 - Delta))) %*%
      t(t(X/Kn) - hd[j,]) / (Rbetd[j] - lambet)
  }
  
  if(any(!is.finite(I)) || any(!is.finite(S)))
    stop("Non-finite I or S produced by lambda_IS().")
  
  return(list(I=I, S=S))
}
