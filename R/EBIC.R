#' Title The extended Bayesian information criterion for selection of the tuning parameter
#'
#' @description Calculate the EBIC score for selection of the tuning parameter.
#'
#' @param X a matrix of covariates that may have an effect on the failure times.
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param Time right censored data which is the follow up time.
#' @param id a vector which identifies the clusters. The length of \code{id}
#' should be the same as the number of observations.
#' @param beta the estimated covariate coefficients.
#' @param stad if it is TRUE, all the covariates in \code{X}, except the first
#' column (intercept), are standardized.
#' @param threshold estimated coefficients below \code{threshold} in absolute
#' value are considered as zero.
#'
#' @export EBIC

EBIC <- function(X, Status, Time, id, beta, stad, threshold){
  

  if(length(Time) != length(Status)){
    stop("Time and Status do not have the same length.")
  }
  
  if(length(Time) != length(id)){
    stop("Time and id do not have the same length.")
  }
  
  if(nrow(X) != length(Time)){
    stop("nrow(X) must equal the number of observations.")
  }
  
  if(length(beta) != ncol(X)){
    stop("length(beta) must equal ncol(X).")
  }
  
  if(any(!Status %in% c(0, 1))){
    stop("Status must contain only 0 and 1.")
  }
  
  if(sum(Status == 1) == 0){
    stop("At least one observed event is required to calculate EBIC.")
  }
  

  K1 <- length(unique(id))
  n <- as.vector(table(id))
  N <- length(id)
  

  event.time <- sort(unique(Time[Status == 1]))
  kk <- length(event.time)
  Repf <- as.numeric(table(factor(Time[Status == 1], levels = event.time)))
  
  beta <- as.double(beta)
  stadX <- X
  
  

  if(stad){
    

    if(ncol(X) > 1){
      
      X.sd <- apply(X[, -1, drop = FALSE], 2, sd)
      

      if(any(!is.finite(X.sd)) || any(X.sd == 0)){
        stop("At least one non-intercept covariate has zero or non-finite standard deviation.")
      }
      
      for(i in 2:ncol(X)){
        stadX[, i] <- (X[, i] - mean(X[, i])) / sd(X[, i])
      }
    }
  }
  
  

  XI <- stadX[Status == 1, , drop = FALSE]
  
  

  eta <- as.vector(stadX %*% beta)
  mu <- exp(eta)
  

  if(any(!is.finite(mu))){
    return(list(
      EBIC = Inf,
      logl = -Inf,
      fhat = rep(NA_real_, sum(Status == 1)),
      Fhat = rep(NA_real_, N)
    ))
  }
  

  dRbet <- sapply(1:N, function(i){
    sum(mu * (Time >= Time[i]))
  })
  
  Rbet.min <- min(dRbet[Status == 1])
  interval <- c(Rbet.min - sum(Status), Rbet.min - 1)

  lambet <- uniroot(
    getlambda,
    interval,
    tol = .Machine$double.eps^0.75,
    dRbet = dRbet,
    Status = Status
  )$root
  
  dRbet.event <- dRbet[Status == 1]
  denom <- dRbet.event - lambet
  

  if(any(!is.finite(denom)) || any(denom <= 0)){
    return(list(
      EBIC = Inf,
      logl = -Inf,
      fhat = rep(NA_real_, sum(Status == 1)),
      Fhat = rep(NA_real_, N)
    ))
  }
  
  

  Time.event <- Time[Status == 1]
  
  Fhat <- sapply(Time, function(tmj){
    sum((Time.event <= tmj) / denom)
  })
  
  

  DeltaF <- sapply(event.time, function(tmj){
    sum((Time.event <= tmj) / denom)
  })
  
  if(kk == 1){
    fhat <- DeltaF[1]
  }else{
    fhat <- c(DeltaF[1], DeltaF[2:kk] - DeltaF[1:(kk - 1)])
  }
  

  fhat <- fhat[match(Time.event, event.time)]
  
  if(any(!is.finite(fhat)) || any(fhat <= 0)){
    return(list(
      EBIC = Inf,
      logl = -Inf,
      fhat = fhat,
      Fhat = Fhat
    ))
  }
  
  

  fun1 <- as.vector(XI %*% beta) + log(fhat)
  fun2 <- -mu * Fhat
  logl <- sum(fun1) + sum(fun2)
  
  

  if(stad){
    
    if(ncol(X) > 1){
      X.mean <- apply(X[, -1, drop = FALSE], 2, mean)
      X.sd   <- apply(X[, -1, drop = FALSE], 2, sd)
      
      beta2 <- c(
        beta[1] - sum(beta[-1] * X.mean / X.sd),
        beta[-1] / X.sd
      )
    }else{
      beta2 <- beta
    }
    
  }else{
    beta2 <- beta
  }
  
  

  if(length(beta2) > 1){
    Vs <- sum(abs(beta2[-1]) >= threshold)
  }else{
    Vs <- 0
  }
  
  gamma <- 1  # gamma in [0, 1]
  

  P <- max(ncol(X) - 1, 0)
  
  EBIC <- -2 * logl +Vs * log(K1) +2 * gamma * lchoose(P, Vs)
  
  
  return(list(
    EBIC = EBIC,
    logl = logl,
    fhat = fhat,
    Fhat = Fhat
  ))
}