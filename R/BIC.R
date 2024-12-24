#' Title The Bayesian information criterion for selection of the tuning parameter
#'
#' @description Calculate the BIC score for selection of the tuning parameter.
#'
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param Time right censored data which is the follow up time.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param beta the estimation of covariates coefficients.
#' @param stad if it is TRUE, all the covariates in the \code{formula} are standardized.
#' @param threshold estimated coefficients below \code{threshold} are considered as zero.
#'
#' @export BIC

BIC<-function(X, Status, Time, id, beta, stad, threshold){
   
  K1 <- length(unique(id) )
  n <- as.vector(table(id))
  N <- length(id)
  kk <- length( table(Time[Status == 1]) )
  Repf <- as.numeric( table(Time[Status == 1]) )
  beta <- as.double(beta)
  stadX <- X
  if(stad){
    for (i in 2:ncol(X)) {
      stadX[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
    }
  }

  XI <- stadX[Status==1,]
 
  mu <- exp(stadX %*% beta )
  
  dRbet <- sapply(1:sum(n),function(i){
    sum( mu*(Time>=Time[i]) )
  })
  Rbet.min <- min(dRbet[Status==1])
  interval <- c(Rbet.min-sum(Status),Rbet.min-1)
  lambet <- uniroot(getlambda,interval,tol = .Machine$double.eps^0.75,
                    dRbet=dRbet,Status=Status)$root
  
  Fhat <- sapply(Time,function(tmj){
    sum( (Time<=tmj)[Status==1]/(dRbet[Status==1]-lambet) )
  })
  
  DeltaF <- sort(unique(Fhat))
  if(DeltaF[1] == 0){ DeltaF <- DeltaF[-1] }
  fhat <- c(DeltaF[1], DeltaF[2:kk] - DeltaF[1:(kk - 1)])

  fhat <- rep(fhat,times=Repf)
  #
  fun1<- log(exp(XI%*%beta)*fhat)
  fun2<- -exp(stadX%*%beta)*Fhat
  logl<- sum(fun1)+sum(fun2)

  if(stad){
    beta2<- c(beta[1] - sum((beta[-1] * apply(X[, -1, drop = FALSE], 2, mean)/apply(X[, -1, drop = FALSE], 2, sd))),
              beta[-1]/apply(X[,  -1, drop = FALSE], 2, sd))
  }else{ beta2 <- beta }
  
  BIC<- ( - logl    +  0.5*log(K1)*sum(abs(beta2[-1])>=threshold) )/N   
  
 
  return(list(BIC = BIC,logl=logl,fhat=fhat,Fhat=Fhat))
}
