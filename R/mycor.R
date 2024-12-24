#' Title Estimation of the working correlation matrix and scale parameter
#'
#' @description Estimate the working correlation matrix and scale parameter, where the structure of working correlation matrix can be specified as \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta_new the estimation of covariates coefficients.
#' @param Time right censored data which is the follow up time.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#'
#' @export mycor

mycor <- function(Status, id, X, beta_new, Time, corstr) {
  K <- length(unique(id))
  nt <- n <- as.vector(table(id))
  eta=X%*%beta_new
  mu=exp(eta)
  sd=sqrt(mu)

  t2 <- Time
  N <- length(Time)
  XI <- X
  t11 <- sort(Time[Status==1])
  kk <- length(table(Time[Status == 1]))

  mu <- exp(X %*% beta_new )

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

  #
  Fhat[Fhat==0] <- 1e-8
  y<-Status/Fhat
  res<-(as.vector(y)-mu)/sd
  maxclsz<-max(nt)

  pphi<-sum(res^2)/(sum(nt)-dim(X)[2])

  rres <- 0
  resm <- matrix(0, ncol = K, nrow = max(nt))
  for (i in 1:K) {
    resm[1:nt[i], i] <- res[id == i]
  }
  res <- resm
  res <- t(res)

  if (corstr=="independence"){
    rho<-0
  }else if(corstr=="exchangeable"){

    for (i in 1:K) {
      if (nt[i] == 1)
        rres <- rres + res[i, 1] else {
          for (j in 1:(nt[i] - 1)) rres <- rres + res[i, j] * sum(res[i, (j + 1):nt[i]])
        }
    }
    rho <- (pphi^(-1)) * rres/(sum(nt * (nt - 1))/2 - dim(X)[2])
  }else if(corstr=="AR1"){
    for (i in 1:K) {
      if (nt[i] == 1)
        rres <- rres + res[i, 1] else {
          for (j in 1:(nt[i] - 1)) rres <- rres + res[i, j] * res[i, j+1]
        }
    }
    rho <- (pphi^(-1)) * rres/(sum((nt - 1)) - dim(X)[2]  )

  }

  Qhat<-array(0,c(maxclsz,maxclsz,K))

  for(i in 1:K){
    cor1<-matrix(0,nt[i],nt[i])
    if (corstr=="independence"){
      cor1<-diag(nt[i]) }else if(corstr=="exchangeable"){
        for (t1 in 1:nt[i]) {
          for (t2 in 1:nt[i]) {
            if (t1!=t2)
            {cor1[t1,t2]<-rho} else
            {cor1[t1,t2]<-1}
          }}}else if(corstr=="AR1"){
            exponent <- abs(matrix(1:nt[i] - 1, nrow = nt[i], ncol = nt[i], byrow = TRUE) - (1:nt[i] - 1))
            cor1 <- rho^exponent
          }
    Qhat[1:nt[i],1:nt[i],i]<-cor1
  }

  return(list("Qhat"=Qhat,"pphi"=pphi,"rho"=rho))

}
