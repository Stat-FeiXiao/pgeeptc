#' Title Estimate the baseline cumulative distribution function in the PTC model using Bernstein polynomial
#'
#' @description The estimation of baseline survival function based on the Bernstein polynomial.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 0 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param beta estimation of the covariate coefficients based on the GEE or QIF approach.
#' @param N the degree of Bernstein polynomial used to approximate the baseline distribution function in model.
#' @param tau the cure threshold. Individuals with survival time greater than \code{tau} are considered cured.
#'
#' @export baseF

baseF<- function(Time, Status, X, beta, N, tau) {
  Kn <- length(Time)
  t2 <- Time
  t11 <- sort(Time)
  c11 <- Status[order(Time)]
  x111 <- as.matrix(X[order(Time), ])
  tt1 <- unique(t11[c11 == 1])
  
  
  
  b.kN <- function(t,k,N){  
    
    y <- choose(N,k) * (t/tau)^k * (1 - t/tau)^(N-k)
    return(y)
  }
  
  D_b.kN = function(t,k,N){
    
    result = (t/tau)^(k - 1) * (k * (1/tau)) * (1 - t/tau)^(N - k) - (t/tau)^k * ((1 - t/tau)^((N - k) - 1) * ((N - k) * (1/tau)))
    y = choose(N,k) * result
    return(y)
  }
  
  f.t.star = function(t.star, gamma){
    
    A.derive <- matrix(0, nrow = length(t.star), ncol = N+1)
    for(k in 0:N){
      A.derive[,k+1] <- D_b.kN(t.star,k,N) 
    }
    
    exp_ga <- exp(gamma)
    PPHI <- cumsum(exp_ga)/sum(exp_ga)
    
    f.t = A.derive %*% PPHI 
    f.t = as.vector(f.t)
    return(f.t)
  }
  
  
  F.t.star = function(t.star, gamma){
    
    A <- matrix(0, nrow = length(t.star), ncol = N+1)
    for(k in 0:N){
      A[,k+1] <- b.kN(t.star,k,N) 
    }
    exp_ga <- exp(gamma)
    PPHI <- cumsum(exp_ga)/sum(exp_ga)
    
    F.t = A %*% PPHI 
    F.t = as.vector(F.t)
    return(F.t)
  }
  
  A <- matrix(0, nrow = length(Time), ncol = N+1)
  for(k in 0:N){
    A[,k+1] <- b.kN(Time,k,N) 
  }
  
  A.derive <- matrix(0, nrow = length(Time), ncol = N+1)
  for(k in 0:N){
    A.derive[,k+1] <- D_b.kN(Time,k,N) 
  }
  
  
  Loglike<-function(gam){  
    
    gamma0 <- gam
    exp_ga <- exp(gamma0)
    PPHI <- cumsum(exp_ga)/sum(exp_ga)
    
    F.Yi = A %*% PPHI 
    F.Yi = as.vector(F.Yi)
    
    f.Yi = A.derive %*% PPHI  
    Like = Status * log(f.Yi) - (exp(X %*% beta) * F.Yi) 
    Like= sum(Like)
    return(Like)
  }
  
  
  Eqtion = function(gam){
    gamma0 <- gam
    exp_ga <- exp(gamma0) 
    PPHI <- cumsum(exp_ga)/sum(exp_ga)
    DPPHI <- eql <- rep(0,(N+1))
    f.Yi = A.derive %*% PPHI  
    
    for(s in 1:(N+1)){
      for(j in 1:(N+1)){
        if(j >= s) {
          DPPHI[j] <- exp_ga[s]* (sum(exp_ga)- sum(exp_ga[1:j]))/((sum(exp_ga))^2)}else{
            DPPHI[j] <-  exp_ga[s]*  (- sum(exp_ga[1:j]) )/((sum(exp_ga))^2)
          }
      }
      DF = A %*% DPPHI
      Df = A.derive %*% DPPHI
      
      
      eql [s]<-sum( Status * (1/f.Yi) * Df - exp(X %*% beta) * DF)
      
      
    }
    
    
    
    return(eql)
  }
  
  
  
  re_esti<-maxLik::maxLik(Loglike, start = rep(0,N+1), method="BFGS", grad = Eqtion )$estimate
  
  re_esti<-as.double(re_esti)
  Fhat <-  F.t.star(Time, re_esti)
  fhat <- f.t.star(Time, re_esti)
  
  list(Fhat = Fhat,fhat=fhat,gam=re_esti)
}