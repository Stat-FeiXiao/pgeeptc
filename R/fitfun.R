#' Title  Selection of the tuning parameter for the penalized generalized estimating equations
#'
#' @description We select the tuning parameters based on the BIC and fit the PTC model using SCAD-based penalized estimating equations.
#'
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param Var if it is TRUE, the program returns Std.Error.
#' @param X a matrix of covariates that may have effect on the failure times.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}.
#' @param stad if it is TRUE, the covariates are standardized.
#' @param beta_int the initial value of the covariate coefficients.
#' @param nopindex a vector indicating which covariates are not penalized.
#' @param nlambda the number of \code{lambda} values.
#' @param lambda.min.ratio If the user does not provide \code{lambda}, the project will calculate the maximum of \code{lambda} and multiply it by \code{lambda.min.ratio} to determine the minimum of \code{lambda}. The default is \code{lambda.min.ratio=1e-4}.
#' @param eps tolerance for convergence. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param lambda a user supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence based on \code{nlambda}. Supplying a value of \code{lambda} overrides this.
#' @param maxiter specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{maxiter} iterations and the estimates will be based on the last iteration.
#' @param tol estimated coefficients below \code{tol} are considered as zero.
#' @param QIC if it is TRUE, the program returns the quasilikelihood information criterion score of the fitted model.
#'
#' @export fitfun

fitfun <- function(Time, Status, id, Var, X, corstr , stad, beta_int, nopindex=NULL, nlambda, lambda.min.ratio, eps, lambda=NULL, maxiter, tol, QIC){

  N <- length(Time)
  K <- length(unique(id))
  p<-ncol(X)
  stadX <- X

  if(stad){
    Ibeta1 <- rep(0,length(beta_int))
    cumIbe <- 0
    for(lln in 2:length(beta_int)){
      cumIbe <- mean(X[,lln]) * beta_int[lln]
      Ibeta1[lln] =  sd(X[,lln ]) * beta_int[lln]
    }
    Ibeta1[1] <- beta_int[1] + cumIbe
    beta_int <- Ibeta1

    for (i in 2:ncol(X)) {
      stadX[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])
    }
  }

  if(is.null(lambda)){
    KM.fit <- survfit(Surv(Time, Status)~1)
    cure.rate <- min(KM.fit$surv)
    bet.init <- c(ifelse(cure.rate==0,0,log(-log(cure.rate))),rep(0,p-1))
    Bthat <- optim(par=bet.init,fn=Loss, control=list(maxit=maxiter,fnscale=-1,reltol=eps), Time=Time, Status= Status, X=X)$par
    IS <- lambda_IS(beta=Bthat,Time=Time,Status=Status,X=X)

    if(is.null(nopindex)){
       lambdaMax <-  max(  abs((diag(IS$I)*Bthat+IS$S)[-1]*Bthat[-1])   )
    }else{
      lambdaMax <-  max(  abs((diag(IS$I)*Bthat+IS$S)[-1]*Bthat[-1])[-nopindex]     )
    }
    lamgrid <- exp(seq(log(lambdaMax*1.5), log(lambdaMax * lambda.min.ratio), length.out = nlambda))

  }else{
    lamgrid <- lambda
    nlambda <- length(lamgrid)
  }
  L_grid <- min(lamgrid)
  R_grid <- max(lamgrid)

  BIC<-rep(0,nlambda)
  res <- list()
  for(s in 1:nlambda){

    lam<- lamgrid[s]
    fit_PGEE <- pgee(Time,Status, id,stadX,corstr ,beta_int,nopindex, lam,eps , maxiter)
    beta1 <- fit_PGEE$beta
    lam1 <- fit_PGEE$lam
    bic <- BIC(X,Status,Time, id,beta1,stad,tol)
    BIC[s] <- fit_PGEE $ BIC <- as.double( bic$BIC )
    fit_PGEE $ logl <-  bic$logl
    fit_PGEE $ fhat <-  bic$fhat
    fit_PGEE $ Fhat <-  bic$Fhat
    res[[s]] <- fit_PGEE
  }

  fit <- res[[which.min(BIC)]]
  fit$lam.minBIC <- fit$lam
  beta1 <- fit$beta

  if(stad){
    beta2 <- c(beta1[1] - sum((beta1[-1] * apply(X[, -1, drop = FALSE], 2, mean)/apply(X[, -1, drop = FALSE], 2, sd))),
              beta1[-1]/apply(X[,  -1, drop = FALSE], 2, sd))
  }else{ beta2 <- beta1 }
  fit$beta <- beta.minBIC <- beta2
  beta2 <- as.matrix(beta2)

  R.fi.hat=mycor(Status,id,X,beta2,Time,corstr)
  fit$  Qhat=R.fi.hat$Qhat
  fit$  pphi=R.fi.hat$pphi
  fit$  rho=R.fi.hat$rho

  for(ss in 1:length(beta.minBIC) ){
  if(abs(beta.minBIC[ss]) < tol) beta.minBIC[ss] <- 0
  }

  fit$beta.minBIC <- beta.minBIC

  if(Var){

    U.H.E.M.val=U_H_E_M(Time,Status,id,X,beta2,fit$Qhat,fit$pphi,fit$lam.minBIC,nopindex)
    H<-U.H.E.M.val$H
    E<-U.H.E.M.val$E
    M<-U.H.E.M.val$M

    fit$V.b<-  MASS::ginv(H+K*E) %*%M%*%  t( MASS::ginv(H+K*E) )

    M1 <- M[beta.minBIC != 0 ,beta.minBIC != 0]
    H1 <- H[beta.minBIC != 0 ,beta.minBIC != 0]
    E1 <- E[beta.minBIC != 0 ,beta.minBIC != 0]

    bet.var <- rep(NA,length(beta.minBIC)); bet.var[beta.minBIC!=0] <- diag(MASS::ginv(H1+K*E1)%*%M1%*%t(MASS::ginv(H1+K*E1))  )
    fit$V.b1<- bet.var

  }

  if(QIC){

    OM <- matrix(0,dim(X)[2],dim(X)[2])

    for(r in 1:dim(X)[2] ){
      for(s in 1:dim(X)[2]){
        OM[r,s] <- - sum( - fit $ Fhat * (exp(X %*% beta2)*X[,r]*X[,s] )   )
      }
    }
    TrM <- sum(diag( OM%*%(fit$V.b) ) )


    fit $QIC.value <- -2*fit$logl+ 2*TrM
  }

  if( fit$lam.minBIC == L_grid |  fit$lam.minBIC ==R_grid ) {
    fit$lam.grad.improper <- TRUE
  }else{
    fit$lam.grad.improper <- FALSE
  }
  if( fit$ lam.minBIC  == R_grid ) {
    fit$indicateR <- TRUE
  }else{
    fit$indicateR <- FALSE
  }
  return(fit)
}


#' @title Loss function for beta
#'
#' @description  Loss function for the promotion cure model.
#'
#' @param beta the estimation of covariate coefficients.
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates that that may have effect on the failure times.
#'
#' @export Loss
Loss <- function(beta,Time,Status,X){

  mu <- exp(X %*% beta)
  N <- length(Time)
  dRbet <- sapply(1:length(Time),function(i){
    sum( mu*(Time>=Time[i]) )
  })
  Rbet.min <- min(dRbet[Status==1])
  interval <- c(Rbet.min-sum(Status),Rbet.min-1)
  lambet <- uniroot(getlambda,interval,tol = .Machine$double.eps^0.75,
                    dRbet=dRbet,Status=Status)$root

  loss <- sum( log(  mu[Status==1]*N / (dRbet[Status==1]-lambet)  ) ) - lambet

  return(loss)
}


#' @title calculate the maximum lambda value
#'
#' @description calculate the maximum lambda value.
#'
#' @param beta the estimation of covariate coefficients.
#' @param Time right censored data which is the follow up time.
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates that that may have effect on the failure times.
#'
#' @export lambda_IS
lambda_IS <-  function(beta,Time,Status,X){

  Time1 <- Time[Status==1]
  tau <- max(Time[Status==1])
  Delta <- 1*(Time<=tau)
  X1 <- X[Status==1,,drop=F]
  Kn <- length(Time)
  mu <- as.vector(exp(X %*% beta))

  dRbet <- sapply(1:length(Time),function(i){
    sum( mu*(Time>=Time[i]) )
  })
  Rbet.min <- min(dRbet[Status==1])

  interval <- c(Rbet.min-sum(Status),Rbet.min-1)
  lambet <- uniroot(getlambda,interval,tol = .Machine$double.eps^0.75,
                    dRbet=dRbet,Status=Status)$root

  Rbetd <- dRbet[Status==1]

  Dd <- t(sapply(1:length(Time1),function(j){apply(X*(mu* (Delta*(Time>=Time1[j])+1-Delta) ),2,mean)}))
  ch2 <- sum(1/(Rbetd*(Rbetd-lambet)))
  ch1 <- apply(Dd/(Rbetd*(Rbetd-lambet)),2,sum)
  ch <- ch1/ch2
  hd <- t(t(Dd)-ch)/Rbetd
  S1 <- apply(X1/ Kn-hd,2,sum)
  S2 <- ch
  S <- (S1 - S2)
  I <- 0
  for(j in 1:length(Time1)){
    I <- I +
      t( t(t(X)-Kn*hd[j,])*(mu*(Delta*(Time>=Time1[j])+1-Delta)) )%*%t(t(X/Kn)-hd[j,])  / (Rbetd[j]-lambet)
  }

  return(list(I=I,S=S))
}

