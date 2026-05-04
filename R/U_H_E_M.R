#' Title The element used in iterative algorithm
#' 
#' @description Generate specific matrices for the iterative algorithm process.
#'
#' @param Time right censored data which is the follow up time. 
#' @param Status the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta_new the estimation of covariates coefficients.
#' @param Qhat the estimation of working correlation matrix.
#' @param pphi the estimation of scale parameter.
#' @param lam the tunning parameter in SCAD penalty.
#' @param nopindex a vector indicating which covariates are not penalized.
#'
#' @export U_H_E_M

U_H_E_M <- function(Time, Status, id, X, beta_new, Qhat, pphi, lam, nopindex) {
  
    K <- length(unique(id))
    n<-nt <- as.vector(table(id))
    aindex=cumsum(nt)
    index=c(0,aindex[-length(aindex)])
    mu=exp(X%*%beta_new)
    nvars <- ncol(X)
    lam<-rep(lam,nvars)
    
    E<- diag( q_scad(abs(as.vector(beta_new)),lam)/(abs(as.vector(beta_new))+1e-6) )
    E[,1]<-0

    if( ! is.null(nopindex) ) E[,nopindex+1]<-0
    
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

    sum201<-matrix(0,(nvars ),1)       
    sum301<-matrix(0,(nvars ),(nvars ))
    sum401<-matrix(0,(nvars ),(nvars ))  
    sum501<-0
    for (i in 1:K) {
      ym<-matrix(0,nt[i],1)
      bigD<-matrix(0,nt[i],(nvars))
      bigA<-matrix(0,nt[i],nt[i])
      bigX<-matrix(0,nt[i],(nvars))
      for (j in 1:nt[i]) {
        
        ym[j]<- Status[j+index[i]]-(Fhat*mu)[j+index[i]]
        bigA[j,j]<-mu[j+index[i]]
        for (k in 1:(nvars )) {
          bigX[j,k]<-X[j+index[i],k]
          
        } 
      } 
      
      bigV<-Qhat[1:nt[i],1:nt[i],i]
      if(nt[i]==1){
        Wla<- Fhat[id==i]
      }else Wla<-diag(Fhat[id==i])

      sum200<-t(bigX)%*%sqrt(bigA)%*%MASS::ginv(pphi*bigV)%*%MASS::ginv(sqrt(bigA))%*%ym    
      sum201<-sum201+sum200
      
      
      sum300<-t(bigX)%*%sqrt(bigA)%*%MASS::ginv(pphi*bigV)%*%MASS::ginv(sqrt(bigA))%*%Wla%*%(bigA)%*%bigX   
      sum301<-sum301+sum300
      
      
      SSA=MASS::ginv(sqrt(bigA))
      SSAym=(SSA%*%ym)
      
      sum400<-t(bigX)%*%sqrt(bigA)%*%MASS::ginv(pphi*bigV)%*%(SSAym)
      sum400  <-       sum400 %*% t(sum400)
      sum401<-sum401+sum400
      
    } 
    
    U<-sum201
    H<-sum301
    E<-E
    M<-sum401

    
    return(list("U"=U,"H"=H,"E"=E,"M"=M))
  }


#' Title The nonconvex SCAD penalty
#'
#' @description Generate the SCAD penalty vector
#'
#' @param beta the estimation of covariates coefficients.
#' @param lambda the specified tuning parameter
#' @param a a parameter in the SCAD penalty. The default is \code{a=3.7}.
#' 
#' @export q_scad

q_scad <-
  function(beta,lambda,a=3.7)
  {
    
    p<-length(beta)
    
    beta<-abs(beta)
    
    b1<-rep(0,p)
    
    b1[beta>lambda]<-1
    
    b2<-rep(0,p)
    
    b2[beta<(lambda*a)]<-1
    
    q<-lambda*(1-b1)+((lambda*a)-beta)*b2/(a-1)*b1
    return(q)
  }

#' Title An auxiliary function for solving lambda given beta
#'
#' @description An auxiliary function for solving lambda given beta
#'
#' @param lam unknown parameter corresponding to the largrane multiplier.
#' @param dRbet a matrix needed in this function.
#' @param Status the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' 
#' @export getlambda
 
getlambda <- function(lam,dRbet,Status){
  
  val <- sum( 1/( dRbet [Status ==1]- lam )) - 1
  
  return(val)
  
}
