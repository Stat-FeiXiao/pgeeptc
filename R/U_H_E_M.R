#' Title The element used in iterative algorithm
#' 
#' @description Generate specific matrices for the iterative algorithm process.
#' 
#' @param Status the censoring indicator, 1 = event of interest happens, and 0 = censoring.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param X a matrix of covariates that that may have effect on the failure times.
#' @param beta_new the given values of covariates coefficients.
#' @param Qhat the estimation of working correlation matrix.
#' @param pphi the estimation of scale parameter.
#' @param Fhat the estimation of baseline cumulative hazard function.
#' @param lam the specific tuning parameter.
#' @param eps a fixed small number used in the iterative algorithm. The default is \code{eps = 1e-6}.
#'
#' @export U_H_E_M

U_H_E_M <-
  function(Status, id, X, beta_new, Qhat, pphi, Fhat, lam, eps=1e-6) {
    K <- length(unique(id))
    n<-nt <- as.vector(table(id))
    aindex=cumsum(nt)
    index=c(0,aindex[-length(aindex)])

    y<-Status/Lambda
    mu=exp(X%*%beta_new)
    nvars <- ncol(X)
    lam<-rep(lam,nvars)
    
    E<- diag( q_scad(abs(as.vector(beta_new)),lam)/(abs(as.vector(beta_new))+eps) )
    E[,1]<-0

    sum201<-matrix(0,(nvars ),1)         #gradient:S
    sum301<-matrix(0,(nvars ),(nvars ))  #naive variance:H
    sum401<-matrix(0,(nvars ),(nvars ))  #a component for robust variance:M
    sum501<-0
    for (i in 1:K) {
      ym<-matrix(0,nt[i],1)
      bigD<-matrix(0,nt[i],(nvars))
      bigA<-matrix(0,nt[i],nt[i])
      bigX<-matrix(0,nt[i],(nvars))
      for (j in 1:nt[i]) {
    
        ym[j]<- y[j+index[i]]-mu[j+index[i]]
        bigA[j,j]<-mu[j+index[i]]
        for (k in 1:(nvars )) {
          bigX[j,k]<-X[j+index[i],k]
        
        } 
      } 

   
      bigV<-Qhat[1:nt[i],1:nt[i],i]

      Wla<-pphi^(-1)*diag(Fhat[id==i])

      sum200<-t(bigX)%*%sqrt(bigA)%*%MASS::ginv(bigV)%*%MASS::ginv(sqrt(bigA))%*%Wla%*%ym    
      sum201<-sum201+sum200


      sum300<-t(bigX)%*%sqrt(bigA)%*%MASS::ginv(bigV)%*%MASS::ginv(sqrt(bigA))%*%Wla%*%(bigA)%*%bigX   
      sum301<-sum301+sum300

 
      SSA=MASS::ginv(sqrt(bigA))
      SSAym=(SSA%*%Wla%*%ym)

      sum400<-t(bigX)%*%sqrt(bigA)%*%MASS::ginv(bigV)%*%(SSAym)
      sum400  <-       sum400 %*% t(sum400)
      sum401<-sum401+sum400

      ##Dev
      sum500<-t(ym)%*%MASS::ginv(sqrt(bigA))%*%MASS::ginv(bigV)%*%MASS::ginv(sqrt(bigA))%*%Wla%*%ym
      sum501<-sum501+sum500

     
    } 

    U<-sum201
    H<-sum301
    E<-E
    M<-sum401
    Dev <-sum501
    
    return(list("U"=U,"H"=H,"E"=E,"M"=M,"Dev"=Dev))
  }


#' Title The nonconvex SCAD penalty
#'
#' @description Generate the SCAD penalty vector
#'
#' @param beta the given values of covariate coefficients.
#' @param lambda the specified tuning parameter
#' @param a  a  parameter in the SCAD penalty.  The default is \code{a = 3.7}.
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

