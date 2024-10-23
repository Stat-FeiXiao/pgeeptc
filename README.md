# pgeeptc
This is an R package for variable selection in marginal promotion time cure model under clustered failure time via SCAD-based penalized generalized estimating equations (PGEE).

We consider three working correlation matrices in the PGEE: *independent*,  *exchangeable* and  *AR-1*.

We establish a robust sandwich variance estimation formula. In addition, a bootstrap procedure for variance estimation of the covariate coefficients and correlation parameters is also provided.

## Package description and included main functions

Installation of this package can be done locally after downloading the package manually from this github website. We will also upload this package to the Comprehensive R Archive Network (CRAN) so that it can be downloaded as a standard R package. Currently, it can be loaded using R command
```R
devtools::install_github("Stat-FeiXiao/pgeeptc")
library(pgeeptc)
```

The main function included in our R package is *pgeeptc()* and there is also a function *print.pgeeptc()* for printing fitted results with a better presentation. To sum up, they can be called via:
- **pgeeptc**: fit the models in various ways with synopsis
```R
pgeeptc(formula, id, data, corstr = "independence", stad=TRUE, Ibeta=NULL, Var=FALSE, lambda=NULL, nopindex=NULL, boots=FALSE, nboot=100, QIC =FALSE, 
        nlambda=100, lambda.min.ratio=1e-4, eps = 1e-5, maxiter = 100, tol = 1e-3){
```
- **print.pgeeptc**: print outputted results from the previous function *pgeeptc()* with syntax
```R
print.pgeeptc(fit)
```

## An example using a tonsil cancer data is shown below:

```R
#### Data preparation
data(tonsil)
Data <- tonsil
Data$x1 <- as.numeric( tonsil$Cond == 1)
Data$x2 <- as.numeric( tonsil$T == 4 )
Data$x3 <- as.numeric( tonsil$N == 3 )
Data$x4 <- as.numeric( tonsil$Grade == 1 )
Data$x5 <- Data$x1 * Data$x2
Data$x6 <- Data$x1 * Data$x3
Data$x7 <- Data$x1 * Data$x4
Data$x8 <- Data$x2 * Data$x3
Data$x9 <- Data$x2 * Data$x4
Data$x10 <- Data$x3 * Data$x4
Data$id <- tonsil$Inst
Data$time <- tonsil$Time
Data$event <- tonsil$Status
```

### Variable selection for the marginal semi-parametric promotion time cure model using PGEE method
- GEE without penalty
```R
set.seed(123)

GEEex_fit <- pgeeptc(formula=Surv(time,event)~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,id=Data$id, 
              data=Data, corstr ='exchangeable',stad=T, Ibeta=NULL, Var=T, lambda=0, 
              nopindex=NULL, boots=T, nboot=100, QIC =T, nlambda=100, lambda.min.ratio=1e-8,eps = 1e-6, maxiter = 100, tol = 1e-3)

print.pgeeptc(GEEex_fit)
```
- PGEE with exchangeable working correlation matrix
```R
PGEEex_fit <- pgeeptc(formula=Surv(time,event)~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,id=Data$id, 
                      data=Data, corstr ='exchangeable',stad=T, Ibeta=NULL, Var=T, lambda=NULL, 
                      nopindex=NULL, boots=T, nboot=100, QIC =T, nlambda=100, lambda.min.ratio=1e-8,eps = 1e-6, maxiter = 100, tol = 1e-3)

print.pgeeptc(PGEEex_fit)
```
- PGEE with independent working correlation matrix
```R
PGEEind_fit <- pgeeptc(formula=Surv(time,event)~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,id=Data$id,
                       data=Data, corstr ='independence',stad=T, Ibeta=NULL, Var=T, lambda=NULL, 
                       nopindex=NULL, boots=FALSE, nboot=100, QIC =T, nlambda=100, lambda.min.ratio=1e-8,eps = 1e-6, maxiter = 100, tol = 1e-3)

print.pgeeptc(PGEEind_fit)
```
