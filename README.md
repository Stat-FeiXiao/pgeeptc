# pgeeptc
This is an R package for variable selection in marginal promotion time cure model under clustered failure time via SCAD-based penalized generalized estimating equations (PGEE).
We consider three working correlation matrices in the PGEE: *independent*,  *exchangeable* and  *AR-1*.
We establish a robust sandwich variance estimation formula. In addition, a bootstrap procedure for variance estimation of the covariate coefficients and correlation parameters is also provided.

## Package description and included main functions

Installation of this package can be done locally after downloading the package manually from this github website. We will also upload this package to the Comprehensive R Archive Network (CRAN) so that it can be downloaded as a standard R package. Currently, it can be loaded using R command
```R
devtools::install_github("Stat-FeiXiao/pgeeptc")
library(qifptc)
```

The main function included in our R package is *pgeeptc()* and there is also a function *print.pgeeptc()* for printing fitted results with a better presentation. To sum up, they can be called via:
```R
pgeeptc(formula, id, data, corstr = "independence", stad=TRUE, Ibeta=NULL, Var=FALSE, lambda=NULL, nopindex=NULL, boots=FALSE, nboot=100, QIC =FALSE, 
        nlambda=100, lambda.min.ratio=1e-4, eps = 1e-5, maxiter = 100, tol = 1e-3){
```
- **print.pgeeptc**: print outputted results from the previous function *pgeeptc()* with syntax
```R
print.pgeeptc(fit)
```

## An example using a periodontal disease data is shown below:

```R
## library
library(survival)
library(survminer)

#### Data preparation
```R
data(teeth)
n <- 9
id1 <- as.numeric(names(table(teeth$id)))[as.numeric(table(teeth$id))==n]
K <- sum(as.numeric(table(teeth$id))==n)
Data <- teeth[teeth$id==id1[1],]
for(i in 2:K){
  Data <- rbind(Data,teeth[teeth$id==id1[i],]) 
}
Data $ id <- rep(1:K,each=n)
Data$Mg <- Data$x10 # 1 for Mucogingival defect
Data$Endo <- Data$x16 # 1 for endo Therapy
Data$Decay <- Data$x21 # 1 for decayed tooth
Data$Gender <- Data$x49 # 1 for female
```

## plot a figure to show the existence of a cure fraction
```R
ggsurvplot(survival::survfit(survival::Surv(time, event) ~ Gender, data = Data), 
           ylim = c(0.6,1),
           ylab = "Survival Probability", xlab = "Survival Time (in Years)", 
           censor.shape="+",
           legend.title = "Gender",
           legend.labs = c("Male","Female")
)
```
![Teeth_KM_Gender](https://github.com/user-attachments/assets/e5fd1984-d3c6-4b55-a40b-43d631ec7b29)


#### Fit the marginal semi-parametric promotion time cure model using GEE method
- exchangeable correlation
```R
teeth.gee.ex <- qifptc(
        formula = Surv(time, event) ~ Gender + Mg + Endo + Decay, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "GEE", corstr="exchangeable", data = Data
)
print.qifptc(teeth.gee.ex)
```
- AR(1) correlation
```R
teeth.gee.ar1 <- qifptc(
        formula = Surv(time, event) ~ Gender + Mg + Endo + Decay, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "GEE", corstr="AR1", data = Data
)
print.qifptc(teeth.gee.ar1)
```
- independence correlation
```R
teeth.gee.ind <- qifptc(
        formula = Surv(time, event) ~ Gender + Mg + Endo + Decay, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "GEE", corstr="independence", data = Data
)
print.qifptc(teeth.gee.ind)
```
#### Fit the marginal semi-parametric promotion time cure model using QIF method
- exchangeable correlation
```R
teeth.qif.ex <- qifptc(
        formula = Surv(time, event) ~ Gender + Mg + Endo + Decay, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "QIF", corstr="exchangeable", data = Data
)
print.qifptc(teeth.qif.ex)
```
- AR(1) correlation
```R
teeth.qif.ar1 <- qifptc(
        formula = Surv(time, event) ~ Gender + Mg + Endo + Decay, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "QIF", corstr="AR1", data = Data
)
print.qifptc(teeth.qif.ar1)
```
- independence correlation
```R
teeth.qif.ind <- qifptc(
        formula = Surv(time, event) ~ Gender + Mg + Endo + Decay, 
        id = Data$id, Var = TRUE, stad=TRUE, method = "QIF", corstr="independence", data = Data
)
print.qifptc(teeth.qif.ind)
```
