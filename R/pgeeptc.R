#============== The main function for variable selection in the PTC model based on PGEE ==============#

#' Title  Variable selection in the marginal semiparametric promotion time cure model based on PGEE
#'
#' @description Variable selection in the marginal semiparametric promotion time cure  (PTC) model based on penalized general estimation equations (PGEE) method.
#' We consider three common correlation structures in this funciton.
#' The variances can be consistently estimated by a sandwich variance estimator. We also present a bootstrap procedure for variance estimation.
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring. The expression to the right of the "~" specifies the effect of covariates on the failure time.
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations.
#' @param data  a data frame in which to interpret the variables named in the \code{formula}.
#' @param corstr a character string specifying the correlation structure. The following are permitted: \code{independence}, \code{exchangeable} and \code{AR1}. The default is \code{corstr = "independence"}.
#' @param stad if it is TRUE, all the covariates in the \code{formula} are standardized. By default, \code{stdz = TRUE}.
#' @param Ibeta the initial value of the covariate coefficients.
#' @param Var if it is TRUE, the program returns Std.Error. By default, \code{Var = TRUE}.
#' @param lambda a user supplied \code{lambda} sequence. Typical usage is to
#' have the program compute its own \code{lambda} sequence based on \code{nlambda}. Supplying a value of \code{lambda} overrides this.
#' @param nopindex a vector indicating which covariates are not penalized.
#' @param boots if it is TRUE, the program returns Std.Error by the bootstrap method. By default, \code{boots = FALSE}.
#' @param nboot the number of bootstrap samples. The default is \code{nboot = 100}.
#' @param QIC if it is TRUE, the program returns the quasilikelihood information criterion score of the fitted model. By default, \code{QIC = FALSE}.
#' @param nlambda the number of \code{lambda} values. The default is \code{nlambda = 100}.
#' @param lambda.min.ratio If the user does not provide \code{lambda}, the project will calculate the maximum of \code{lambda} and multiply it by \code{lambda.min.ratio} to determine the minimum of \code{lambda}. The default is \code{lambda.min.ratio=1e-4}.
#' @param eps tolerance for convergence. Iteration stops once the relative change in deviance is less than \code{eps}. The default is \code{eps = 1e-5}.
#' @param maxiter specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{maxiter} iterations and the estimates will be based on the last iteration. The default \code{maxiter = 100}.
#' @param tol estimated coefficients below \code{tol} are considered as zero. The default is \code{tol = 1e-3}.
#'
#' @return An object of class \code{pgeeptc} is returned. It can be examined by \code{print.pgeeptc()}.
#'
#' @export pgeeptc

pgeeptc<-function(formula, id, data, corstr = "independence", stad=TRUE, Ibeta=NULL, Var=FALSE, lambda=TRUE, nopindex=NULL, boots=FALSE, nboot=100, QIC =FALSE, nlambda=100, lambda.min.ratio=1e-4, eps = 1e-5, maxiter = 100, tol = 1e-3){
  call <- match.call()

  uid <- sort(unique(id))
  newid <- rep(0, length(id))
  for (i in 1:length(id)) {
    j <- 1
    repeat {
      if (id[i] != uid[j])
        j <- j + 1 else {
          newid[i] <- j
          break
        }
    }
  }
  data$id <- newid
  data1 <- data[data$id == 1, ]
  for (i in 2:length(uid)) {
    data1 <- rbind(data1, data[data$id == i, ])
  }
  data <- data1
  id <- data$id
  nobs <- length(id)
  K <- length(unique(id))
  nt <- as.vector(table(id))
  maxclsz <-max(nt)
  nobs<-sum(nt)
  mf <- stats::model.frame(formula, data)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  X <- as.matrix(X)
  colnames(X) <- colnames(stats::model.matrix(attr(mf, "terms"), mf))
  beta_name <- colnames(X)
  beta_length <- ncol(X)
  Y <- stats::model.extract(mf, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  Time <- Y[, 1]
  Status <- Y[, 2]

  if(length(id) != length(Time))  stop("Id and Time do not have the same length!")

  if(!(is.double(X)))  X <- as.double(X)
  if(!(is.double(Time)))  Time <- as.double(Time)
  if(!(is.double(id))) id <- as.double(id)

  if(is.null(Ibeta)) {
    cvfit <- cv.glmnet(X[,-1],  Surv(Time,Status) , family = "cox", alpha = 1)
    coef<-as.numeric( coef(cvfit,s="lambda.min") )
    KMfit <- survfit(Surv(Time,Status)~1)
    cure.rate <- min(KMfit$surv)
    Ibeta <- c(ifelse(cure.rate==0,0,log(-log(cure.rate))),coef)
  }

  p<-ncol(X)
  if(beta_length != p) {stop("Dimension of beta != ncol(X)!")}

  if(!(is.double(K)))      K <- as.double(K)
  if(!(is.double(maxclsz)))  maxclsz <- as.double(maxclsz)
  if(!(is.double(nobs)))   nobs <- as.double(nobs)

  corstrs <- c("independence", "exchangeable", "AR1")

  corstrv <- as.integer(match(corstr, corstrs, -1))
  if(corstrv < 1) stop("unknown corstr!")

  if (!Var & QIC) stop("Calculating QIC needs variance estimation!")

  fitf <- fitfun(Time,Status, id,Var,X,corstr ,stad,Ibeta, nopindex,nlambda,lambda.min.ratio,eps,lambda,maxiter, tol,QIC)

  fitf$num_of_clusters <- K
  fitf$max_cluster_size <- maxclsz

  fitf$call <- call
  fitf$Var <- Var
  fitf$QIC <- QIC
  fitf$boots <- boots
  fitf$beta_name <- beta_name

  if(Var){

    fitf$beta_var <- fitf$V.b1
    fitf$beta_sd <-  sqrt(fitf$beta_var)
    fitf$beta_zvalue <- fitf$beta/fitf$beta_sd
    fitf$beta_pvalue <- (1 - stats::pnorm(abs(fitf$beta_zvalue))) * 2
  }

  if(boots){

    boots_sd<-boots( formula,data,  id, nboot, corstr, Ibeta , nopindex,  fitf$lam.minBIC, eps, maxiter)
    fitf $ boots_rho_sd <- boots_sd$sd_rho_boots
    fitf $ boots_beta_sd <- boots_sd$sd_beta_boots
    fitf$boots_beta_zvalue <- fitf$beta/fitf $ boots_beta_sd
    fitf$boots_beta_pvalue <- (1 - stats::pnorm(abs(fitf$boots_beta_zvalue))) * 2
    fitf$boots_rho_zvalue <- fitf$rho/fitf$boots_rho_sd
    fitf$boots_rho_pvalue <- (1 - pnorm(abs(fitf$boots_rho_zvalue))) * 2
  }


  class(fitf) <- "pgeeptc"
  return(fitf)
}

#' Title Print pgeeptc object
#'
#' @param x an object of \code{pgeeptc}.
#' @param ... further arguments to be added in the \code{print.pgeeptc} function.
#'
#' @export print.pgeeptc

print.pgeeptc <- function(x, ...) {
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nEstimation of Marginal Promotion Time Cure Model Based on PGEE:\n")
  if (x$Var) {
    bt <- array(c(x$beta.minBIC), c(length(x$beta.minBIC), 4))
    rownames(bt) <- c(x$beta_name)
    colnames(bt) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
    if(x$boots){
      bt[, 2] <- c(x$boots_beta_sd)
      bt[, 3] <- c(x$boots_beta_zvalue)
      bt[, 4] <- c(x$boots_beta_pvalue)
    }else{
      bt[, 2] <- c(x$beta_sd)
      bt[, 3] <- c(x$beta_zvalue)
      bt[, 4] <- c(x$beta_pvalue)
    }
  }
  else {
    bt <- array(c(x$beta), c(length(x$beta), 1))
    rownames(bt) <-  c(x$beta_name)
    colnames(bt) <- "Estimate"
  }
  print(bt)
  cat("\n")
  cat("\nEstimation of Correlation Parameters:\n")
  if (x$boots) {
    hatr <- array(x$rho, c(1, 4))
    rownames(hatr) <- "rho"
    colnames(hatr) <- c("Estimate", "Std.Error", "Z value", "Pr(>|Z|)")
    hatr[, 2] <- x$boots_rho_sd
    hatr[, 3] <- x$boots_rho_zvalue
    hatr[, 4] <- x$boots_rho_pvalue
  }
  else {
    hatr <- array(x$rho, c(1, 1))
    rownames(hatr) <- "rho"
    colnames(hatr) <- "Estimate"
  }
  print(hatr)
  cat("\n")
  if(x$QIC){
    cat("QIC:", x$QIC.value)
    cat("\n")
  }
  cat("Number of clusters:", x$num_of_clusters)
  cat("        Maximum cluster size:", x$max_cluster_size)
  cat("\n")
  invisible(x)
}


#==== teeth - A Periodontal Disease Data ====#
#' @title teeth - A Periodontal Disease Data
#'
#' @aliases teeth
#' @description Survival of teeth with various predictors.
#' @usage data(teeth)
#'
#' @format A data frame with 65,890 teeth on the following 56 variables.
#'  \describe{
#'    \item{x1}{numeric. \emph{mobil} Mobility score (on a scale 0--5).}
#'    \item{x2}{numeric. \emph{bleed} Bleeding on Probing (percentage).}
#'    \item{x3}{numeric. \emph{plaque} Plaque Score (percentage).}
#'    \item{x4}{numeric. \emph{pocket_mean} Periodontal Probing Depth (tooth-level mean).}
#'    \item{x5}{numeric. \emph{pocket_max} Periodontal Probing Depth (tooth-level mean).}
#'    \item{x6}{numeric. \emph{cal_mean} Clinical Attachment Level (tooth-level mean).}
#'    \item{x7}{numeric. \emph{cal_max} Clinical Attachment Level (tooth-level max).}
#'    \item{x8}{numeric. \emph{fgm_mean} Free Gingival Margin (tooth-level mean).}
#'    \item{x9}{numeric. \emph{fgm_max} Free Gingival Margin (tooth-level max).}
#'    \item{x10}{numeric. \emph{mg} Mucogingival Defect.}
#'    \item{x11}{numeric. \emph{filled} Filled Surfaces.}
#'    \item{x12}{numeric. \emph{decay_new} Decayed Surfaces -- new.}
#'    \item{x13}{numeric. \emph{decay_recur} Decayed Surfaces -- recurrent.}
#'    \item{x14}{numeric. \emph{dfs} Decayed and Filled Surfaces.}
#'    \item{x15}{numeric. \emph{crown} Crown.}
#'    \item{x16}{numeric. \emph{endo} Endodontic Therapy.}
#'    \item{x17}{numeric. \emph{implant} Tooth Implant.}
#'    \item{x18}{numeric. \emph{pontic} Bridge Pontic.}
#'    \item{x19}{numeric. \emph{missing_tooth} Missing Tooth.}
#'    \item{x20}{numeric. \emph{filled_tooth} Filled Tooth.}
#'    \item{x21}{numeric. \emph{decayed_tooth} Decayed Tooth.}
#'    \item{x22}{numeric. \emph{furc_max} Furcation Involvement for Molars.}
#'    \item{x23}{numeric. \emph{bleed_ave} Bleeding on Probing (mean percentage).}
#'    \item{x24}{numeric. \emph{plaque_ave} Plaque Index (mean percentage).}
#'    \item{x25}{numeric. \emph{pocket_mean_ave} Periodontal Probing Depth (mean of tooth mean).}
#'    \item{x26}{numeric. \emph{pocket_max_ave} Periodontal Probing Depth (mean of tooth max).}
#'    \item{x27}{numeric. \emph{cal_mean_ave} Clinical Attachment Level (mean of tooth mean).}
#'    \item{x28}{numeric. \emph{cal_max_ave} Clinical Attachment Level (mean of tooth max).}
#'    \item{x29}{numeric. \emph{fgm_mean_ave} Free Gingival Margin (mean of tooth max).}
#'    \item{x30}{numeric. \emph{fgm_max_ave} Free Gingival Margin (mean of tooth max).}
#'    \item{x31}{numeric. \emph{mg_ave} Mucogingival Defect (mean).}
#'    \item{x32}{numeric. \emph{filled_sum} Filled Surfaces (total).}
#'    \item{x33}{numeric. \emph{filled_ave} Filled Surfaces (mean).}
#'    \item{x34}{numeric. \emph{decay_new_sum} New Decayed Surfaces (total).}
#'    \item{x35}{numeric. \emph{decay_new_ave} New Decayed Surfaces (mean).}
#'    \item{x36}{numeric. \emph{decay_recur_sum} Recurrent Decayed Surfaces (total).}
#'    \item{x37}{numeric. \emph{decay_recur_ave} Recurrent Decayed Surfaces (mean).}
#'    \item{x38}{numeric. \emph{dfs_sum} Decayed and Filled Surfaces (total).}
#'    \item{x39}{numeric. \emph{dfs_ave} Decayed and Filled Surfaces (mean).}
#'    \item{x40}{numeric. \emph{filled_tooth_sum} Number of Filled Teeth.}
#'    \item{x41}{numeric. \emph{filled_tooth_ave} Percentage of Filled Teeth.}
#'    \item{x42}{numeric. \emph{decayed_tooth_sum} Number of Decayed Teeth.}
#'    \item{x43}{numeric. \emph{decayed_tooth_ave} Percentage of Decayed Teeth.}
#'    \item{x44}{numeric. \emph{missing_tooth_sum} Number of Missing Teeth.}
#'    \item{x45}{numeric. \emph{missing_tooth_ave} Percentage of Missing Teeth.}
#'    \item{x46}{numeric. \emph{total_tooth} Number of Teeth.}
#'    \item{x47}{numeric. \emph{dft} Number of Decayed and Filled Teeth.}
#'    \item{x48}{numeric. \emph{baseline_age} Patient Age at Baseline (years).}
#'    \item{x49}{numeric. \emph{gender} Gender.}
#'    \item{x50}{numeric. \emph{diabetes} Diabetes Mellitus.}
#'    \item{x51}{numeric. \emph{tobacco_ever} Tobacco Use.}
#'    \item{x52}{numeric. \emph{molar} Molar.}
#'    \item{id}{numeric. Patient ID.}
#'    \item{tooth}{numeric. Tooth ID.}
#'    \item{event}{numeric. Tooth Loss Status.}
#'    \item{time}{numeric. Follow Up Time.}
#'  }
#'
#' @details  The original data consist of 65228 people enrolled in a study to investigate the association between the time of
#' tooth loss from patients with periodontal disease and its relative covariates. The data is collected from patients treated at
#' Creighton University School of Dentistry from  August 2007 until March 2013.
#'
#' @keywords datasets
#'
"teeth"


#==== tonsil - A Tonsil Cancer data ====#
#' @title tonsil - Multi-Center Clinical Trial of Tonsil Carcinoma
#'
#' @aliases tonsil
#' @description A tonsil cancer clinical trial study conducted by the Radiation Therapy Oncology Group in the United States. The survival time is defined as the time (in days) from diagnosis to death. In this study, patients in one institution were randomly assigned to one of two treatment groups: radiation therapy alone or radiation therapy together with a chemotherapeutic agent. A part of the data from the study is available in Kalbfleisch and Prentice (2002).
#' @usage data(tonsil)
#'
#' @format A part of the data from the study is available in Kalbfleisch and Prentice (2002), which includes times (in days) from diagnosis to death of 195 patients with squamous cell carcinoma of three sites in the oropharynx between 1968 and 1972 in six participating institutions. Other variables include
#'  \describe{
#'      \item{\code{Inst}}{institution code, from 1 to 6, represents six participating institutions}
#'      \item{\code{Sex}}{1 = male, 2 = female.}
#'      \item{\code{Trt}}{treatment: 1 = standard, 2 = test.}
#'      \item{\code{Grade}}{1 = well differentiated, 2 = moderately differentiated, 3 = poorly differentiated.}
#'      \item{\code{Age}}{in years at time of diagnosis.}
#'      \item{\code{Cond}}{condition: 1 = no disability, 2 = restricted work, 3 = requires assistance with self care, 4 = bed confined.}
#'      \item{\code{Site}}{1 = faucial arch, 2 = tonsillar fossa, 3 = posterior pillar, 4 = pharyngeal tongue, 5 = posterior wall.}
#'      \item{\code{T}}{T staging: 1 = primary tumor measuring 2 cm or less in largest diameter; 2 = primary tumor measuring 2 to 4 cm in largest diameter, minimal infiltration in depth; 3 = primary tumor measuring more than 4 cm; 4 = massive invasive tumor.}
#'      \item{\code{N}}{N staging: 0 = no clinical evidence of node metastases; 1 = single positive node 3 cm or less in diameter, not fixed; 2 = single positive node more than 3 cm in diameter, not fixed; 3 = multiple positive nodes or fixed positive nodes.}
#'      \item{\code{EntryDate}}{Date of entry: Day of year and year.}
#'      \item{\code{Status}}{0 = censored, 1 = dead.}
#'      \item{\code{Time}}{in days from day of diagnosis.}
#'  }
#'
#' @references  Kalbfleisch, J. D. and Prentice, R. L. (2002) \emph{The Statistical Analysis of Failure Time Data}. John Wiley &  Sons, New York, 2nd edition.
#'
#' @keywords datasets
#'
"tonsil"
