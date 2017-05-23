#' Peeling-Function
#'
#' @import survival
#'
#' @description This function is an implementation of the (singular) Peeling-Algorithm as suggested by Friedman and Fisher (1999). In each iteration the fraction alpha is peeled from one edge of the current box.
#'
#' @param formula an object of class "\code{\link{formula}}" with a response but no interaction terms.
#' It indicates the response over which the target function should be maximized and the covariates that are used for the later box definitions.
#' If this argument is missing, the argument data is used as the \code{\link{model.frame}}.
#' @param data an object of class \code{\link{data.frame}} containing the variables named in the formula.
#' @param peel_alpha alpha-fraction used for the peeling.
#' @param beta_min minimum support that one Box should have (stop-criterion).
#' @param target target-function to be maximized.
#' @param alter_crit logical. If \code{TRUE} the alternative criterion is used for peeling.
#' @param use_NAs logical. If \code{TRUE} observations with missing values are included in the analysis.
#'
#' @details The outcome of the \code{formula} can either be numeric, logical or a survival object (see \code{\link{Surv}}). If it is a survival object the \code{target} is set to the number of events per amount of time.
#'
#' @details This function is the main part of the multiple Version \code{PRIM_peel_bs}.
#'
#' @return \code{PRIM_peel} returns an object of class "\code{peel}", which is a list containing at least the following components:
#' @return \item{f}{vector of the target functions evaluated on the box at each peeling step.}
#' @return \item{beta}{vector of the supports beta of the boxes at each peeling step.}
#' @return \item{box}{a \code{\link{data.frame}} defining the borders of the boxes. Each row belongs to one peeling step. The columns with "\code{min.}" and "\code{max.}" describe the lower and upper boundaries of the at least ordinal covariates. Therefore the value taken is the last one that is \bold{not} included in the current box.
#'
#' For the nominal variables there are columns for every category they can take. If the category is removed from the box the value \code{FALSE} is taken. The names of these columns are structured like: \code{<variable name>.<category>}
#'
#' For each variable with missing values (only if \code{use_NAs = TRUE}) there is also a column taking the value \code{FALSE} if the \code{NA}s of this variable are removed from the current box. The names of these columns are structured like: \code{<variable name>.NA}
#' }
#' @return \item{box_metric, box_nom, box_na}{easier to handle definitions of the boxes for other functions}
#' @return \item{subsets}{\code{list} of logical vectors indicating the subsets at each peeling step (i.e. the observations that lie in the box)}
#' @return \item{data_orig}{original dataset that is used for the peeling.}
#'
#'
#' @references Friedman, J. H. and Fisher, N. I., 'Bump hunting in high-dimensional data', Statistics and Computing \bold{9} (2) (1999), 123-143
#'
#' @examples
#' # generating random data:
#' set.seed(123)
#' n <- 500
#' x1 <- runif(n = n, min = -1)
#' x2 <- runif(n = n, min = -1)
#' x3 <- runif(n = n, min = -1)
#' cat <- as.factor(sample(c("a","b","c", "d"), size = n, replace = TRUE))
#' wsk <- (1-sqrt(x1^2+x2^2)/sqrt(2))
#' y <- as.logical(rbinom(n = n, prob = wsk, size = 1))
#' dat <- cbind.data.frame(y, x1, x2, x3, cat)
#' #plot(dat$x1, dat$x2, col=dat$y+1, pch=16)
#' remove(x1, x2, x3, y, wsk, cat, n)
#'
#' # apply the PRIM_peel function:
#' prim <- PRIM_peel(formula=y ~ ., data=dat, beta_min = .01, peel_alpha = .1)
#' plot(prim) # trajectory
#' prim$box # box definitions
#'
#' @export

PRIM_peel <- function(formula, data, peel_alpha = .05, beta_min = .01, target = mean, alter_crit = TRUE, use_NAs = TRUE){

  # generate the data.frame X the function is working with
  if(use_NAs) {if(!missing(formula)) X <- model.frame(formula, data=data, na.action = na.pass) else X <- data}
  if(use_NAs==FALSE) {if(!missing(formula)) X <- model.frame(formula, data=data, na.action = na.omit) else X <- na.omit(data)}

  # use the events per time as target, if the response is a survival object
  if(is.Surv(X[,1])) target <- function(k) sum(k[,2])/sum(k[,1])

  # mark the ordinal covariates and set them into numeric
  ord_cov <- sapply(1:ncol(X),function(k) is.ordered(X[,k]))[-1]
  for (i in 1:ncol(X)) if (is.ordered(X[,i])) X[,i] <- as.numeric(X[,i])

  # mark nominal variables
  nom_var <- sapply(1:ncol(X),function(k) is.factor(X[,k]))

  # which covartiates from all non nominal covariates are ordinal
  ord_cov_short <- ord_cov[!nom_var[-1]]

  # which variables have at least one missing value
  na_vars <- apply(apply(X, 2, is.na),2,any)

  if(na_vars[1]) {
    X <- X[!is.na(X[,1]),]
    na_vars[1] <- FALSE
  }

  # define start values
  beta <- 1
  f <- target(X[,1])
  sub <- list(X)
  crit <- NULL
  obs_in_box <- list(rep(TRUE, times=nrow(X)))
  if(sum(!nom_var[-1])>0){
    box <- rbind.data.frame(c(rep(-Inf, times=sum(!nom_var[-1])), rep(Inf, times=sum(!nom_var[-1]))))
    names(box) <- paste(rep(c("min", "max"), each=ncol(X[,!nom_var])-1), names(X[,!nom_var])[-1], sep=".")
  }
  if (sum(nom_var)>0){
    b_n <- rbind.data.frame(rep(TRUE, times=sum(sapply(1:ncol(subset(X, select=nom_var)), function(k) length(levels(subset(X, select=nom_var)[,k]))))))
    names(b_n) <- paste(rep(names(subset(X, select=nom_var)), times=sapply(1:ncol(subset(X, select=nom_var)), function(k) length(levels(subset(X, select=nom_var)[,k])))), unlist(lapply(1:ncol(subset(X, select=nom_var)), function(k) levels(subset(X, select=nom_var)[,k]))), sep=".")

    if(sum(!nom_var[-1])>0) box <- cbind(box, b_n) else box <- b_n
  }
  nomcols_in_box <- ifelse(sum(nom_var)>0, ncol(b_n), 0)
  if (sum(na_vars)>0){
    b_na <- rbind.data.frame(rep(TRUE, times=sum(na_vars)))
    names(b_na) <- paste(names(na_vars)[na_vars==TRUE], "NA", sep = ".")
    box <- cbind(box, b_na)
  }

  # loop
  i <- 1
  while(TRUE){
    # define the current data.frame (i.e. observations in the box)
    sub <- X[obs_in_box[[i]], ]

    # generate all candidate subsets resulting from the ordinal and skalar variables
    if(sum(!nom_var[-1])>0){
      sub_metric <- sub[,!nom_var, drop=FALSE][,-1, drop=FALSE]
      box_mins <- apply(sub_metric, 2, function(k) quantile(k, probs = peel_alpha, na.rm = TRUE))
      box_maxs <- apply(sub_metric, 2, function(k) quantile(k, probs = 1 - peel_alpha, na.rm = TRUE))
      subsets_low <- lapply(1:length(box_mins), function(k) sub_metric[,k] > box_mins[k])
      subsets_high <- lapply(1:length(box_maxs), function(k)  subset = sub_metric[,k] < box_maxs[k])
      # replace NAs with TRUE, so that obs with NAs stay in the current box
      subsets_low <- lapply(subsets_low, function(k) {k[is.na(k)] <- TRUE; k})
      subsets_high <- lapply(subsets_high, function(k) {k[is.na(k)] <- TRUE; k})

   } else subsets_low <- subsets_high <- NULL

    # generate all candidate subsets resulting from the nominal variables
    subsets_nom <- NULL
    if (sum(nom_var)>0){
      sub_nom <- subset(sub, select=nom_var)
      for (j in 1:sum(nom_var)){
        subsets_nom <- c(subsets_nom, lapply(levels(sub_nom[,j]), function(k)  sub_nom[,j]!=k))
      }
      # replace NAs with TRUE, so that obs with NAs stay in the current box
      subsets_nom <- lapply(subsets_nom, function(k) {k[is.na(k)] <- TRUE; k})
    }

    # generate all candidate subsets resulting from the missing values of variables
    subsets_na <- NULL
    if (sum(na_vars)>0){
      sub_na <- subset(sub, select=na_vars)
      subsets_na <- lapply(1:sum(na_vars), function(k) !is.na(sub_na[,k]))
    }

    # put all candidate subsets in one list object
    subsets <- c(subsets_low, subsets_high, subsets_nom, subsets_na)

    # put all possible updates of the box limits in one list
    if(sum(!nom_var[-1])>0) update <- as.list(c(box_mins, box_maxs)) else update <- NULL
    if (sum(nom_var)>0) nom_update <- as.list(rep(FALSE, times=nomcols_in_box)) else nom_update <- NULL
    if (sum(na_vars)>0) na_update <- as.list(rep(FALSE, times=ncol(b_na))) else na_update <- NULL
    b <- c(update, nom_update, na_update)

    # evaluate the target function on all data sets
    t <- sapply(subsets, function(k) target(sub[,1][k]))
    # if one category of a nominal variable is already removed, the target function of this subset becomes NA
    t[box[i,] == FALSE & sapply(1:ncol(box), function(k) is.logical(box[,k]))] <- NA

    # alternative criterion
    if(alter_crit) cr <- (t - f[i]) / (nrow(sub) - sapply(subsets, sum)) else cr <- t

    # find out which data set maximizes the target function
    which_max <- which.max(cr)


    # end the loop if there is no maximum (i.e. no observations left) or the support falls below beta_min
    if(length(which_max) == 0) break
    if(sum(subsets[[which_max]])/nrow(X) < beta_min) break

    # update the next value of all outputs
    obs_in_box[[i+1]] <- obs_in_box[[i]]
    obs_in_box[[i+1]][obs_in_box[[i+1]]]  <- subsets[[which_max]]
    f[i+1] <- t[which_max]
    crit[i+1] <- cr[which_max]
    beta[i+1] <- sum(obs_in_box[[i+1]])/nrow(X)
    box[i+1,] <- box[i,]
    box[i+1, which_max] <- b[[which_max]]
    if(sum(obs_in_box[[i+1]])<=1) break
    i<-i+1
  }

  # construct a data frame, which contains only the box limits of the "metric" (i.e. skalar or ordinal) Variables
  if(sum(!nom_var[-1])>0){
    box_metric <- box[, 1:(2*sum(!nom_var[-1]))]
  } else box_metric <- NULL

  # construct a list, which contains only the box limits of the nominal Variables
  if (sum(nom_var)>0){
   if(sum(!nom_var[-1])>0) box_n <- box[, (2*sum(!nom_var[-1])+1):(2*sum(!nom_var[-1])+nomcols_in_box)] else box_n <- box
   lev <-  c(1, sapply(1:sum(nom_var), function(k) length(levels(subset(X, select=nom_var)[,k]))))
   box_nom <- lapply(1:nrow(box_n), function(k) lapply(1:(length(lev)-1), function(l) setNames(as.vector(t(box_n[k, sum(lev[1:l]):(sum(lev[1:(l+1)])-1) ])), levels(subset(X, select=nom_var)[,l]))))
   for(j in 1:length(box_nom)) names(box_nom[[j]]) <- names(X)[nom_var]
  } else box_nom <- NULL

  # construct a data.frame, which contains only the box limits given by missing values
  if (sum(na_vars)>0){
    b <- box[, (2*sum(!nom_var[-1])+nomcols_in_box+1):ncol(box), drop=FALSE]
    box_na <- lapply(1:nrow(b), function(k) {
      a <- as.vector(t(b[k,]))
      names(a) <- names(na_vars[na_vars])
      a
      })
  } else box_na <- NULL

  # define the function output
  ret <- list(f=f, beta=beta, box=box, box_metric=box_metric, box_nom=box_nom, box_na=box_na, subsets=obs_in_box, data_orig=X, target=target, nom=nom_var)
  class(ret) <- "peel"
  return(ret)
}

