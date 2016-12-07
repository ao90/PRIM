#' Pasting-Function
#'
#' @import survival
#'
#' @description This function is an implementation of the Pasting-Algorithm as suggested by Friedman and Fisher (1999). In each iteration the fraction alpha is pasted to one edge of the current box.
#'
#' @param fixbox an object of class \code{fixbox}, which was defined after the peeling function and now should be used for pasting.
#' @param paste_alpha alpha-fraction that is pasted to the box at each iteration.
#' @param max_steps maximum number of pasting steps the function should make.
#' @param stop_by_dec logical. If \code{TRUE} the pasting stops if the target at one step is lower than the target of the last step.
#'
#' @details The outcome of this function is also a "\code{peel}"-object, because it has basically the same structure as the outcome of the peeling functions. The only difference is, that pasting goes from small supports to bigger ones, while by peeling its the other way round.
#'
#' @return \code{PRIM_paste} returns an object of class "\code{peel}", which is a list containing at least the following components:
#' @return \item{f}{vector of the target functions evaluated on the box at each pasting step.}
#' @return \item{beta}{vector of the supports beta of the boxes at each pasting step.}
#' @return \item{box}{a \code{\link{data.frame}} defining the borders of the boxes. Each row belongs to one pasting step. The columns with "\code{min.}" and "\code{max.}" describe the lower and upper boundaries of the at least ordinal covariates. Therefore the value taken is the last one that is \bold{not} included in the current box.
#'
#' For the nominal variables there are columns for every category they can take. If the category is removed from the box the value \code{FALSE} is taken. The names of these columns are structured like: \code{<variable name>.<category>}
#'
#' For each variable with missing values (only if \code{use_NAs = TRUE}) there is also a column taking the value \code{FALSE} if the \code{NA}s of this variable are removed from the current box. The names of these columns are structured like: \code{<variable name>.NA}
#' }
#' @return \item{box_metric, box_nom, box_na}{easier to handle definitions of the boxes for other functions}
#' @return \item{subsets}{\code{list} of logical vectors indicating the subsets at each pasting step (i.e. the observations that lie in the box)}
#' @return \item{data_orig}{original dataset that is used (extracted from \code{fixbox}).}
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
#' plot(dat$x1, dat$x2, col=dat$y+1, pch=16)
#' remove(x1, x2, x3, y, wsk, cat, n)
#'
#' # apply the PRIM_peel function:
#' prim <- PRIM_peel(y ~ ., data = dat, beta_min = .01, peel_alpha = .1)
#' plot(prim)
#' abline(h=prim$f[17], v=prim$beta[17]) # box decided to paste
#' fix <- define_fixbox(prim, 17) # define fixbox
#'
#' # apply the PRIM_paste function:
#' paste <- PRIM_paste(fix, stop_by_dec = FALSE)
#' head(cbind(paste$box, paste$f, paste$beta))
#'
#' @export


PRIM_paste <- function(fixbox, paste_alpha=0.01, max_steps=50, stop_by_dec=TRUE){

  # help function that returns the next higher or lower border of a "box", so that at least N new observations are included
  next_value <- function(vec, border, N = 1, lower = FALSE){
    if(lower) {vec <- vec * (-1); border <- -border}
    sort_vec <- sort(vec)
    if(border > sort_vec[length(sort_vec)]) if(lower) return(-Inf) else return(Inf)
    if(!(border %in% vec)) border <- vec[which.min(ifelse(vec-border > 0, vec-border, Inf))]
    if(lower) res <- -sort_vec[sort_vec > border][N] else res <- sort_vec[sort_vec > border][N]
    if(lower & is.na(res)) res <- -Inf
    if(!lower & is.na(res)) res <- Inf
    return(res)
  }

  # define start values
  target <- fixbox$target
  data <- fixbox$data_orig
  obs_in_box <- list(fixbox$subset)
  f <-  fixbox$f
  beta <- fixbox$beta
  box_metric <- fixbox$fixbox_metric
  if(!is.null(box_metric)) row.names(box_metric) <- 1
  box_nom <- list(fixbox$fixbox_nom)
  box_na <- list(fixbox$fixbox_na)

  # which covariables are metric; which variables are nominal
  metric_cov <- ifelse(!is.null(fixbox$fixbox_metric), ncol(fixbox$fixbox_metric)/2, 0)
  nom_var <- sapply(1:ncol(fixbox$data_orig), function(k) is.factor(fixbox$data_orig[,k]))

  # types of "variables" defining the subboxes in one iteration
  types <- rep("metric", 2*metric_cov)
  if(!is.null(box_nom[[1]])) types <- c(types, rep("nom", sum(sapply(box_nom[[1]], length))))
  types <- c(types, rep("na", length(box_na[[1]])))

  # loop
  for(i in 1:max_steps){
    # define the current data.frame (i.e. observations in the box)
    sub <- data[obs_in_box[[i]],]

    # generate all candidate subsets resulting from the ordinal and skalar variables
    if(!is.null(fixbox$fixbox_metric)){
      help_box <- lapply(1:metric_cov, function(k) box_metric[i,])
      help_box <- lapply(1:metric_cov, function(k) {help_box[[k]][1,k] <- -Inf; help_box[[k]]})
      lower_boxes <- lapply(1:metric_cov, function(k) {
        b <- box_metric[i,]
        b[1,k] <- next_value(vec=data[inbox(data, help_box[[k]], box_nom[[length(box_nom)]], box_na[[length(box_na)]]),][, !nom_var, drop=FALSE][,k+1], border=as.numeric(box_metric[i,k]), N=ceiling(paste_alpha*nrow(sub)), lower=TRUE)
        b
        })
      subsets_low <- lapply(1:metric_cov, function(k) inbox(data, lower_boxes[[k]], box_nom[[length(box_nom)]], box_na[[length(box_na)]]))

      help_box <- lapply(1:metric_cov, function(k) box_metric[i,])
      help_box <- lapply(1:metric_cov, function(k) {help_box[[k]][metric_cov+k] <- Inf; help_box[[k]]})
      higher_boxes <- lapply(1:metric_cov, function(k) {
        b <- box_metric[i,]
        b[1,metric_cov+k] <- next_value(vec=data[inbox(data, help_box[[k]], box_nom[[length(box_nom)]], box_na[[length(box_na)]]),][, !nom_var, drop=FALSE][,k+1], border=as.numeric(box_metric[i,metric_cov+k]), N=ceiling(paste_alpha*nrow(sub)), lower=FALSE)
        b
      })
      subsets_high <- lapply(1:metric_cov, function(k) inbox(data, higher_boxes[[k]], box_nom[[length(box_nom)]], box_na[[length(box_na)]]))
    } else subsets_low <- subsets_high <- NULL

    # generate all candidate subsets resulting from the nominal variables
    nom_boxes <- NULL
    if(!is.null(fixbox$fixbox_nom)){
      for(j in 1:length(box_nom[[i]])){
        nom_boxes <- c(nom_boxes, lapply(1:length(box_nom[[i]][[j]]), function(k) {a<-box_nom[[i]]; a[[j]][k]<-TRUE; a}))
      }
      subsets_nom <- lapply(1:length(nom_boxes), function(k) inbox(data, box_metric[i,], nom_boxes[[k]], box_na[[length(box_na)]]))
    } else subsets_nom <- NULL

    # generate all candidate subsets resulting from the missing values of variables
    na_boxes <- NULL
    if(!is.null(fixbox$fixbox_na)){
      na_boxes <- lapply(1:length(box_na[[i]]), function(k) {a<-box_na[[i]]; a[k]<-TRUE; a})
      subsets_na <- lapply(1:length(na_boxes), function(k) inbox(data, box_metric[i,], box_nom[[length(box_nom)]], na_boxes[[k]]))
    } else subsets_na <- NULL



    # put all candidate subsets in one list object
    subsets <- c(subsets_low, subsets_high, subsets_nom, subsets_na)

    if(!is.null(fixbox$fixbox_metric)) b <- c(lower_boxes, higher_boxes)

    # evaluate the target function on all data sets
    t <- sapply(subsets, function(k) target(data[k,1]))
    # if one box limit is already infinite, the target function therefor is set to NA
    # in this case no further data point could be added
    if(!is.null(fixbox$fixbox_metric)) t[types=="metric"][!is.finite(as.numeric(box_metric[i,]))] <- NA
    # if one category of a nominal variable is already included, the target function of this subset becomes NA
    t[types=="nom"][unlist(box_nom[[length(box_nom)]]) == TRUE] <- NA
    # if the NAs of one variable are already included, the target function of this subset becomes NA
    t[types=="na"][unlist(box_na[[length(box_na)]]) == TRUE] <- NA


    # find out which data set maximizes the target function
    which_max <- which.max(t)

    # end the loop if there is no maximum (i.e. no observations left) or if the new target is smaller than the last one (if stop_by_dec=TRUE)
    if(length(which_max) == 0) break
    if(stop_by_dec & t[which_max] < f[i]) break

    # update the next value of all outputs
    obs_in_box[[i+1]]  <- subsets[[which_max]]
    f[i+1] <- t[which_max]
    beta[i+1] <- sum(obs_in_box[[i+1]])/nrow(data)
    if (metric_cov!=0) if((types=="metric")[which_max]) box_metric[i+1,] <- b[[which_max]] else box_metric[i+1,] <- box_metric[i,]
    if((types=="nom")[which_max]) box_nom[[i+1]] <- nom_boxes[[which_max-sum(types=="metric")]] else box_nom[[i+1]] <- box_nom[[length(box_nom)]]
    if((types=="na")[which_max]) box_na[[i+1]] <- na_boxes[[which_max-sum(types %in% c("metric", "nom"))]] else box_na[[i+1]] <- box_na[[length(box_na)]]
    }

  # create the element "box" with all box limits in one data.frame
  if(!is.null(fixbox$fixbox_nom)){
    b_n <- as.data.frame(matrix(unlist(box_nom), ncol = length(unlist(box_nom[[1]])), byrow = TRUE))
    names(b_n) <- names(unlist(box_nom[[1]]))
  } else b_n <- NULL
  if(!is.null(fixbox$fixbox_na)){
    b_na <- as.data.frame(matrix(unlist(box_na), ncol = length(unlist(box_na[[1]])), byrow = TRUE))
    names(b_na) <- paste("NA", names(unlist(box_na[[1]])), sep = ".")
  } else b_na <- NULL

  if((!is.null(box_metric)) & (!is.null(b_n))) box <- cbind(box_metric, b_n) else
    if(!is.null(box_metric)) box <- box_metric else box <- b_n
  if(!is.null(b_na)) box <- cbind(box, b_na)

  # define the function output
  ret <- list(f=f, beta=beta, box=box, box_metric=box_metric, box_nom=box_nom, box_na=box_na, subsets=obs_in_box, data_orig=fixbox$data_orig, target=target)
  class(ret) <- "peel"
  return(ret)
}


