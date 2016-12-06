#' Interbox Dissimilarity
#'
#' @description This Function calculates the interbox dissimilarity between two boxes defined on one data set.
#'
#' @param fixbox1 first object of class "\code{fixbox}".
#' @param fixbox2 second object of class "\code{fixbox}".
#' @param data original data used to calculate the supports needed for the interbox dissimilarity. If this argument is missing \code{data_orig} from \code{fixbox1} is used as \code{data}.
#'
#' @details The interbox dissimilarity is a diagnostic tool of PRIM measuring the dissimilarity between two boxes \eqn{B_k} and \eqn{B_l}.
#' It is defined as the difference between the smallest box \eqn{B_{kl}} that covers both boxes and the support of their union.
#'
#' The interbox dissimilarity can assume values between 0 and 1.
#' While nested boxes have a dissimilarity of 0, it gets bigger the more different the two boxes are.
#'
#'
#'
#' @export


inter_diss <- function(fixbox1, fixbox2, data){
  if(missing(data)) data <- fixbox1$data_orig

  if(!is.null(fixbox1$fixbox_metric)){
    b <- rbind(fixbox1$fixbox_metric, fixbox2$fixbox_metric)
    bigbox_metric <- t(as.data.frame(c(apply(b[,1:(ncol(b)/2), drop=FALSE], 2, min), apply(b[,(ncol(b)/2+1):ncol(b), drop=FALSE], 2, max))))
    row.names(bigbox_metric) <- 1
  } else bigbox_metric <- NULL

  if(!is.null(fixbox1$fixbox_nom)){
    bigbox_nom <- lapply(1:length(fixbox1$fixbox_nom), function(k) fixbox1$fixbox_nom[[k]] | fixbox2$fixbox_nom[[k]])
    names(bigbox_nom) <- names(fixbox1$fixbox_nom)
  } else bigbox_nom <- NULL

  if(!is.null(fixbox1$fixbox_na)){
    bigbox_na <- sapply(1:length(fixbox1$fixbox_na), function(k) fixbox1$fixbox_na[[k]] | fixbox2$fixbox_na[[k]])
    names(bigbox_na) <- names(fixbox1$fixbox_na)
  } else bigbox_na <- NULL

  beta_bb <- mean(inbox(data, fixbox_metric = bigbox_metric, fixbox_nom = bigbox_nom, fixbox_na = bigbox_na))

  beta_v <- mean(inbox(data, fixbox_metric = fixbox1$fixbox_metric, fixbox_nom = fixbox1$fixbox_nom, fixbox_na = fixbox1$fixbox_na) | inbox(data, fixbox_metric = fixbox2$fixbox_metric, fixbox_nom = fixbox2$fixbox_nom, fixbox_na = fixbox2$fixbox_na))

  return(beta_bb - beta_v)
}


