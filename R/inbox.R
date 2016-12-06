#' Defining Subset lying in a Box
#'
#' @description Help function that returns the subset of the whole data set which is contained in one box (\code{fixbox}). This function is needed in other functions like \code{\link{PRIM_peel_bs}}.
#'
#' @param X \code{\link{data.frame}} from which the subset should be taken.
#' @param fixbox_metric one row of a \code{box_metric} from a "\code{peel}"-object.
#' @param fixbox_nom one element of a \code{box_nom} list from a "\code{peel}"-object.
#' @param fixbox_na one element of a \code{box_na} list from a "\code{peel}"-object.
#'
#' @details The arguments \code{fixbox_metric}, \code{fixbox_nom}, \code{fixbox_na} can also come from a "\code{fixbox}"-object.
#'
#'
#' @return \code{inbox} returns a logical vector indicating the observations that lie in the box.
#'
#' @export


inbox <- function(X, fixbox_metric = NULL, fixbox_nom = NULL, fixbox_na = NULL) {
  inbox_metric_index <- function(X, box) {
    apply(sapply(2:ncol(X), function(k){
      a <- X[,k] > box[,k-1] & X[,k] < box[,(ncol(X)-1)+k-1]
      a[is.na(a)] <- TRUE
      a
    }),1,all )
  }
  inbox_nom_index <- function(X, box){
    apply(sapply(1:length(box), function(k){
      a <- X[, match(names(box)[k], names(X))]  %in% names(box[[k]])[as.vector(t(box[[k]]))]
      a[is.na(X[, match(names(box)[k], names(X))])] <- TRUE
      a
    }),1,all )
  }
  inbox_na_index <- function(X, box){
    apply(sapply(1:length(box), function(k) {if(box[k]==FALSE) a <- !is.na(X[, match(names(box)[k], names(X))])  else a <- rep(TRUE, nrow(X)); a}), 1, all)
  }

  metric <- sapply(1:ncol(X), function(k) is.numeric(X[,k]) | is.ordered(X[,k]) | is.logical(X[,k]))
  if(sum(metric[-1])!=0) ind_metric <- inbox_metric_index(X[,metric], fixbox_metric) else ind_metric <- TRUE
  if(!is.null(fixbox_nom)) ind_nom <- inbox_nom_index(X, fixbox_nom) else ind_nom <- TRUE
  if(!is.null(fixbox_na)) ind_na <- inbox_na_index(X, fixbox_na) else ind_na <- TRUE
  ind <- ind_metric & ind_nom & ind_na
  ind
}

