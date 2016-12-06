#' Defining a "fixbox"-Object
#'
#' @description Function for determining on one box out of a "\code{peel}"-object.
#'
#' @param prim object of class "\code{peel}" from which the box should be defined.
#' @param step number of the step (i.e. row in the argument \code{box} of the "\code{peel}"-object.) of the requested box.
#'
#' @return \code{define_fixbox} returns a object of class "\code{fixbox}", which is the basically same as a "\code{peel}"-object, but only contains one box.
#' It is a list containing at least the following elements:
#'
#' @return \item{f}{target function of the selected box.}
#' @return \item{beta}{support of the selected box.}
#' @return \item{box}{a \code{\link{data.frame}} with one row defining the borders of the box. The columns with "\code{min.}" and "\code{max.}" describe the lower and upper boundaries of the at least ordinal covariates. Therefore the value taken is the last one that is \bold{not} included in the box.
#'
#' For the nominal variables there are columns for every category they can take. If the category is removed from the box the value \code{FALSE} is taken. The names of these columns are structured like: \code{<variable name>.<category>}
#'
#' For each variable with missing values (only if \code{use_NAs = TRUE}) there is also a column taking the value \code{FALSE} if the \code{NA}s of this variable are removed from the box. The names of these columns are structured like: \code{<variable name>.NA}
#' }
#' @return \item{box_metric, box_nom, box_na}{easier to handle definitions of the box for other functions}
#' @return \item{subset}{logical vector indicating the subset (i.e. the observations that lie in the box).}
#' @return \item{data_orig}{original dataset that was used for the peeling.}
#'
#' @export

define_fixbox <- function(prim, step){
  f <- prim$f[step]
  beta <- prim$beta[step]
  fixbox <- prim$box[step,]
  fixbox_metric <- prim$box_metric[step,]
  subset <- prim$subsets[[step]]
  if(!is.null(prim$box_nom[[1]])) fixbox_nom <- prim$box_nom[[step]] else fixbox_nom <- NULL
  if(!is.null(prim$box_na[[1]])) fixbox_na <- prim$box_na[[step]] else fixbox_na <- NULL

  ret <- list(f=f, beta=beta, fixbox=fixbox, fixbox_metric=fixbox_metric, fixbox_nom=fixbox_nom, fixbox_na=fixbox_na, subset=subset, data_orig=prim$data_orig, target=prim$target)
  class(ret) <- "fixbox"
  return(ret)
}


