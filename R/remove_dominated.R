#' Remove Dominated Boxes from a Peeling Output
#'
#' @import survival
#'
#' @description This function removes all boxes from one "\code{peel}"-object that are dominated by another box.
#' I.e. there is at least one other box with a better output and a bigger support beta.
#'
#'
#' @param prim object of class "\code{peel}".
#' @param sort logical. If \code{TRUE} the boxes in the output are sorted decreasing by their box supports.
#' @param dup.rm logical. Tf \code{TRUE} duplicated boxes in the output are removed. Boxes are duplicated if another box has exactly the same target and support.
#'
#' @details Dominated boxes mainly occur in outputs of the multiple peeling function (\code{\link{PRIM_peel_bs}}).
#' So this function is practically only useful for "\code{peel}"-objects beeing the results of such a function.
#' Without the dominated boxes the object gets much smaller by keeping all the relevant boxes.
#' Also the trajectory (see \code{\link{plot.peel}}) without the not useful dominated boxes looks much clearer.
#'
#' @return This function returns another "\code{peel}"-object having the same structure as the input object, just without the dominated boxes.
#'
#' @seealso \code{\link{PRIM_peel_bs}}, \code{\link{plot.peel}}
#'
#' @export


remove_dominated <- function(prim, sort=TRUE, dup.rm=TRUE){

  # which observations are dominated by at least one other
  dominated <- sapply(1:length(prim$f), function(k) any(prim$f[-k] >= prim$f[k] & prim$beta[-k] >= prim$beta[k] &!
                                                        (prim$f[-k] == prim$f[k] & prim$beta[-k] == prim$beta[k])))

  # remove dominated boxes from the result
  result <- list(f=prim$f[!dominated], beta=prim$beta[!dominated], box=prim$box[!dominated,], box_metric=prim$box_metric[!dominated,], box_nom=prim$box_nom[!dominated], box_na=prim$box_na[!dominated], subsets=prim$subsets[!dominated], data_orig=prim$data_orig, target=prim$target)

  # remove boxes having the same target and support than one other (duplicates)
  if(dup.rm==TRUE){
    dup <- duplicated(cbind(result$f,result$beta))
    result <- list(f=result$f[!dup], beta=result$beta[!dup], box=result$box[!dup,], box_metric=result$box_metric[!dup,], box_nom=result$box_nom[!dup], box_na=result$box_na[!dup], subsets=result$subsets[!dup], data_orig=prim$data_orig, target=prim$target)
  }

  # sort the boxes decreasing by beta
  if(sort==TRUE){
    ord <- order(result$beta, decreasing = TRUE)
    result$f <- result$f[ord]
    result$beta <- result$beta[ord]
    result$subsets <- result$subsets[ord]
    result$box <- result$box[ord,]
    result$box_metric <- result$box_metric[ord,]
    result$box_nom <- result$box_nom[ord]
    result$box_na <- result$box_na[ord]
  }

  class(result) <- "peel"
  return(result)
}

