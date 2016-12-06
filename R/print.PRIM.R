#' Printing a "PRIM"-Object
#'
#'
#' @description This function prints a \code{PRIM}-object. It is a method for the generic function \code{print} of class "\code{PRIM}".
#'
#' @param x object of class "\code{PRIM}".
#' @param ... further arguments of the \code{\link{print}}-function.
#'
#' @export


print.PRIM <- function(x, ...){
  
  for(i in 1:(length(x$f)-1)){
    cat("Box ", i, ":", sep = "")
    cat("   target =", x$f[i])
    cat("   support =", x$beta[i])
    
    cat("\nDefinition:\n\n")
    
    if(!is.null(x$box_metric)){
      cat("Quantitative variables:\n")
      print(data.frame(lower=t(x$box_metric)[1:(ncol(x$box_metric)/2),i], upper=t(x$box_metric)[(ncol(x$box_metric)/2+1):ncol(x$box_metric),i], row.names = substring(names(x$box_metric[1,]), first = 5)[1:(ncol(x$box_metric)/2)]))
      cat("\n")
    }
    
    if(!is.null(x$box_nom)){
      cat("Qualitative variables:\n")
      for(j in 1:length(x$box_nom[[i]])) {
        cat(paste(names(x$box_nom[[i]][j]), ":\n", sep = ""))
        print(x$box_nom[[i]][[j]])
      }
      cat("\n")
    }
      
    if(!is.null(x$box_na)){
      cat("Missing values of variables:\n")
      print(x$box_na[[i]])
    }
    
    cat("\n-------------------------------------------\n\n")
  }
  
  
  cat("Data not included in a box:")
  cat("\ntarget =", x$f[length(x$f)])
  cat("\nsupport =", x$beta[length(x$f)])
  
}

