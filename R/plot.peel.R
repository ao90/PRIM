#' Plotting a Trajectory
#'
#'
#' @description This funcion plots the trajectory of a "\code{peel}"-object. In the trajectory the target functions evaluated on the data in the box are plotted against the supports beta.
#'
#' @param x object of class "\code{peel}".
#' @param ... further arguments of the \code{\link{plot}}-function.
#'
#' @export


plot.peel <- function(x, ...){
  f <- x$f
  beta <- x$beta
  plot(f ~ beta, ...)
}
