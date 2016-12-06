#' Relative Frequency Plots
#'
#' @description This Function creates the relative frequency plots of a box defined by \code{\link{define_fixbox}}.
#'
#' @param fixbox object of class "\code{fixbox}" to be used to create the relative frequency plots.
#' @param plot_together logical. If \code{TRUE} all plots are drawn together in one figure by \code{par(mfrow = ...)}.
#' @param ... further arguments of the functions \code{\link{plot}} and \code{\link{barplot}}.
#'
#' @details Relative frequency plots are a diagnostic tool for a result of PRIM (\code{fixbox}).
#' These are histograms or barplots illustrating ratios of distributions \eqn{r_{jk}} for each variable \eqn{x_j}.
#' This ratio is defined by the distribution of \eqn{x_j} in one box \eqn{B_k} divided by the distribution of this variable over the whole data set.
#'
#'
#' @seealso \code{\link{define_fixbox}}, \code{\link{PRIM}}
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
#' # apply the PRIM function to find the best box with a support of at least 0.1:
#' p <- PRIM(y~., data=dat, beta_min = 0.1, max_boxes = 1, print_position = FALSE)
#'
#' # relative frequency plots:
#' relfreq(p$fixboxes[[1]])
#'
#' @export


# Function to create the relative frequency plots of a (fix-)box
relfreq <- function(fixbox, plot_together=TRUE, ...){
  X <- fixbox$data_orig
  if(is.Surv(X[,1])) X <- X[,-1]
  if(plot_together) par(mfrow=c(ceiling(sqrt(ncol(X))), ceiling(sqrt(ncol(X)))))
  for (i in 1:ncol(X)){
    if(is.factor(X[,i]) | is.logical(X[,i])){
      t <- table(X[,i])/nrow(X)
      t2 <- table(X[fixbox$subset,i])/sum(fixbox$subset)
      relfreq <- ifelse(t==0, 0, t2/t)
      barplot(relfreq, ylab="rel. freq.", main=names(X)[i], ...)
    } else{
      h <- hist(X[,i], plot = FALSE)
      h2 <- hist(X[fixbox$subset,i], breaks=h$breaks, plot = FALSE)
      h$density <- ifelse(h$density==0, 0, h2$density/h$density)
      plot(h, freq = F, ylab="rel. freq.", col = "grey", main=names(X)[i], xlab="", ...)
    }
  }
  if(plot_together) par(mfrow=c(1,1))
}

