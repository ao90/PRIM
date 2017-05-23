#' Multiple Peeling-Function
#'
#' @import survival
#'
#' @description This function is an implementation of the multiple Peeling-Algorithm as suggested by Friedman and Fisher (1999). The singular peeling function \code{\link{PRIM_peel}} is repeated for different alpha's and bootstrap samples out of the original data.
#'
#'
#' @param formula an object of class "\code{\link{formula}}" with a response but no interaction terms.
#' It indicates the response over which the target function should be maximized and the covariates that are used for the later box definitions.
#' @param data an object of class \code{\link{data.frame}} containing the variables named in the formula.
#' @param peel_alpha vector of a sequence of different alpha-fractions used for the peelings.
#' @param B number of bootstrap samples on which the peeling is applied to for each alpha. For \code{B = 0} no bootstraps are created.
#' @param beta_min minimum support that one Box should have (stop-criterion).
#' @param target target-function to be maximized.
#' @param alter_crit logical. If \code{TRUE} the alternative criterion is used for peeling.
#' @param use_NAs logical. If \code{TRUE} observations with missing values are included in the analysis.
#' @param seed seed to be set before the first iteration. Only useful for \code{B > 0}.
#' @param print_position logical. If \code{TRUE} the current position of the algorithm is printed out.
#'
#' @details The outcome of the \code{formula} can either be numeric, logical or a survival object (see \code{\link{Surv}}). If it is a survival object the \code{target} is set to the number of events per amount of time.
#'
#' @details The output of this function can become very large because all outputs of the singular peel function \code{PRIM_peel} are put together in one output.
#' Therefore it is usefull to remove all the dominated boxes (see \code{\link{remove_dominated}}).
#'
#' @return \code{PRIM_peel_bs} returns an object of class "\code{peel}", which is a list containing at least the following components:
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
#' @seealso \code{\link{remove_dominated}}, \code{\link{PRIM_peel}}, \code{\link{PRIM_paste}}, \code{\link{PRIM}}
#'
#' @references Friedman, J. H. and Fisher, N. I., 'Bump hunting in high-dimensional data', Statistics and Computing 9 (2) (1999), 123-143
#'
#' Ott, A. and Hapfelmeier, A., 'Nonparametric Subgroup Identification by PRIM and CART: A Simulation and Application Study', Computational and Mathematical Methods in Medicine, vol. 2017 (2017), 17 pages, Article ID 5271091
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
#' # apply the PRIM_peel_bs function:
#' prim <- PRIM_peel_bs(formula=y ~ ., data=dat, beta_min = .01)
#' plot(prim) # multiple trajectory
#' head(prim$box) # box definitions
#'
#' @export


PRIM_peel_bs <- function(formula, data, peel_alpha = seq(0.01, 0.4, 0.03), B = 0, beta_min= .01, target = mean, alter_crit=TRUE, use_NAs=TRUE, seed, print_position=TRUE){

  # generate the data.frame X the function is working with
  if(use_NAs) {if(!missing(formula)) X <- model.frame(formula, data=data, na.action = na.pass) else X <- data}
  if(use_NAs==FALSE) {if(!missing(formula)) X <- model.frame(formula, data=data, na.action = na.omit) else X <- na.omit(data)}

  # use the events per time as target, if the response is a survival object
  if(is.Surv(X[,1])) target <- function(k) sum(k[,2])/sum(k[,1])

  if(!missing(seed)) set.seed(seed)
  box <- NULL
  box_metric <- NULL
  box_nom <- NULL
  box_na <- NULL

  # generate B bootstrap samples for each alpha
  # and evaluate the PRIM function on each to get the box limits
  for(i in peel_alpha){
    if(print_position) cat("\nalpha =", i, "\n Iteration: " )
    for (j in 0:B){
      if(j==0) X_new <- X else X_new <- X[sample(1:nrow(X), replace = TRUE), ]
      prim <- PRIM_peel(data=X_new, peel_alpha=i, beta_min = beta_min, alter_crit = alter_crit, target = target)
      # sample again if all NAs of one variable are dropped out because of the bootstrap
      if(j!=0 & length(box_na[[1]])!=length(prim$box_na[[1]])) {
        repeat{
          if(j==0) X_new <- X else X_new <- X[sample(1:nrow(X), replace = TRUE), ]
          prim <- PRIM_peel(data=X_new, peel_alpha=i, beta_min = beta_min, alter_crit = alter_crit, target = target)
          if(length(box_na[[1]])==length(prim$box_na[[1]])) break
        }
      }
      # put all possible box definitions together
      box <- rbind(box, prim$box)
      if(!is.null(prim$box_metric)) box_metric <- rbind(box_metric, prim$box_metric)
      if(!is.null(prim$box_nom)) box_nom <- c(box_nom, prim$box_nom)
      if(!is.null(prim$box_na)) box_na <- c(box_na, prim$box_na)
      if(print_position) cat(j, " ")
    }
  }

  names(box) <- names(prim$box)

  # ordered factors have to be treaten like numeric variables
  for (i in 1:ncol(X)) if (is.ordered(X[,i])) X[,i] <- as.numeric(X[,i])

  # which variables have at least one missing value
  na_vars <- apply(apply(X, 2, is.na),2,any)
  if(na_vars[1]) {
    X <- X[!is.na(X[,1]),]
    na_vars[1] <- FALSE
  }

  # define subsets of the whole data that lie in the boxes
  sub <- lapply(1:nrow(box), function(k) {
    inbox(X, fixbox_metric=box_metric[k,], fixbox_nom=box_nom[[k]], fixbox_na = box_na[[k]])
  })

  # calculate target and support on the subsets
  f <- sapply(sub, function(k) target(X[k,1]))
  beta <- sapply(1:length(f), function(k) sum(sub[[k]])/nrow(X))

  # define the function output
  ret <- list(f=f, beta=beta, box=box, box_metric=box_metric, box_nom=box_nom, box_na=box_na, subsets=sub, data_orig=X, target=target)
  class(ret) <- "peel"
  return(ret)
}

