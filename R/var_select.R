#' Sequential Relevance of Box-defining Variables
#'
#' @description This function creates a table of the sequential relevance of the variables used in the definition of a box defined by \code{\link{define_fixbox}}.
#'
#' @param fixbox object of class "\code{fixbox}" to be used to create the table of the sequential variable relevances.
#'
#' @details This function does a sequential removal off all variables defining one box.
#' In each iteration every variable left in the box definition is tried to be completley left out of the definition.
#' From all variables this one is removed, which causes the least decrease of the target function evaluated on all observations lying in the box.
#'
#'
#' @return \code{var_select} returns a list with the following components:
#' \item{tab}{table showing the sequential relevances of the variables used in the definition of the \code{fixbox}.}
#' \item{fixboxes}{\code{list} of the boxes defined at every iteration step (by removing of one variable).}
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
#' # sequential variable relevances:
#' var_select(p$fixboxes[[1]])$tab
#'
#' @export


var_select <- function(fixbox){

  target <- fixbox$target
  X <- fixbox$data_orig
  fixboxes <- list(fixbox)
  n_met <- ifelse(!is.null(fixbox$fixbox_metric) ,ncol(fixbox$fixbox_metric)/2, 0)
  rem_vars <- NA
  vars <- c(substring(names(fixbox$fixbox_metric[1:n_met]), first = 5), names(fixbox$fixbox_nom))

  # loop: define the "new_boxes" that way, that every variable is tried to be completly left out of the box definition
  i <- 1
  while(TRUE){
    if(!is.null(fixboxes[[i]]$fixbox_metric)){
      new_boxes_metric <- lapply(1:n_met, function(k) {
        b <- fixboxes[[i]]$fixbox_metric
        b_na <- fixboxes[[i]]$fixbox_na
        b[,k] <- -Inf
        b[,n_met+k] <- Inf
        if(!is.null(b_na)) b_na[match(substring(names(b)[k], first = 5), names(b_na))] <- TRUE
        c(list(b), list(fixboxes[[i]]$fixbox_nom), list(b_na))
        })
    } else new_boxes_metric <- NULL

    if(!is.null(fixboxes[[i]]$fixbox_nom)){
      new_boxes_nom <- lapply(1:length(fixboxes[[i]]$fixbox_nom), function(k){
        b <- fixboxes[[i]]$fixbox_nom
        b_na <- fixboxes[[i]]$fixbox_na
        b[[k]] <- rep(TRUE, times=length(fixboxes[[i]]$fixbox_nom[[k]]))
        names(b[[k]]) <- names(fixboxes[[i]]$fixbox_nom[[k]])
        if(!is.null(b_na))b_na[match(names(b)[k], names(b_na))] <- TRUE
        c(list(fixboxes[[i]]$fixbox_metric), list(b), list(b_na))
      })
    } else new_boxes_nom <- NULL

    new_boxes <- c(new_boxes_metric, new_boxes_nom)

    # define the subsets of observations, that lie in each candidate box
    subsets <- lapply(new_boxes, function(k) inbox(X = X, fixbox_metric = k[[1]], fixbox_nom = k[[2]], fixbox_na = k[[3]]))

    # calculate the target functions on each subset and search the maximum target
    crit <- sapply(subsets, function(k) target(X[k,1]))
    crit[crit==fixboxes[[i]]$f] <- NA
    which_max <- which.max(crit)

    # end the loop if there is no maximum (i.e. there is no variable left in the box definition)
    if(length(which_max) == 0) break

    # define the new fixbox and save the name of the removed variable
    sub <- subsets[[which_max]]
    box <- new_boxes[[which_max]]
    fixboxes[[i+1]] <- list(f=crit[which_max], beta=sum(sub)/nrow(X), fixbox_metric=box[[1]], fixbox_nom=box[[2]], fixbox_na=box[[3]], subset=sub, data_orig=X, target=target)
    rem_vars[i+1] <- vars[which_max]

    i <- i+1
  }

  # create a data.frame showing the sequential relevances for the output
  f <- sapply(fixboxes, function(k) k$f)
  beta <- sapply(fixboxes, function(k) k$beta)
  a <- cbind.data.frame(rem_vars, f, beta)
  names(a) <- c("variable", "f", "beta")
  return(list(tab=a, fixboxes=fixboxes))
}


