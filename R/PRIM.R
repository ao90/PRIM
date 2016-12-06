#' Combined Function for the Patient Rule Induction Method (PRIM)
#'
#' @import survival
#' @importFrom graphics barplot hist par plot
#' @importFrom stats model.frame na.omit na.pass quantile setNames
#'
#' @description This function is a automated implementation of PRIM as suggested by Friedman and Fisher (1999).
#' It includes multiple peeling (\code{\link{PRIM_peel_bs}}), pasting (\code{\link{PRIM_paste}}) and the covering stretegy to find more than one box.
#'
#' @param formula an object of class "\code{\link{formula}}" with a response but no interaction terms.
#' It indicates the response over which the target function should be maximized and the covariates that are used for the later box definitions.
#' @param data an object of class \code{\link{data.frame}} containing the variables named in the formula.
#' @param f_min minimum target the final box must have. From all boxes, that fulfill this criterion, the one with the biggest support is taken after the peeling. If this argument is missing the box with the biggest target having at least a support of \code{beta_min} is taken.
#' @param beta_min minimum support that one Box must have.
#' @param max_boxes maximum number of boxes to be found.
#' @param peel_alpha vector of a sequence of different alpha-fractions used for the peelings.
#' @param B number of bootstrap samples on which the peeling is applied to for each alpha. For \code{B = 0} no bootstraps are created.
#' @param target target-function to be maximized.
#' @param alter_crit logical. If \code{TRUE} the alternative criterion is used for peeling.
#' @param use_NAs logical. If \code{TRUE} observations with missing values are included in the analysis.
#' @param seed seed to be set before the first iteration. Only useful for \code{B > 0}.
#' @param print_position logical. If \code{TRUE} the current position of the algorithm is printed out.
#' @param paste_alpha alpha-fraction that is pasted to the box at each pasting step
#' @param max_steps maximum number of pasting steps the function should make.
#' @param stop_by_dec logical. If \code{TRUE} the pasting stops if the target at one step is lower than the target of the last step.
#'
#' @details This function repeats the peeling and pasting algorithm for the same settings of the metaparameters until a stop ctiterion is reached.
#' After each iteration the observations already included in a box are removed from the data, on which the next box is built. This strategy is called covering.
#' This iteration stops if either \code{max_boxes} is reached or if the target function of the "best" box is lower than the overall target.
#'
#' In each iteration step this function does a multiple peeling characterized by the sequenz \code{alpha_peel} and \code{B}.
#' From the peeling output the box defined by \code{beta_min} and \code{f_min} is chosen.
#' After that the pasting function seeks for boxes with bigger supports and bigger targets and takes the one with the highest target function within the box.
#'
#' @return \code{PRIM} returns an object of class "\code{prim}", which is a list containing at least the following components:
#' \item{f}{vector of the target functions evaluated on each box. The last element is the target of all observations not lying in a box.}
#' \item{beta}{vector of the supports of each box. The last element is the fraction of observations not lying in a box.}
#' \item{box}{a \code{\link{data.frame}} defining the borders of the boxes. Each row belongs to one box. The columns with "\code{min.}" and "\code{max.}" describe the lower and upper boundaries of the at least ordinal covariates. Therefore the value taken is the last one that is \bold{not} included in the current box.
#'
#' For the nominal variables there are columns for every category they can take. If the category is removed from the box the value \code{FALSE} is taken. The names of these columns are structured like: \code{<variable name>.<category>}
#'
#' For each variable with missing values (only if \code{use_NAs = TRUE}) there is also a column taking the value \code{FALSE} if the \code{NA}s of this variable are removed from the current box. The names of these columns are structured like: \code{<variable name>.NA}
#' }
#' \item{box_metric, box_nom, box_na}{easier to handle definitions of the boxes for other functions}
#' \item{subsets}{\code{list} of logical vectors indicating the subsets (i.e. the observations that lie in each box)}
#' \item{fixboxes}{\code{list} of all \code{fixbox}'es defining the final boxes.}
#' \item{data_orig}{original dataset that is used.}
#'
#' @references Friedman, J. H. and Fisher, N. I., 'Bump hunting in high-dimensional data', Statistics and Computing \bold{9}, 123-143
#'
#' @seealso \code{\link{PRIM_peel_bs}}, \code{\link{PRIM_paste}}, \code{\link{define_fixbox}}
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
#' # apply the PRIM function to find the best boxes with a support of at least 0.1:
#' p <- PRIM(y~., data=dat, beta_min = 0.1, max_boxes = 3)
#' p
#'
#'
#'
#' @export



PRIM <- function(formula, data, f_min, beta_min = 0.2, max_boxes = Inf, peel_alpha = seq(0.01, 0.4, 0.03), B = 0, target = mean, alter_crit = TRUE, use_NAs = TRUE, seed, print_position = TRUE, paste_alpha = 0.01, max_steps = 50, stop_by_dec = TRUE){

  if(use_NAs==FALSE) {data <- na.omit(data)}
  row.names(data) <- 1:nrow(data)
  # use the events per time as target, if the response is a survival object
  if(is.Surv(model.frame(formula, data=data)[,1])) target <- function(k) sum(k[,2])/sum(k[,1])

  if(!missing(seed)) set.seed(seed)

  fixboxes <- list(NULL)
  i <- 1
  remove_obs <- list(NULL)

  while(TRUE){
    if(print_position) cat("\n\nBox", i)

    # run the peeling on the current data (first step = whole data set; any other step = data without observations included in an earlier box)
    if(i==1) {
      peel <- PRIM_peel_bs(formula = formula, data = data, B = B, peel_alpha = peel_alpha, beta_min = beta_min, alter_crit = alter_crit, print_position = print_position, use_NAs=use_NAs)
    } else {
      peel <- PRIM_peel_bs(formula = formula, data = data[-(remove_obs[[i]]), ], B = B, peel_alpha = peel_alpha, beta_min = beta_min, alter_crit = alter_crit, print_position = print_position, use_NAs=use_NAs)
    }

    # set all values for f and beta to NA, which belong to boxes with supports smaller than beta_min (these boxes must not be used)
    if(i==1) {
      peel$f[peel$beta < beta_min] <- NA
      peel$beta[peel$beta < beta_min] <- NA
    } else {
      peel$f[peel$beta < (beta_min*nrow(data)/nrow(data[-remove_obs[[i]], ]))] <- NA
      peel$beta[peel$beta < (beta_min*nrow(data)/nrow(data[-remove_obs[[i]], ]))] <- NA
    }

    # search for the highest beta out of all boxes with f > f_min, if f_min is given
    # search for the highest f out of all boxes, if the argument f_min is missing
    if(!missing(f_min)) {
      peel$beta[peel$f < f_min] <- NA
      which_max_peel <- which.max(peel$beta)
    } else which_max_peel <- which.max(peel$f)
    # if no box fulfills the constraints take the first one so that no error is produced
    if(length(which_max_peel) == 0) which_max_peel <- 1

    # define the fixbox (result of the peeling)
    fixbox_peel <- define_fixbox(peel, which_max_peel)

    # run the pasting on the fixbox and select the box with the highest f from the output
    Paste <- PRIM_paste(fixbox = fixbox_peel, paste_alpha = paste_alpha, max_steps = max_steps, stop_by_dec = stop_by_dec)
    which_max_paste <- which.max(Paste$f)
    # define the new fixbox
    fixboxes[[i]] <- define_fixbox(Paste, which_max_paste)

    # break the loop if the target of the found box is lower than the overall target or if max_boxes is reached
    # "last" is an indicator, whether the last box found should be included in the output or not
    if(fixboxes[[i]]$f <= target(fixboxes[[1]]$data_orig[,1]) ){
      last <- 1
      break
    }
    if(i == max_boxes ) {
      last <- 0
      break
    }

    # which observations should be removed after each step (covering)
    remove_obs[[i+1]] <- c(remove_obs[[i]], as.numeric(row.names(fixboxes[[i]]$data_orig))[fixboxes[[i]]$subset])

    i <- i+1
  }

  # extract the values for "box", "box_metric", "box_nom", "box_na" and "subsets" out of the list of fixboxes
  X <- fixboxes[[1]]$data_orig
  f  <- box <- subsets <- box_metric <- box_nom <- box_na <- NULL
  if(length(fixboxes) > last){
    for(i in 1:(length(fixboxes)-last)){
      f <- c(f,fixboxes[[i]]$f)
      box <- rbind(box,fixboxes[[i]]$fixbox)
      if(!is.null(fixboxes[[1]]$fixbox_metric)) box_metric <- rbind(box_metric,fixboxes[[i]]$fixbox_metric)
      if(!is.null(fixboxes[[1]]$fixbox_nom)) box_nom <- c(box_nom, list(fixboxes[[i]]$fixbox_nom))
      if(!is.null(fixboxes[[1]]$fixbox_na)) box_na <- c(box_na, list(fixboxes[[i]]$fixbox_na))
      if(i==1){subsets[[i]] <- fixboxes[[i]]$subset} else{
        subsets[[i]] <- rep(FALSE, nrow(X))
        subsets[[i]][-remove_obs[[i]]] <- fixboxes[[i]]$subset
      }

    }
    # if no box is found under these constraints, give back a warning
  } else {
    warning("no box found")
    ret <- list(subsets=list(rep(FALSE, times=nrow(X))), data_orig=X, target=target)
    class(ret) <- "PRIM"
    return(ret)
  }

  # define the last element of "subsets" as all observations, that lie in no box
  v <- rep(FALSE, nrow(X))
  for(i in 1:length(subsets)){
    v <- v | subsets[[i]]
  }
  subsets[[length(subsets)+1]] <- !v

  # calculate the target and beta on the subsets
  f[[length(f)+1]] <- target(X[subsets[[length(subsets)]],1])
  beta <- sapply(1:length(subsets), function(k) sum(subsets[[k]])/nrow(X))

  row.names(box) <- 1:nrow(box)
  if(!is.null(box_metric)) row.names(box_metric) <- 1:nrow(box_metric)

  # define the function output
  ret <- list(f=f, beta=beta, box=box, box_metric=box_metric, box_nom=box_nom, box_na=box_na, subsets=subsets, fixboxes=fixboxes, data_orig=X, target=target, remove_obs=remove_obs)
  class(ret) <- "PRIM"
  return(ret)
}

