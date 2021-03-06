% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define_fixbox.R
\name{define_fixbox}
\alias{define_fixbox}
\title{Defining a "fixbox"-Object}
\usage{
define_fixbox(prim, step)
}
\arguments{
\item{prim}{object of class "\code{peel}" from which the box should be defined.}

\item{step}{number of the step (i.e. row in the argument \code{box} of the "\code{peel}"-object.) of the requested box.}
}
\value{
\code{define_fixbox} returns a object of class "\code{fixbox}", which is the basically same as a "\code{peel}"-object, but only contains one box.
It is a list containing at least the following elements:

\item{f}{target function of the selected box.}

\item{beta}{support of the selected box.}

\item{box}{a \code{\link{data.frame}} with one row defining the borders of the box. The columns with "\code{min.}" and "\code{max.}" describe the lower and upper boundaries of the at least ordinal covariates. Therefore the value taken is the last one that is \bold{not} included in the box.

For the nominal variables there are columns for every category they can take. If the category is removed from the box the value \code{FALSE} is taken. The names of these columns are structured like: \code{<variable name>.<category>}

For each variable with missing values (only if \code{use_NAs = TRUE}) there is also a column taking the value \code{FALSE} if the \code{NA}s of this variable are removed from the box. The names of these columns are structured like: \code{<variable name>.NA}
}

\item{box_metric, box_nom, box_na}{easier to handle definitions of the box for other functions}

\item{subset}{logical vector indicating the subset (i.e. the observations that lie in the box).}

\item{data_orig}{original dataset that was used for the peeling.}
}
\description{
Function for determining on one box out of a "\code{peel}"-object.
}

