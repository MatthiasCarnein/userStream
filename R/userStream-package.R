#' userStream
#'
#' Stream Clustering algorithm applicable to custoemr segmentation. The algorithm allows to insert new but also update already processed observations.
#'
#' @author Matthias Carnein \email{Matthias.Carnein@@uni-muenster.de}
#'
#' @name userStream
#' @docType package
#' @useDynLib userStream
#' @import Rcpp
NULL

loadModule("MOD_userStream", TRUE)
