#' Genetic Association Tests for Arbitrarily Related Subjects
#'
#' Test the association between a specified set of \code{MMM} target markers and a binary
#' or quantitative trait using subjects with any genealogical relationship.
#' 
#' \strong{The three functions featured in \url{http://stanford.edu/~ggong/gatars/} are:}
#' 
#' \code{gatars_sampling_set}: Create the sampling sets that are used in genome resampling.
#' 
#' \code{gatars_test_size}: Uses \code{davies} to calculate the p-values of the basic statistics,
#'  and genome resampling to calculate the p-values of the optimized statistics.       
#'   
#' \code{matrix_image_fn}: A little helper function to display the Psi matrix.
#' 
#' \strong{The datasets included are:}
#' 
#' \code{gatars_example}: A list of the 6 objects you need to provide to run \code{gatars}. Use this example to learn how to run \code{gatars}.
#' 
#' \code{hotspot}: The hotspot data used to create the independent segments that are used by \code{gatars_sampling_set}.
#' 
#' 
#' @docType package
#' @name gatars
NULL
