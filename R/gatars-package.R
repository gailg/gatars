#' Genetic Association Tests for Arbitrarily Related Subjects
#'
#' Calculate the p-values of a family of statistics that test the association between a specifed set of genetic markers and a binary or quantitative trait.
#' 
#' \strong{The three functions used in \url{http://stanford.edu/~ggong/gatars/} are:}
#' 
#' \code{gatars_sampling_set}: Create the sampling sets that are used in genome resampling.
#' 
#' \code{gatars}: Uses \code{davies} to calculate the p-values of the basic statistics,
#'  and genome resampling to calculate the p-values of the optimal statistics.       
#'   
#' \code{matrix_image_fn}: A little helper function to display the Psi matrix.
#' 
#' \strong{The datasets used are:}
#' 
#' \code{gatars_example}: A list of the 6 objects you need to provide to run \code{gatars}.
#' 
#' \code{hotspot}: The hotspot data used to create the independent segments that are used by \code{gatars_sampling_set}.
#' 
#' The \code{hotspot} data are from “A Fine-Scale Map of Recombination Rates and Hotspots Across the Human Genome”, Myers et.al., \strong{Science} 
#' October 2005.
#' 
#' @docType package
#' @name gatars
NULL
