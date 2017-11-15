#' @title Create sampling sets for genome resampling
#' 
#' @description Testing the null hypothesis for the optimized statistics requires
#' sampling sets which are sets of markers located throughout the autosomal
#' genome satisfying two requirements. \code{gatars_sampling_set} 
#' creates \code{MMM} such sampling sets,  one for 
#' for each target marker, such that (1) the sampling set for the \code{mmm}-th target 
#' marker contains markers all of whose minor allele frequencies match that 
#' of the \code{mmm}-th target marker and (2) markers in the sampling 
#' sets are independent of the target markers and of any markers known to be 
#' associated with the trait.
#' 
#' @details Each column in \code{genotype} corresponds to a row in \code{bim}, and both 
#' correspond to a marker. Assuming that markers within two consecutive hotspots 
#' are independent of those within any other two consecutive hotspots, the recombination 
#' hotspots divide the autosomal genome contained in the columns of \code{genotype}
#' into independent segments.  \code{gatars_sampling_set} removes
#'  all segments containing any target markers or 
#' any markers defined by the \code{exclusion_region} data set. On the markers in the 
#' remaining segments, \code{gatars_sampling_set} calculates the empirical 
#' minor allele frequencies and says that 
#' a marker's minor allele frequency  matches the minor allele frequency \code{pi[mmm]} 
#' of the \code{mmm}-th target marker if it falls in the closed interval
#' \code{pi[mmm] * [1 - epsilon, 1 + epsilon]}.  (If the number of markers satisfying the 
#' matching requirement exceeds \code{1000}, \code{gatars_sampling_set} 
#' randomly chooses \code{1000}.)
#' 
#' @inheritParams params_sampling_set_fn
#' 
#' @return A list containing the following objects
#' \itemize{
#' \item{\code{params_sampling_set}: } {
#' A list containing objects that will be useful in the calculations of 
#' \code{gatars_test_size}.
#' The result of a call to \code{params_sampling_set}.  
#' See \code{\link{params_sampling_set_fn}}.
#' }
#' \item{\code{report}: } {
#' A \code{data.frame} containing \code{MMM} rows and the columns
#' \code{min}, \code{pi}, \code{max}, and \code{set_size}.
#' \code{min}/\code{max} contains the smallest/largest 
#' minor allele frequency among
#' the columns in the \code{mmm}-th sampling set, \code{pi}
#' the minor allele frequency of the \code{mmm}-th target marker, and \code{set_size}
#' the number of columns in the \code{mmm}-th sampling set.
#' }
#' \item{\code{sampling_set}: }{
#' A list of \code{MMM} matrices, one matrix for each target marker.
#' The \code{mmm}-th matrix is the sampling set for the \code{mmm}-th
#' target marker and has \code{NNN} rows and up to \code{1000} columns, each 
#' column containing a column from \code{genotype}.  These columns
#' do not intersect with any of the target markers or exclusion regions
#' and the minor allele frequencies of the columns in \code{mmm}-th sampling set match
#' the minor allele frequency of the \code{mmm}-th target marker.
#' }
#' }
#' 
#' @template params_sampling_set_examples
#' @examples
#' sampling_set = gatars_sampling_set(
#'     bim, epsilon, exclusion_region,
#'     genotype, hotspot, target_markers)
#' print(sampling_set)
#' names(sampling_set)
#' names(sampling_set$params_sampling_set)
#' sampling_set$report
#' str(sampling_set$sampling_set)
#' @export
gatars_sampling_set = function(
  bim,
  epsilon,
  exclusion_region,
  genotype,
  hotspot,
  target_markers
){
  params_sampling_set = params_sampling_set_fn(
    bim,
    epsilon,
    exclusion_region,
    genotype,
    hotspot,
    target_markers
    )
  sss = sampling_set_fn(params_sampling_set)
  answer = list(params_sampling_set = params_sampling_set,
                report = sss$report,
                sampling_set = sss$sampling_set)
  class(answer) = c("gatars_sampling_set", class(answer))
  answer
}
