#' @title Create sampling sets for genome resampling
#' 
#' @description Each column in \code{genotype} corresponds to a row in \code{bim}, and both 
#' correspond to a marker or snp. After removing the target markers and exclusion regions,
#' form \code{MMM} sampling sets, one for each of the target markers.  
#' The two conditions we require of the sampling sets are \strong{the matching requirement:}
#' the snps in the sampling set for a target marker has minor allele frequencies that 
#' match closely with that of the target marker, and \strong{the independence requirement:}
#' any snp from any sampling set is statistically independent of the target markers and any 
#' marker in the exclusion regions.
#' 
#' @details We assume that hotspots from Myers et.al. cut the genome into independent segments, 
#' and so snps residing within two consecutive hotspots are independent of snps residing 
#' within another two consecutive hotspots. Defining a segment to be a set of snps that 
#' all lie within two consecutive hotspots, we obtain a (large) set of independent segments. 
#' To satisfy \strong{the independence requirement} we need only remove any segments that
#' contain target markers or any defined by the \code{exclusion region}.
#' 
#' @inheritParams params_sampling_set_fn
#' 
#' @return A list containing the following 11 objects
#' \itemize{
#' \item{\code{params_sampling_set}: } {
#' The values of a call to \code{params_sampling_set}.  
#' See \code{\link{params_sampling_set_fn}}.
#' }
#' \item{\code{report}: } {
#' A \code{data.frame} containing \code{MMM} rows and the columns
#' \code{min}, \code{p_target}, \code{max}, and \code{set_size}.
#' \code{min}/\code{max} contains the smallest/largest maf among
#' the columns in the \code{mmm}-th sampling set, \code{p_target}
#' the maf of the \code{mmm}-th target marker, and \code{set_size}
#' the number of columns in the \code{mmm}-th sampling set.
#' }
#' \item{\code{sampling_set}: }{
#' A list of \code{MMM} matrices, one matrix for each target snp.
#' The \code{mmm}-th matrix is the sampling set for the \code{mmm}-th
#' target snp and has \code{NNN} rows and up to 1000 columns, each 
#' column containing a column from \code{genotype}.  These columns
#' do not intersect with any of the target snps or exclusion regions
#' and the mafs of the columns in \code{mmm}-th sampling set match
#' the maf of the \code{mmm}-th target snp.
#' }
#' }
#' 
#' @template params_sampling_set_examples
#' @examples
#' sampling_set = gatars_sampling_set(
#'     bim, epsilon_on_log_scale, exclusion_region,
#'     genotype, hotspot, target_markers)
#' print(sampling_set)
#' str(sampling_set$sampling_set)
#' @export
gatars_sampling_set = function(
  bim,
  epsilon_on_log_scale,
  exclusion_region,
  genotype,
  hotspot,
  target_markers
){
  params_sampling_set = params_sampling_set_fn(
    bim,
    epsilon_on_log_scale,
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
