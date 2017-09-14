#' @title Create sampling sets for genome resampling
#' 
#' @description Each column in \code{genotype} corresponds to a row in \code{bim}, and both 
#' correspond to a marker or snp. After removing the target markers, we want to form, 
#' from the remaining snps, \code{MMM} sampling sets, one for each of the target markers.  
#' The two conditions we require of the sampling sets are \strong{the matching requirement:}
#' the snps in the sampling set for a target marker has minor allele frequencies that 
#' match closely with that of the target marker, and \strong{the independence requirement:}
#' any snp from any sampling set is statistically independent of the target markers and any 
#' marker in the \code{exclusion_region}.
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
#' \item{\code{bim}: } {
#' An echo of the input \code{bim}.
#' }
#' \item{\code{genotype}: } {
#' An echo of the input \code{bim}.
#' }
#' \item{\code{e_g_target}: }{
#' A matrix with `NNN` rows and `MMM` columns.  All rows are the same.  
#' The \code{mmm}-th column contains the mean of the \code{mmm}-th
#' column of \code{g_target}.
#' }
#' }
#' 
#' @template gatars_sampling_set_examples
#' @export
gatars_sampling_set = function(
  bim,
  genotype,
  target_markers,
  exclusion_region,
  hotspot,
  epsilon_on_log_scale = 0.02
){
  params_sampling_set = params_sampling_set_fn(
    bim,
    epsilon_on_log_scale,
    exclusion_region,
    genotype,
    hotspot,
    target_markers
    )
  answer = sampling_set_fn(params_sampling_set)
  class(answer) = c("gatars_sampling_set", class(answer))
  answer
}
