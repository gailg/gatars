#' @title Create sampling sets
#' 
#' @description For each target marker, gatars requires a sampling set,
#' a collection of columns from \code{genotype}.  Each sampling set
#' must contain snps that (1) do not intersect with any of the target snps or 
#' exclusion regions, and (2) have mafs (minor allele frequencies)
#' that match the maf of its
#' corresponding target snp.
#' 
#' @param params_sampling_set, and in particular the objects 
#' \code{epsilon}, \code{MMM}, and \code{p_target}.
#' See \code{\link{params_sampling_set_fn}}.
#' 
#' @return A \code{list} containing the following two objects
#' \itemize{
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
#' @examples 
#' bim = gatars_example$bim
#' epsilon = 0.01
#' exclusion_region = gatars_example$exclusion_region
#' genotype = gatars_example$genotype
#' target_markers = gatars_example$target_markers[3:5]
#' set.seed(1)
#' params_sampling_set = params_sampling_set_fn(
#'   bim, epsilon, exclusion_region,
#'   genotype, hotspot, target_markers)
#' names(params_sampling_set)
#' sampling_set = sampling_set_fn(params_sampling_set)
#' names(sampling_set)
#' sampling_set$report
#' 
#' @export
sampling_set_fn = function(params_sampling_set){
  genotype = params_sampling_set$genotype
  MMM = params_sampling_set$MMM
  epsilon = params_sampling_set$epsilon
  # p_target are the mafs of the target snps and the mafs we want to match
  p_target = params_sampling_set$p_target
  # Remove from bim and genotype the regions specified by exclusion_region
  # The working files are then bim_with_target_and_exclusion_region_removed
  # and genotype_with_stuff_removed.
  # p_sim is the maf genotype_with_stuff_removed
  bim_with_target_and_exclusion_regions_removed = 
    bim_with_target_and_exclusion_regions_removed_fn(params_sampling_set)
  genotype_with_stuff_removed = genotype[, bim_with_target_and_exclusion_regions_removed$snp]
  p_sim = unname(colMeans(genotype_with_stuff_removed)/2)
  # For sampling_set I march through each mmm = 1, ...M, 
  # for each mmm, I find the good_ones, the columns of genotype_with_stuff_removed
  # whose mafs match.  
  # I use epsilon to parametrize the matching.
  # Notice that if I get more than 1000 good ones, I sample to get at most 1000.
  sampling_set = lapply(1:MMM, function(mmm){
    p_this = p_target[mmm]
    lower = p_this * (1 - epsilon)
    upper = p_this * (1 + epsilon)
    good_ones_q = (lower < p_sim & p_sim < upper)
    good_ones = which(good_ones_q)
    N_good_ones = length(good_ones)
    if(N_good_ones > 1000){
      good_ones = sort(sample(good_ones, 1000))
    }
    answer =  genotype_with_stuff_removed[, good_ones, drop = FALSE]
    answer
  })
  # report tells me how many good ones I got in each sampling set.
  report = do.call(rbind, lapply(1:MMM, function(mmm){
    this_genotype = sampling_set[[mmm]]
    p_sim = unname(colMeans(this_genotype)/2)
    data.frame(min = round(min(p_sim), 5),
               p_target = round(p_target[mmm], 5),
               max = round(max(p_sim), 5),
               set_size = ncol(this_genotype))
  }))
  rownames(report) = NULL
  report
  answer = list(
    report = report,
    sampling_set = sampling_set)
  answer
}
