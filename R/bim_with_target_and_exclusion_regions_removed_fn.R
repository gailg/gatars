#' @title Chop out those independent segments of bim that contain 
#' target markers or exclusion regions
#' 
#' @description For each chromosome, create independent segments of bim
#' (using hotspots), and throw out any segments that contain any target
#' snps or any part of an exclusion region.  Areas within hotspots are
#' also rejected.
#' 
#' @param params_sampling_set, and in particular the objects \code{bim},
#' \code{hotspot}, \code{target_markers}, and \code{exclusion_region}.
#' See \code{\link{params_sampling_set_fn}}.
#' 
#' @return A \code{data.frame} equal to \code{bim} with some lines removed.
#' The removed lines correspond to those snps that fall in those
#' independent segments containing \code{target_markers}
#' or any part of \code{exclusion_region} or any snps that fall inside
#' any of the hotspots described by \code{hotspot}
#' 
#' @examples 
#' #--------------- Get ready to call the function
#' bim = alternative_example$bim
#' epsilon_on_log_scale = 0.02
#' table(bim$chromosome)
#' bim_for_chromosome_21 = bim[bim$chromosome == 21, ]
#' head(bim_for_chromosome_21, 20)
#' target_markers = c("exm-rs1006899",  "exm1564039")
#' bim[bim$snp %in% target_markers, ]
#' exclusion_region = NULL
#' genotype = alternative_example$genotype
#' params_sampling_set = params_sampling_set_fn(
#'   bim, epsilon_on_log_scale, exclusion_region,
#'   genotype, hotspot, target_markers)
#' names(params_sampling_set)
#' #--------------- Call the function
#' bim_with_target_and_exclusion_regions_removed = 
#'   bim_with_target_and_exclusion_regions_removed_fn(params_sampling_set)
#' #--------------- Notice that some lines have been removed
#' str(bim_with_target_and_exclusion_regions_removed)
#' str(bim)
#' #--------------- Focus on chromosome 21 because it contains the fewest snps
#' before = bim_for_chromosome_21 = bim[bim$chromosome == 21, ]
#' str(before)
#' after = bim_with_target_and_exclusion_regions_removed_for_chromosome_21 = 
#'   bim_with_target_and_exclusion_regions_removed[
#'     bim_with_target_and_exclusion_regions_removed$chromosome == 21, ]
#' str(after)
#' #--------------- These are the bp that have been removed
#' setdiff(before$bp, after$bp)
#' #--------------- Does this make sense?
#' hhh = hotspot_for_chromosome_21 = hotspot[hotspot$chromosome == 21, ]
#' bim[bim$snp %in% target_markers, ]
#' hhh[15600000 < hhh$center & hhh$center < 16000000, ]
#' bim_for_target_markers = bim[bim$snp %in% target_markers, ]
#' to_the_left = hhh[hhh$end < bim_for_target_markers[1, ]$bp, ]
#' to_the_left[nrow(to_the_left), ]
#' to_the_right = hhh[bim_for_target_markers[1, ]$bp < hhh$start, ]
#' to_the_right[1, ]
#' bim[bim$snp %in% target_markers, ]
#' bim_remove_due_to_first_target_marker = bim_for_chromosome_21[
#'   to_the_left[nrow(to_the_left), ]$end <= bim_for_chromosome_21$bp &
#'   bim_for_chromosome_21$bp <= to_the_right[1, ]$start, ]
#' #--------------- The part of bim that are removed due to the 
#----------------first target marker
#' bim_remove_due_to_first_target_marker  
#' #--------------- Repeat the above for the second_target_marker.
#' #--------------- This shows that sometimes not just the target
#' #--------------- marker is removed.
#' to_the_left = hhh[hhh$end < bim_for_target_markers[2, ]$bp, ]
#' to_the_left[nrow(to_the_left), ]
#' to_the_right = hhh[bim_for_target_markers[2, ]$bp < hhh$start, ]
#' to_the_right[1, ]
#' bim_remove_due_to_second_target_marker = bim_for_chromosome_21[
#'   to_the_left[nrow(to_the_left), ]$end <= bim_for_chromosome_21$bp &
#'   bim_for_chromosome_21$bp <= to_the_right[1, ]$start, ]
#' bim_remove_due_to_second_target_marker  
#' #--------------- The function also removes the parts of bim
#' #--------------- that fall inside the hotspots
#' bim_removed_because_they_are_in_hotspots = bim_for_chromosome_21[
#'   sapply(bim_for_chromosome_21$bp, function(this){
#'     any(hhh$start <= this & this <= hhh$end)
#'   }), ]
#' bim_removed_because_they_are_in_hotspots$bp
#' #--------------- Check to see if the parts removed because of the 
#' #--------------- target markers and those removed because they fall
#' #--------------- within hotspots match those removed in the 
#' #--------------- original function call
#' all.equal(sort(c(bim_remove_due_to_first_target_marker$bp,
#'                  bim_remove_due_to_second_target_marker$bp,
#'                  bim_removed_because_they_are_in_hotspots$bp)),
#'           setdiff(before$bp, after$bp))
#'           
#' @export
bim_with_target_and_exclusion_regions_removed_fn = function(params_sampling_set){
  bim = params_sampling_set$bim
  hotspot = params_sampling_set$hotspot
  target_markers = params_sampling_set$target_markers
  target_bim = bim[bim$snp %in% target_markers, ]
  target_part = with(target_bim, data.frame(chromosome, begin = bp, end = bp))
  exclusion_region = rbind(params_sampling_set$exclusion_region, target_part)
  exclusion_chromosome = unique(exclusion_region$chromosome)
  chromosomes_with_target_and_exlusion_regions_removed = lapply(1:22, function(chromosome){  
    answer = if(chromosome %in% exclusion_chromosome){
      segments_0 = independent_segment_fn(bim, chromosome, hotspot)
      eee = exclusion_region[exclusion_region$chromosome == chromosome, ]
      intersect_with_genotype_q = sapply(segments_0, function(this){
        any(unlist(lapply(1:nrow(eee), function(kkk){
          any(eee[kkk, ]$begin : eee[kkk, ]$end %in% this$bp)
        })))
      })
      do.call(rbind, segments_0[!intersect_with_genotype_q])
    } else {
      bim[bim$chromosome == chromosome, ]
    }
    answer
  })
  do.call(rbind, chromosomes_with_target_and_exlusion_regions_removed)
}

# "2016-10-19 12:35:20 PDT" handles exclusion_region