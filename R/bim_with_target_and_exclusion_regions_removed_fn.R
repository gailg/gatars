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
#' 
#' #--------------- bim and target_markers
#' bim = gatars_example$bim
#' target_markers = c("exm1055449", "exm1562514")
#' target_markers
#' bim[bim$snp %in% target_markers, ]
#' #--------------- independent_segment_containing_first_target_marker
#' independent_segments_13 = independent_segment_fn(bim, 13, hotspot)
#' is_this_the_one = sapply(independent_segments_13, function(this){
#'   any(target_markers %in% this$snp)
#' })
#' which(is_this_the_one)
#' independent_segment_containing_first_target_marker = 
#'   independent_segments_13[is_this_the_one][[1]]
#' #--------------- independent_segment_containing_second_target_marker
#' independent_segments_21 = independent_segment_fn(bim, 21, hotspot)
#' is_this_the_one = sapply(independent_segments_21, function(this){
#'   any(target_markers %in% this$snp)
#' })
#' which(is_this_the_one)
#' independent_segment_containing_second_target_marker = 
#'   independent_segments_21[is_this_the_one][[1]]
#' #--------------- inside_a_hotspot
#' inside_a_hotspot = bim[bim$snp == "exm1562489", ]
#' the_hotspot_line_containing_iah = 
#'   hotspot[hotspot$chromosome == inside_a_hotspot$chromosome &
#'             hotspot$start < inside_a_hotspot$bp & inside_a_hotspot$bp < hotspot$end, ]
#' list(inside_a_hotspot = inside_a_hotspot,
#'      the_hotspot_line_containing_iah = the_hotspot_line_containing_iah)
#' #--------------- throw_in_some
#' throw_in_some = do.call(rbind, list(
#'   independent_segments_13[[9]],
#'   independent_segments_21[[10]]))
#' #--------------- add pretty row.names to make things clearer
#' row.names(independent_segment_containing_first_target_marker) = 
#'   paste0("t1_", 1:nrow(independent_segment_containing_first_target_marker))
#' row.names(independent_segment_containing_second_target_marker) = 
#'   paste0("t2_", 1:nrow(independent_segment_containing_second_target_marker))
#' row.names(inside_a_hotspot) = 
#'   paste0("h_", 1:nrow(inside_a_hotspot))
#' row.names(throw_in_some) = 
#'   paste0("throw_", 1:nrow(throw_in_some))
#' #--------------- bim_baby
#' bim_baby = do.call(rbind, list(
#'   independent_segment_containing_first_target_marker,
#'   independent_segment_containing_second_target_marker,
#'   inside_a_hotspot,
#'   throw_in_some))
#' #--------------- show and tell
#' list(independent_segment_containing_first_target_marker =
#'        independent_segment_containing_first_target_marker,
#'      independent_segment_containing_second_target_marker =
#'        independent_segment_containing_second_target_marker,
#'      inside_a_hotspot = inside_a_hotspot,
#'      bim_baby = bim_baby)
#' #--------------- a fake params_sampling_set
#' params_sampling_set = list(
#'   bim = bim_baby,
#'   exclusion_region = NULL,
#'   hotspot = hotspot,
#'   target_markers = target_markers)
#' #--------------- the answer from bim_with... should be the same as throw_some_in
#' bim_with_target_and_exclusion_regions_removed_fn(params_sampling_set)
#' 
#' 
#' @export
bim_with_target_and_exclusion_regions_removed_fn = function(params_sampling_set){
  bim = params_sampling_set$bim
  hotspot = params_sampling_set$hotspot
  target_markers = params_sampling_set$target_markers
  target_bim = bim[bim$snp %in% target_markers, ]
  target_part = with(target_bim, data.frame(chromosome, start = bp, end = bp))
  exclusion_region = rbind(params_sampling_set$exclusion_region, target_part)
  exclusion_chromosome = unique(exclusion_region$chromosome)
  chromosomes_with_target_and_exlusion_regions_removed = lapply(1:22, function(chromosome){  
    answer = if(chromosome %in% exclusion_chromosome){
      segments_0 = independent_segment_fn(bim, chromosome, hotspot)
      if( is.null(segments_0) ) {
        NULL
      } else {
        eee = exclusion_region[exclusion_region$chromosome == chromosome, ]
        intersect_with_genotype_q = sapply(segments_0, function(this){
          any(unlist(lapply(1:nrow(eee), function(kkk){
            any(eee[kkk, ]$start : eee[kkk, ]$end %in% this$bp)
          })))
        })
        do.call(rbind, segments_0[!intersect_with_genotype_q])
      }
    } else {
      bim[bim$chromosome == chromosome, ]
    }
    answer
  })
  do.call(rbind, chromosomes_with_target_and_exlusion_regions_removed)
}

# "2016-10-19 12:35:20 PDT" handles exclusion_region