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