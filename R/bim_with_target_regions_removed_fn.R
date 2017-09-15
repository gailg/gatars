#' @export
bim_with_target_regions_removed_fn = function(params_sampling_set){
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