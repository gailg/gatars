#' @export
independent_segment_0_fn = function(chromosome, bim, hotspot){
  bim_this_chromosome = bim[bim$chromosome == chromosome, ]
  hot = hotspot_this_chromosome = hotspot[hotspot$chromosome == chromosome, ]
  # bim_just_these removes snps that are inside hotspots
  bim_just_these = bim_this_chromosome[!sapply(bim_this_chromosome$bp, function(this){
    any(hot$start <= this & this <= hot$end)
  }), ]
  bp = bim_just_these$bp
  # hot_just_these removes head and tail rows that do no cutting
  hot_just_these = hot[bp[1] <= hot$center & hot$center <= bp[length(bp)], ]
  center = hot_just_these$center
  # bp_index is the index of the last snp left of this hot$center
  bp_index = sapply(center, function(this){
    max(which(bp < this))
  })
  # breaks will be used to separate the bp index numbers
  # segments is the answer
  breaks = unique(bp_index)
  numbers = 1:length(bp)
  breaks_with_ends = unique(c(0, breaks, max(numbers)))
  left = breaks_with_ends[-length(breaks_with_ends)]
  right = breaks_with_ends[-1]
  segment_numbers = lapply(1:length(left), function(kkk){
    numbers[left[kkk] < numbers & numbers <= right[kkk]]
  })
  segments = lapply(segment_numbers, function(this){
    bim_just_these[this, ]
  })
  segments
}

# "2016-10-19 12:35:20 PDT" get remove tester_of_segments and return just 
# segments instead of a list containing segments