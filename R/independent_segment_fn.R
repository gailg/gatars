#' @title Chop a chromosome into independent segments
#' 
#' @description (a) Subset \code{bim} to just those
#' inside \code{chromosome}, (b) discard the head and tail (those bp falling 
#' outside the first and last hotspots of \code{chromosome}),
#' and (c) using \code{hotspot} cut up the remaining \code{bim} 
#' into independent segments.
#' 
#' @param bim The bim object described in \code{link{gatars_sampling_set}}.
#' 
#' @param chromosome An integer equal to the chromosome number (1 through 22)
#' of the chromosome you want to chop up.
#' 
#' @param hotspot The hotspot object described in \code{link{gatars_sampling_set}}.
#' 
#' @return A list of \code{data.frame}'s, one \code{data.frame} for each
#' indpendent segment.  Each \code{data.frame} contains one or few rows from
#' \code{bim}, and has columns \code{chromosomse}, \code{snp}, and \code{bp}.
#' 
#' @examples 
#' bim = gatars_example$bim
#' chromosome = 21 # I chose 21 because in our example it is the shortest
#' table(bim$chromosome)
#' independent_segment_fn_answer = independent_segment_fn(bim, chromosome, hotspot)
#' length(independent_segment_fn_answer)
#' head(independent_segment_fn_answer)
#' 
#' @export
independent_segment_fn = function(bim, chromosome, hotspot){
  bim_this_chromosome = bim[bim$chromosome == chromosome, ]
  hot = hotspot_this_chromosome = hotspot[hotspot$chromosome == chromosome, ]
  in_a_hotspot_q = sapply(bim_this_chromosome$bp, function(this){
    any(hot$start <= this & this <= hot$end)
  })
  # v-0.2.8 I need to worry about not having any bim left after all the snps
  # inside hotspots are removed
  if(all(in_a_hotspot_q)){
    NULL
  } else {
    # bim_just_these removes snps that are inside hotspots
    bim_just_these = bim_this_chromosome[!in_a_hotspot_q, ]
    bp = bim_just_these$bp
    # hot_just_these removes head and tail rows that do no cutting
    hot_just_these = hot[bp[1] <= hot$center & hot$center <= bp[length(bp)], ]
    center = hot_just_these$center
    # bp_index is the index of the last snp left of this hot$center
    # Saying it with more details, for each center of a hotspot, 
    # give me the index of the bp immediately to the center's left.
    # Repeats indicate that there were no bp's between hotspots.
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
}

# "2016-10-19 12:35:20 PDT" get remove tester_of_segments and return just 
# segments instead of a list containing segments