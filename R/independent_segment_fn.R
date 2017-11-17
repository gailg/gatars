#' @title Chop a chromosome into independent segments
#' 
#' @description (a) Subset \code{bim} to just those
#' inside \code{chromosome}, (b) discard any rows of sim
#' that fall inside a hotspot,
#' and (c) using \code{hotspot} cut up the remaining \code{bim} 
#' into independent segments.
#' 
#' @param bim The bim object described in \code{link{gatars_sampling_set}}.
#' 
#' @param chromosome An integer equal to the chromosome number (1 through 22)
#' of the chromosome you want to chop up.
#' 
#' @param hotspot The hotspot object described in \code{\link{gatars_sampling_set}}.
#' 
#' @return A list of \code{data.frame}'s, one \code{data.frame} for each
#' indpendent segment.  Each \code{data.frame} contains one or few rows from
#' \code{bim}, and has columns \code{chromosomse}, \code{snp}, and \code{bp}.
#' 
#' @examples 
#' # Calling independent_segment_fn on a "real" example".  
#' # I will do it on chromosome = 21 because it is the shortest chromosome
#' bim = gatars_example$bim
#' chromosome = 21 # 
#' table(bim$chromosome)
#' independent_segment_fn_answer = independent_segment_fn(bim, chromosome, hotspot)
#' length(independent_segment_fn_answer)
#' head(independent_segment_fn_answer)
#' #
#' # ---------------Use bim_teeny and hot_teeny for a more useful example
#' # -------------- Get ready to call the independent_segment_fn
#' bim = gatars_example$bim
#' bim_for_chromosome_21 = bim[bim$chromosome == 21, ]
#' hot_for_chromosome_21 = hotspot[hotspot$chromosome == 21, ]
#' head(bim_for_chromosome_21, 20)
#' head(hot_for_chromosome_21, 50)
#' bim_teeny = bim_for_chromosome_21[1:13, ]
#' hot_teeny = hot_for_chromosome_21[hot_for_chromosome_21$center %in% c(13227919, 15561921, 16361420), ]
#' row.names(bim_teeny) = NULL
#' row.names(hot_teeny) = NULL
#' list(bim_teeny = bim_teeny, 
#'      hot_teeny = hot_teeny)
#' # -------------- hot_q tells me which bim row falls in which hot row
#' hot_q = do.call(rbind, lapply(1:nrow(bim_teeny), function(jjj){
#'   hot_teeny$start <=  bim_teeny$bp[jjj] & bim_teeny$bp[jjj] <= hot_teeny$end
#' }))
#' dimnames(hot_q) = list(
#'   paste0("bim_", 1:nrow(hot_q)), paste0("hot_", 1:ncol(hot_q)))
#' hot_q
#' # -------------- and I see that bim_teeny row 3 falls in hot_teeny row 2
#' bim_teeny_rows_in_hotspot = apply(hot_q, 1, any)
#' bim_teeny_rows_in_hotspot
#' bim_teeny_not_in_hotspot = bim_teeny[!bim_teeny_rows_in_hotspot, ]
#' list(bim_teeny_not_in_hotspot = bim_teeny_not_in_hotspot,
#'      hot_teeny = hot_teeny)
#' # -------------- From inspection, I see that the independent segments should be
#' # -------------- 1:2, 3, 5:12, and 13
#' independent_segment_fn(bim_teeny, 21, hot_teeny)
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