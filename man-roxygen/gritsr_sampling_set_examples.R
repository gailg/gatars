#' @examples
#' bim = alternative_example$bim
#' genotype = alternative_example$genotype
#' fam = alternative_example$fam
#' target_markers = alternative_example$target_markers[3:5]
#' exclusion_region = data.frame(chromosome = integer(), begin = integer(), end = integer())
#' Psi = alternative_example$Psi
#' target_markers = target_markers
#' 
#' NNN = nrow(fam)
#' first_ten = 1:10
#' last_ten = NNN - (9:0)
#' matrix_image_fn(Psi[c(first_ten, last_ten), c(first_ten, last_ten)],
#'                 main = "First and last 10 rows and columns of Psi")
#'                 library(Matrix)
#'                 
#' genotype_target_markers = genotype[, target_markers]
#' list(target_markers = target_markers,
#'      rank = as.numeric(rankMatrix( genotype_target_markers)))
#'      
#' set.seed(2)
#' sampling_set = gritsr_sampling_set(
#'   bim, genotype, target_markers, exclusion_region, hotspot, 
#'   epsilon_on_log_scale = 0.02
#' )
#' print(sampling_set)
#' names(sampling_set)
#' str(sampling_set$sampling_set)
