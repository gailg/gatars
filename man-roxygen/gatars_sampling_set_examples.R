#' @examples
#' bim = alternative_example$bim
#' genotype = alternative_example$genotype
#' fam = alternative_example$fam
#' target_markers = alternative_example$target_markers[3:5]
#' exclusion_region = alternative_example$exclusion_region
#' Psi = alternative_example$Psi
#' target_markers = target_markers
#' 
#' NNN = nrow(fam)
#' first_ten = 1:10
#' last_ten = NNN - (9:0)
#' matrix_image_fn(Psi[c(first_ten, last_ten), c(first_ten, last_ten)],
#'                 main = "First and last 10 rows and columns of Psi")
#'                 
#' genotype_target_markers = genotype[, target_markers]
#' library(Matrix)
#' list(target_markers = target_markers,
#'      rank = as.numeric(rankMatrix( genotype_target_markers)))
#'      
#' set.seed(2)

