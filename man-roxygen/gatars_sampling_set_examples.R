#' @examples
#' bim = gatars_example$bim
#' genotype = gatars_example$genotype
#' phenotype = gatars_example$phenotype
#' target_markers = gatars_example$target_markers[3:5]
#' exclusion_region = gatars_example$exclusion_region
#' Psi = gatars_example$Psi
#' target_markers = target_markers
#' 
#' NNN = nrow(phenotype)
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

