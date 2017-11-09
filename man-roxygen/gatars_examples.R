#' @examples
#' bim = gatars_example$bim
#' genotype = gatars_example$genotype
#' phenotype = gatars_example$phenotype
#' target_markers = gatars_example$target_markers[3:5]
#' Psi = gatars_example$Psi
#' exclusion_region = gatars_example$exclusion_region
#' 
#' NNN = nrow(phenotype)
#' first_ten = 1:10
#' last_ten = NNN - (9:0)
#' matrix_image_fn(Psi[c(first_ten, last_ten), c(first_ten, last_ten)],
#'                 main = "First and last 10 rows and columns of Psi")
#' library(Matrix)
#'                 
#' genotype_target_markers = genotype[, target_markers]
#' list(target_markers = target_markers,
#'      rank = as.numeric(rankMatrix( genotype_target_markers)))
#'      
#' set.seed(2)
#' epsilon = 0.01
#' sampling_set = gatars_sampling_set(
#'     bim, epsilon, exclusion_region,
#'     genotype, hotspot, target_markers)
#' print(sampling_set)
#' 
#' N_simulated_nulls = 10
#' gatars(phenotype, Psi, sampling_set, N_simulated_nulls)
