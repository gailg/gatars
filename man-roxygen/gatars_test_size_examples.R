#' @examples
#' # Preparing the data
#' bim = gatars_example$bim
#' exclusion_region = gatars_example$exclusion_region
#' genotype = gatars_example$genotype
#' phenotype = gatars_example$phenotype
#' Psi = gatars_example$Psi
#' target_markers = gatars_example$target_markers[3:5]
#' 
#' # Show figure to illustrate Psi
#' NNN = nrow(phenotype)
#' first_ten = 1:10
#' last_ten = NNN - (9:0)
#' matrix_image_fn(Psi[c(first_ten, last_ten), c(first_ten, last_ten)],
#'                 main = "First and last 10 rows and columns of Psi")
#'                 
#' # Check the rank of the genotype_target_markers matrix
#' library(Matrix)
#' genotype_target_markers = genotype[, target_markers]
#' list(target_markers = target_markers,
#'      rank = as.numeric(rankMatrix( genotype_target_markers)))
#'      
#' # Call gatars_sampling_set to create sampling_set
#' epsilon = 0.01
#' sampling_set = gatars_sampling_set(
#'     bim, epsilon, exclusion_region,
#'     genotype, hotspot, target_markers)
#' print(sampling_set)
#' 
#' # Call gatars_test_size using an N_simulated_nulls which is way too small
#' N_simulated_nulls = 10
#' gatars_test_size(phenotype, Psi, sampling_set, N_simulated_nulls)
