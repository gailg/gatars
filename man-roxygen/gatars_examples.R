#' @examples
#' bim = alternative_example$bim
#' genotype = alternative_example$genotype
#' fam = alternative_example$fam
#' target_markers = alternative_example$target_markers[3:5]
#' Psi = alternative_example$Psi
#' exclusion_region = alternative_example$exclusion_region
#' 
#' NNN = nrow(fam)
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
#' epsilon_on_log_scale = 0.02
#' sampling_set = gatars_sampling_set(
#'     bim, epsilon_on_log_scale, exclusion_region,
#'     genotype, hotspot, target_markers)
#' print(sampling_set)
#' 
#' N_sim_reps = 10
#' gatars(fam, Psi, sampling_set, N_sim_reps = 10)
