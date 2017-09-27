#' @export
crocodile = function(fam, Psi, sampling_set, N_sim_reps, weights = NULL){
  params_sampling_set = sampling_set$params_sampling_set
  sampling_set = sampling_set$sampling_set
  params = params_fn(params_sampling_set, fam, Psi, sampling_set, N_sim_reps, weights)
  ooo = uno_experimento_fn(params, calculate_optimized = TRUE)
  ooo
}