#' @export
p_value_optimized_fn = function(
  adaptive_conf_level,
  alpha_uni, 
  calculate_optimized,
  MMM,
  N_sim_reps_interval,
  N_sim_reps_limit,
  Phi,
  sampling_set,
  test_size,
  theta,
  WWW,
  x_observed,
  y_1,
  y_2
){
  optimized_names = names(x_observed)
  if(calculate_optimized){
    successes = rep(0, length(x_observed))
    names(successes) = names(x_observed)
    N_sim_reps_required = 0
    still_looking = TRUE
    while(still_looking){
      sss = genome_resampling_fn(alpha_uni, MMM, N_sim_reps_interval, optimized_names, Phi, 
                                 sampling_set, theta, WWW, y_1, y_2)
      
      simulated = sss$simulated
      so_far_so_good = sss$so_far_so_good
      uuu = rejuvenate_successes_fn(adaptive_conf_level, N_sim_reps_limit, N_sim_reps_required,
                                    optimized_names, simulated, so_far_so_good, successes, test_size, x_observed)
      N_sim_reps_required = uuu$N_sim_reps_required
      successes = uuu$successes
      still_looking = uuu$still_looking
    }
    p_value = if(so_far_so_good){
      successes / N_sim_reps_required
    } else {
      list(message = "could not obtain full-rank genotype_sim matrices")
    }
  } else {
    so_far_so_good = TRUE
    N_sim_reps_required = 0
    p_value = sapply(optimized_names, function(dummy) -1)
  }
  list(
    so_far_so_good = so_far_so_good,
    N_sim_reps_required = N_sim_reps_required,
    p_value_optimized = p_value)
}