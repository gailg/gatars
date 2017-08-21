#' @export
p_value_fancy_fn = function(params,
  calculate_fancy,
  fancy_names,
  MMM,
  observed_fancy,
  Phi,
  rho_uni,
  theta_lu,
  WWW,
  www,
  y_1,
  y_2
){
  if(calculate_fancy){
    successes = rep(0, length(fancy_names))
    names(successes) = fancy_names
    N_sim_reps_required = 0
    still_looking = TRUE
    while(still_looking){
      sss = simulated_fn(params,
                         fancy_names, MMM, Phi, rho_uni, theta_lu, WWW, www, y_1, y_2)
      simulated = sss$simulated
      so_far_so_good = sss$so_far_so_good
      uuu = update_successes_fn(params,
                                fancy_names, N_sim_reps_required, observed_fancy, simulated, so_far_so_good, successes)
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
    p_value = sapply(fancy_names, function(dummy) -1)
  }
  list(
    so_far_so_good = so_far_so_good,
    N_sim_reps_required = N_sim_reps_required,
    p_value_fancy = p_value)
}
