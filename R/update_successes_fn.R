#' @export
update_successes_fn = function(params,
  fancy_names,
  N_sim_reps_required,
  observed_fancy,
  simulated,
  so_far_so_good,
  successes
){
  adaptive_conf_level = params$adaptive_conf_level
  N_sim_reps_limit = params$N_sim_reps_limit
  test_size = params$test_size
  if(so_far_so_good){
    more = unlist(sapply(fancy_names, function(this){
      sum(simulated[[this]] > observed_fancy[[this]])
    }, simplify = FALSE))
    successes = successes + more
    N_sim_reps_required = N_sim_reps_required + nrow(simulated)
    ambiguous = any(sapply(successes, function(xxx){
      ci = prop.test(xxx,
                     N_sim_reps_required,
                     conf.level = adaptive_conf_level)$conf.int
      ci[1] <= test_size && test_size <= ci[2]
    }))
    still_looking = ambiguous && N_sim_reps_required < N_sim_reps_limit
  } else { # prepare to jetison
    still_looking = FALSE
  }
  list(
    N_sim_reps_required = N_sim_reps_required,
    successes = successes,
    still_looking = still_looking)
}
