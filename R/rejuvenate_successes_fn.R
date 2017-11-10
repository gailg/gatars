#' @export
rejuvenate_successes_fn = function(
  adaptive_conf_level,
  N_simulated_nulls_limit,
  N_simulated_nulls_required,
  optimized_names,
  simulated,
  so_far_so_good,
  successes,
  test_size,
  x_observed
){
  if(so_far_so_good){
    more = unlist(sapply(optimized_names, function(this){
      sum(simulated[[this]] > x_observed[[this]])
    }, simplify = FALSE))
    successes = successes + more
    N_simulated_nulls_required = N_simulated_nulls_required + nrow(simulated)
    ambiguous = any(sapply(successes, function(xxx){
      ci = prop.test(xxx,
                     N_simulated_nulls_required,
                     conf.level = adaptive_conf_level)$conf.int
      ci[1] <= .10
    }))
    still_looking = ambiguous && N_simulated_nulls_required < N_simulated_nulls_limit
  } else { # prepare to jetison
    still_looking = FALSE
  }
  list(
    N_simulated_nulls_required = N_simulated_nulls_required,
    successes = successes,
    still_looking = still_looking)
}