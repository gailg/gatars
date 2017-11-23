#' @export
uno_experimento_fn = function(
  adaptive_conf_level, calculate_optimized, g_target, MMM, 
  N_simulated_nulls_interval, N_simulated_nulls_limit, Phi, sampling_set, theta_init, WWW, y_1, y_2
  ){
  answer = if (rankMatrix(g_target) < MMM) {
    list(message = "error--g_target not full rank")
  } else {
    bo = basic_and_optimized_lu_fn(g_target, Phi, theta_init, WWW, y_1, y_2)
    p_value_basic = bo$p_value_basic
    p_value_basic = bo$p_value_basic
    q_basic = bo$q_basic
    q_optimized = bo$q_optimized
    qqq = c(q_basic, q_optimized)
    theta = bo$theta
    x_observed = bo$xxx
    x_observed
    #--------------------------------------------------- p_value_optimized
    ooo = p_value_optimized_fn(
      adaptive_conf_level, calculate_optimized, MMM, 
      N_simulated_nulls_interval, N_simulated_nulls_limit, Phi, sampling_set,
      theta, WWW, x_observed, y_1, y_2
    )
    so_far_so_good = ooo$so_far_so_good
    N_simulated_nulls_required = ooo$N_simulated_nulls_required
    p_value_optimized = ooo$p_value_optimized
    #-------------------------------------------------------------- output
    if(so_far_so_good){
      p_value = c(p_value_basic, p_value_optimized)
      list(
        N_simulated_nulls_required = N_simulated_nulls_required,
        p_value = p_value,
        q = qqq,
        x = x_observed)
    } else {
      list(message = "could not obtain full-rank g_target_sim matrices")
    }
  }
  answer
}