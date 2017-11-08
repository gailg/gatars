#' @export
uno_experimento_fn = function(params, calculate_optimized){
  adaptive_conf_level = params$adaptive_conf_level
  alpha_uni = params$alpha_uni
  g_target = params$g_target
  MMM = params$MMM
  N_sim_reps_interval = params$N_sim_reps_interval
  N_sim_reps_limit = params$N_sim_reps_limit
  theta_init = params$theta_init
  WWW = params$WWW
  y_1 = params$yyy
  y_2 = params$mu
  Phi = params$Phi
  sampling_set = params$sampling_set
  test_size = params$test_size
  theta_init = params$theta_init
  answer = if (rankMatrix(g_target) < MMM) {
    list(message = "error--g_target not full rank")
  } else {
    bo = basic_and_optimized_lu_fn(alpha_uni, g_target, Phi, theta_init, WWW, y_1, y_2)
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
      adaptive_conf_level, alpha_uni, calculate_optimized, MMM, 
      N_sim_reps_interval, N_sim_reps_limit, Phi, sampling_set,
      test_size, theta, WWW, x_observed, y_1, y_2
    )
    so_far_so_good = ooo$so_far_so_good
    N_sim_reps_required = ooo$N_sim_reps_required
    p_value_optimized = ooo$p_value_optimized
    #-------------------------------------------------------------- output
    if(so_far_so_good){
      p_value = c(p_value_basic, p_value_optimized)
      list(
        N_sim_reps_required = N_sim_reps_required,
        p_value = p_value,
        q = qqq,
        x = x_observed)
    } else {
      list(message = "could not obtain full-rank g_target_sim matrices")
    }
  }
  answer
}