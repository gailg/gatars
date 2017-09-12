#' @export
one_experiment_fn = function(
  params,
  calculate_fancy = TRUE
){
  theta_init = params$theta_init
  g_target = params$g_target
  MMM = params$MMM
  WWW = params$WWW
  www = params$www
  y_1 = params$yyy
  y_2 = params$e_y
  rho_uni = params$rho_uni
  Phi = params$Phi
  answer = if (rankMatrix(g_target) < MMM) {
    list(message = "error--g_target not full rank")
  } else {
    xxx = davies_not_lu_depends_on_g_target_fn(
      g_target, MMM, rho_uni, Phi, theta_init, WWW, www, y_1, y_2)
    p_value_straight = xxx$p_value
    observed = xxx$qqq
    fancy_names = names(observed)[!(names(observed) %in% names(p_value_straight))]
    observed_fancy = observed[, fancy_names]
    #----------------------------------------------------------------p_value_fancy
    fff = p_value_fancy_fn(params,
                           calculate_fancy, fancy_names, MMM, observed_fancy, Phi, rho_uni,
                           theta_init, WWW, www, y_1, y_2)
    so_far_so_good = fff$so_far_so_good
    N_sim_reps_required = fff$N_sim_reps_required
    p_value_fancy = fff$p_value_fancy
    #----------------------------------------------------------------p_value_fancy
    if(so_far_so_good){
      p_value = cbind(p_value_straight, t(p_value_fancy))
      names(observed) = short_names_fn(names(observed))
      names(p_value) = short_names_fn(names(p_value))
      list(
        N_sim_reps_required = N_sim_reps_required,
        observed = observed,
        p_value = p_value)
    } else {
      list(message = "could not obtain full-rank g_target_sim matrices")
    }
  }
}
