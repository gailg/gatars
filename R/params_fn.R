#' @export
params_fn = function(
  params_sampling_set,
  fam,
  Psi,
  sampling_set,
  N_sim_reps,
  weights = NULL
){
  adaptive_conf_level = 0.99
  N_sim_reps_interval = N_sim_reps
  N_sim_reps_limit = N_sim_reps
  rho_uni = seq(0, 1, by = 0.1)
  test_size = 0.01
  theta_init = rep(pi/4, 3)
  y_1 = yyy = fam$y
  y_2 = e_y = fam$e_y
  Phi = Phi_fn(Psi, y_1, y_2)
  MMM = params_sampling_set$MMM
  www_num = if(!is.null(weights)){
    weights
  } else {
    rep(1, MMM)
  }
  www = www_num/sum(www_num) * MMM
  WWW = diag(www)
  www = t(t(www))
  list(
    adaptive_conf_level = adaptive_conf_level,
    e_y = e_y,
    g_target = params_sampling_set$g_target,
    MMM = params_sampling_set$MMM,
    NNN = params_sampling_set$NNN,
    N_sim_reps_interval = N_sim_reps_interval,
    N_sim_reps_limit = N_sim_reps_limit,
    Phi = Phi,
    rho_uni = rho_uni,
    sampling_set = sampling_set,
    test_size = test_size,
    theta_init = theta_init,
    WWW = WWW,
    www = www,
    yyy = yyy)
}
