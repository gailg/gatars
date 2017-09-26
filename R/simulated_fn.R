#' @import Matrix
#' @export
simulated_fn = function(params,
  fancy_names,
  MMM,
  Phi,
  rho_uni,
  WWW,
  www,
  y_1,
  y_2
){
  sampling_set = params$sampling_set
  N_sim_reps_interval = params$N_sim_reps_interval
  simulated = data.frame(standard_lee = rep(NA, N_sim_reps_interval),
                         optim_skat = rep(NA, N_sim_reps_interval),
                         odd = rep(NA, N_sim_reps_interval),
                         optim = rep(NA, N_sim_reps_interval))
  so_far_so_good = TRUE
  n_sim = 1
  while(so_far_so_good & n_sim <= N_sim_reps_interval){
    still_looking_for_genotype_sim = TRUE
    N_bad_genotype_sim = 0
    while(still_looking_for_genotype_sim){
      genotype_sim = do.call(cbind, lapply(sampling_set, function(this_snp){
        this_snp[, sample(1:ncol(this_snp), 1)]
      }))
      if(rankMatrix(genotype_sim) == MMM) {
        good_genotype_sim = TRUE
        still_looking_for_genotype_sim = FALSE
      } else {
        print("genotype_sim not full rank")
        N_bad_genotype_sim = N_bad_genotype_sim + 1
        if(N_bad_genotype_sim > 1000) {
          good_genotype_sim = FALSE
          still_looking_for_genotype_sim = FALSE
        }
      }
    }
    if(good_genotype_sim){
      xxx = davies_not_lu_depends_on_g_target_fn(
        genotype_sim, MMM, rho_uni, Phi, WWW, www, y_1, y_2)
      one_row_in_simulated = xxx$qqq[, fancy_names]
      simulated[n_sim, ] = one_row_in_simulated
      n_sim = n_sim + 1
    } else {
      so_far_so_good = FALSE
    }
  }
  list(
    simulated = simulated,
    so_far_so_good = so_far_so_good
  )
}
