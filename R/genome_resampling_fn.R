#' @export
genome_resampling_fn = function(
  alpha_uni,
  MMM,
  N_sim_reps_interval,
  optimized_names,
  Phi,
  sampling_set,
  theta, 
  WWW,
  y_1,
  y_2
){
  simulated = data.frame(rep(NA, N_sim_reps_interval), 
                         rep(NA, N_sim_reps_interval), 
                         rep(NA, N_sim_reps_interval), 
                         rep(NA, N_sim_reps_interval))
  names(simulated) = optimized_names
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
      bo = basic_and_optimized_lu_fn(alpha_uni, genotype_sim, Phi, theta, WWW, y_1, y_2)
      one_row_in_simulated = bo$xxx
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