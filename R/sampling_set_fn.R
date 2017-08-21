#' @export
sampling_set_fn = function(params_sampling_set){
  bim_with_target_regions_removed = bim_with_target_regions_removed_fn(params_sampling_set)
  p_target = params_sampling_set$p_target
  MMM = params_sampling_set$MMM
  epsilon_on_log_scale = params_sampling_set$epsilon_on_log_scale
  big_genotype = genotype[, bim_with_target_regions_removed$snp]
  p_sim = unname(colMeans(big_genotype)/2)
  sampling_set = lapply(1:MMM, function(mmm){
    p_this = p_target[mmm]
    lower = p_this * (1 - epsilon_on_log_scale)
    upper = p_this * (1 + epsilon_on_log_scale)
    good_ones_q = (lower < p_sim & p_sim < upper)
    good_ones = which(good_ones_q)
    N_good_ones = length(good_ones)
    if(N_good_ones > 1000){
      good_ones = sort(sample(good_ones, 1000))
    }
    answer =  big_genotype[, good_ones, drop = FALSE]
    answer
  })
  report = do.call(rbind, lapply(1:MMM, function(mmm){
    this_genotype = sampling_set[[mmm]]
    p_sim = unname(colMeans(this_genotype)/2)
    data.frame(min = round(min(p_sim), 5),
               p_target = round(p_target[mmm], 5),
               max = round(max(p_sim), 5),
               set_size = ncol(this_genotype))
  }))
  rownames(report) = NULL
  report
  list(
    params_sampling_set = params_sampling_set,
    report = report,
    sampling_set = sampling_set)
}
