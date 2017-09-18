

larrys_example_genotype_etc_fn = function(params_example){
  N_hap_rows = params_example$N_hap_rows
  
  p_sampling_set_generator = params_example$p_sampling_set_generator
  NNN = params_example$NNN
  p_subset_genotype = params_example$p_subset_genotype
  bim_original = params_example$bim_original
  genotype_bim_original = params_example$genotype_bim_original
  dimnames_genotype_2 = bim_original$snp
  target_markers = genotype_bim_original$snp
  
  haplotypes = do.call(cbind, lapply(p_sampling_set_generator, function(ppp){
    rbinom(N_hap_rows, 1, ppp)
  }))
  null_model_chromosomes = larrys_one_null_model_chromosome_fn(params_example)
  genotype_0 = do.call(rbind, lapply(1:nrow(null_model_chromosomes), function(kkk){
    this_row = null_model_chromosomes[kkk, ]
    haplotypes[this_row[1], ] + haplotypes[this_row[2], ]
  }))
  dimnames(genotype_0) = list(1:NNN, dimnames_genotype_2)
  genotype_1 = genotype_0[, colSums(genotype_0) > 0]
  choose_these = dimnames(genotype_1)[[2]] %in% target_markers |
    as.logical(rbinom(1:ncol(genotype_1), 1, p_subset_genotype))
  genotype = genotype_1[, choose_these]
  bim = bim_original[bim_original$snp %in% dimnames(genotype)[[2]], ]
  # genotype_dosage = dosage[, dimnames(dosage)[[2]] %in%  dimnames_genotype_dosage_2]
  # genotype_bim = genotype_bim_original[genotype_bim_original$snp %in% dimnames(genotype_dosage)[[2]], ]
  list(bim = bim,
       genotype = genotype,
       target_markers = target_markers)
}







larrys_one_null_model_chromosome_fn = function(params_example, y_needed = FALSE){
  N_affected_sib_pairs = params_example$N_affected_sib_pairs
  N_independent_cases = params_example$N_independent_cases
  N_independent_controls = params_example$N_independent_controls
  N_unaffected_sib_pairs = params_example$N_unaffected_sib_pairs
  affected_sib_pairs = if(N_affected_sib_pairs == 0){
    NULL
  } else {
    sib_pairs_zero_fn(params_example, N_affected_sib_pairs)$haplotype_indices
  }
  independent_cases = if(N_independent_cases == 0){
    NULL
  } else {
    independent_people_zero_fn(params_example, N_independent_cases)$haplotype_indices
  }
  independent_controls = if(N_independent_controls == 0){
    NULL
  } else {
    independent_people_zero_fn(params_example, N_independent_controls)$haplotype_indices
  }
  unaffected_sib_pairs = if(N_unaffected_sib_pairs == 0){
    NULL
  } else {
    sib_pairs_zero_fn(params_example, N_unaffected_sib_pairs)$haplotype_indices
  }
  haplotype_indices = rbind(
    affected_sib_pairs, 
    independent_cases, 
    unaffected_sib_pairs,
    independent_controls)
  rownames(haplotype_indices) = NULL
  haplotype_indices
} 

larrys_example_fam_fn = function(params_example, g_target){
  N_affected_sib_pairs = params_example$N_affected_sib_pairs
  N_independent_controls = params_example$N_independent_controls
  age_generator = params_example$age_generator
  beta_nongenetic = params_example$beta_nongenetic
  causal_snps = params_example$causal_snps
  beta_causal_snps = params_example$beta_causal_snps
  NNN = params_example$NNN
  
  N_FAMILY = N_affected_sib_pairs * 2
  N_CASE_CONTROL = N_independent_controls
  mecc = c(rep("FAMILY", N_FAMILY), rep("CASE_CONTROL", N_CASE_CONTROL))
  family = ifelse(mecc == "FAMILY", 1, 0)
  
  age = round(c(
    rnorm(N_FAMILY, mean = age_generator["FAMILY", ]$mean, sd = age_generator["FAMILY", ]$sd),
    rnorm(N_CASE_CONTROL, mean = age_generator["CASE_CONTROL", ]$mean, sd = age_generator["CASE_CONTROL", ]$sd)))
  
  beta_1 = c(beta_nongenetic, beta_causal_snps)
  xxx = as.matrix( data.frame(intercept = 1, age, family, interaction = family * age) )
  x_1 = cbind(xxx, g_target[, causal_snps])
  p_logistic = exp(x_1 %*% beta_1) / (1 + exp(x_1 %*% beta_1))
  yyy = rbinom(NNN, 1, p_logistic)
  df = data.frame(age, family, y = yyy)
  logistic_regression = glm(y ~ age * family, family = binomial, data = df)
  
  e_y = logistic_regression$fitted.values
  if(FALSE){
    summary(logistic_regression)
    col = ifelse(yyy == 1, "black", "gray")
    plot(age, p_logistic, ylim = c(0, 1), pch = substr(mecc, 1, 1), cex = .6, col = col,
         main = "F = family, M = CASE_CONTROL; black signifies affected, gray unaffected")
    col = ifelse(yyy == 1, "red", "pink")
    points(age, e_y,  ylim = c(0, 1), pch = substr(mecc, 1, 1), cex = .6, col = col)
    data.frame(mean_y_FAMILY = mean(yyy[mecc == "FAMILY"]),
               mean_y_CASE_CONTROL = mean(yyy[mecc == "CASE_CONTROL"]))
  }
  fam = data.frame(id = 1:NNN, y = yyy, e_y)
  fam
}


