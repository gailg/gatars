#' @title Estimate the p-values of the optimized statistics
#' 
#' @description At its most basic level, \code{p_value_optimized_fn} generates a 
#' "large number" of simulated null statistics to compare to the observed.
#' The p-value is the proportion of these simulated nulls that exceed the observed.
#' The user can specify how big to make "a large number".  For the simulations, I 
#' added another level of control:  If after generating a moderate number of
#' simulated nulls, it becomes obvious that the observed of all four optimized
#' statistics is not going to be significant, there is no longer any need to
#' continue to the big number that might be needed if at least one of the 
#' significance levels reported by the moderate number of simulated nulls is close
#' to the desired significance level.
#' 
#' To specify "a large number, use \code{N_simulated_nulls_limit}.
#' To allow \code{p_value_optimized_fn} to abort when it becomes obvious that the 
#' four optimized statistics will not be significant, use
#' \code{N_simulated_nulls_interval} to specify a moderate number on which to judge
#' whether it is safe to abort.
#' \code{gatars_test_size} allows the user to set \code{N_simulated_nulls} to be
#' the very large number and internally sets 
#' \code{N_simulated_nulls_interval} equal to \code{N_simulated_nulls_limit}
#' equal to \code{N_simulated_nulls}.
#' 
#' @param adaptive_conf_level params_fn sets this real number to 0.99.  This if the 
#' level of confidence used to to decide if the number of simulated nulls obtained
#' so far will not be significant.  
#' 
#' @param calculate_optimized A logical.  Set this equal to \code{TRUE} if genome 
#' resampling is desired.  I set \code{calculate_optimized} to \code{FALSE} when
#' calculating the power in the simulations where I used simulations rather than
#' genome resampling to obtain simulated nulls.
#' 
#' @param MMM A positive integer equal to the length of \code{target_markers}
#' in \code{gatars_sampling_set}
#' 
#' @param N_simulated_nulls_interval A positive integer equal to the number
#' of rows of simulated optimized statistics requested from each call to 
#' \code{genome_resampling_fn}.  This is equal to "the moderate number of simulated
#' nulls" to use to get an idea of the significance level before deciding if
#' bigger numbers of simulated nulls will be required.
#' 
#' @param N_simulated_nulls_limit A positive integer equal to the absolute largest
#' number of simulated nulls to generate to estimate the p-values of the observed
#' statistics.
#' 
#' @param Phi A numerical matrix of dimension \code{2} by \code{2}.
#' \code{Phi_{k_1, k_2} = y_{k_1} Psi y_{k_2}}.  
#' This  matrix is a useful intermediate calculation for getting
#' \code{V_z}: \code{V_z = kronecker(Phi, W_VG_W)}. 
#' It is of dimension \code{2} by \code{2} because there are two entities 
#' \code{y_1} and \code{y_2}.
#' 
#' @param sampling_set #' A list of \code{MMM} matrices, one matrix for each target marker.
#' The \code{mmm}-th matrix is the sampling set for the \code{mmm}-th
#' target marker and has \code{NNN} rows and up to \code{1000} columns, each 
#' column containing a column from \code{genotype}.  These columns
#' do not intersect with any of the target markers or exclusion regions
#' and the minor allele frequencies of the columns in \code{mmm}-th sampling set match
#' the minor allele frequency of the \code{mmm}-th target marker. One of the objects
#' returned by \code{gatars_sampling_set}.
#' 
#' 
#' 
#' @export
p_value_optimized_fn = function(
  adaptive_conf_level,
  calculate_optimized,
  MMM,
  N_simulated_nulls_interval,
  N_simulated_nulls_limit,
  Phi,
  sampling_set,
  theta,
  WWW,
  x_observed,
  y_1,
  y_2
){
  optimized_names = names(x_observed)
  if(calculate_optimized){
    # successes counts for each of the four optimized statistics the number of simulated nulls
    # that exceeded the corresponding observed statistic.
    # N_simulated_nulls_required counts the number of simulated nulls generated so far
    successes = rep(0, length(x_observed))
    names(successes) = names(x_observed)
    N_simulated_nulls_required = 0
    still_looking = TRUE
    # Each iteration through the while loop calls genome_resampling_fn to generate 
    # for each of the four optimized statistics N_simulated_nulls_interval simulated nulls 
    # and records their x-values in sss.
    # Then it calls rejuvenate_successes_fn to decide which simulated nulls exceeded the
    # observed statistics, adds these to the successes vector.
    # N_simulated_nulls_required is also updated.
    # The loop continues until it is obvious that all four statistics are not significant
    # or until N_simulated_nulls_required exceeds N_simulated_nulls_limit
    while(still_looking){
      sss = genome_resampling_fn(MMM, N_simulated_nulls_interval, optimized_names, Phi, 
                                 sampling_set, theta, WWW, y_1, y_2)
      
      simulated = sss$simulated
      so_far_so_good = sss$so_far_so_good
      uuu = rejuvenate_successes_fn(adaptive_conf_level, N_simulated_nulls_limit, N_simulated_nulls_required,
                                    optimized_names, simulated, so_far_so_good, successes, x_observed)
      N_simulated_nulls_required = uuu$N_simulated_nulls_required
      successes = uuu$successes
      still_looking = uuu$still_looking
    }
    # The p-value is the proportion of the simulated nulls that were a success.
    p_value = if(so_far_so_good){
      successes / N_simulated_nulls_required
    } else {
      list(message = "could not obtain full-rank genotype_sim matrices")
    }
  } else {
    so_far_so_good = TRUE
    N_simulated_nulls_required = 0
    p_value = sapply(optimized_names, function(dummy) -1)
  }
  list(
    so_far_so_good = so_far_so_good,
    N_simulated_nulls_required = N_simulated_nulls_required,
    p_value_optimized = p_value)
}