#' @title Given a new set of simulated nulls, update \code{successes} and \code{N_simulated_nulls_required}
#' 
#' @description The innards of \code{p_value_optimized_fn} contains a \code{while}
#' loop that calls \code{genome_resampling_fn} to obtain \code{N_simulated_nulls_interval}
#' replications of simulated nulls.  \code{rejuvenate_successes_fn} counts for each statistic 
#' the number of successes in the new replications and undates \code{successes} and
#' \code{N_simulated_nulls_required}. In the end, the p-value is \code{successes} divided
#' by \code{N_simulated_nulls_required}.
#' 
#' @param adaptive_conf_level \code{params_fn} sets this real number to 0.99.  This is the 
#' level of confidence used to to decide if the number of simulated nulls obtained
#' so far will not be significant. 
#' 
#' @param N_simulated_nulls_limit A positive integer equal to the largest
#' number of simulated nulls to generate to estimate the p-values of the observed
#' statistics.
#' 
#' @param N_simulated_nulls_required A positive integer equal to the number of simulated nulls
#' before.
#' 
#' @param simulated A \code{data.frame} with one column for each of the optimized statistics
#' and rows containing replications of simulated nulls.  It is these replications that 
#' \code{rejuvenate_successes_fn} will be using to update \code{successes} and 
#' \code{N_simulated_nulls_required}.
#' 
#' @param so_far_so_good A logical equal to \code{TRUE} if each simulated genotype matrix
#' used to calculate \code{simulated} was obtained with full rank before 1000 bad tries.
#' 
#' @param successes A named vector containing the number of successes of each optimized
#' statistic before.  A simulated null is a success if it exceeds the corresponding
#' observed optimized statistic.
#' 
#' @param x_observed A named numerical vector containing the transformed x of the observed optimized
#' statistics.
#' 
#' @return A list containing the following three objects
#' \itemize{
#' \item{\code{N_simulated_nulls_required}: } {
#' An integer equal to \code{N_simulated_nulls_required} before plus the number of rows
#' in \code{simulated}.
#' }
#' \item{\code{successes}: } {A named numerical vector equal to the \code{successes}
#' before plus the number of successes counted in \code{simulated}.
#' }
#' \item{\code{still_looking}: } {A logical equal to \code{TRUE} if after updating there
#' is still some reasonable probability that at least one of the optimizing statistics
#' is significant and \code{N_simulated_nulls_required} is still lower than
#' \code{N_simulated_nulls_limit}.
#' }
#' }
#' 
#' @examples 
#' library(Matrix)
#' bim = gatars_example$bim
#' genotype = gatars_example$genotype
#' phenotype = gatars_example$phenotype
#' Psi = gatars_example$Psi
#' target_markers = gatars_example$target_markers[3:5]
#' g_target = genotype[, target_markers]
#' MMM = ncol(g_target)
#' NNN = nrow(g_target)
#' e_g_target_1 = colMeans(g_target)
#' p_target = e_g_target_1/2
#' e_g_target = matrix(rep(e_g_target_1, nrow(g_target)), nrow = nrow(g_target), byrow = TRUE)
#' y_1 = yyy = phenotype$y
#' y_2 = mu = phenotype$mu
#' Phi = Phi_fn(Psi, y_1, y_2)
#' www_num = rep(1, MMM)
#' www = www_num/sum(www_num) * MMM
#' WWW = diag(www)
#' zzz_etc = zzz_and_first_two_moments_fn(g_target, Phi, WWW, y_1, y_2)
#' zzz = zzz_etc$zzz
#' mu_z = zzz_etc$mu_z
#' V_z = zzz_etc$V_z
#' AAA = AAA_fn(1, 0, 0, MMM)
#' theta_init = rep(pi/3, 2)
#' statistics = c("BS", "BT", "ST", "BST")
#' bo = basic_and_optimized_lu_fn(g_target, Phi, theta_init, WWW, y_1, y_2, statistics)
#' bo$xxx
#' bo$theta
#' theta = bo$theta
#' x_observed = bo$xxx
#' adaptive_conf_level = 0.01
#' calculate_optimized = TRUE
#' epsilon = 0.01
#' exclusion_region = NULL
#' sampling_set = gatars_sampling_set(
#'   bim, epsilon, exclusion_region,
#'   genotype, hotspot, target_markers)
#' print(sampling_set)
#' sampling_set = sampling_set$sampling_set
#' str(sampling_set)
#' set.seed(1)
#' N_simulated_nulls_limit = 100
#' N_simulated_nulls_interval = 20
#' optimized_names = names(x_observed)
#' 
#' successes = rep(0, length(x_observed))
#' names(successes) = names(x_observed)
#' N_simulated_nulls_required = 0
#' still_looking = TRUE
#' # inside a while loop that continues until still_looking is set to FALSE
#' sss = genome_resampling_fn(MMM, N_simulated_nulls_interval, optimized_names, Phi, 
#'                           sampling_set, theta, WWW, y_1, y_2, statistics)
#' simulated = sss$simulated
#' simulated
#' x_observed
#' x_observed["BST"] < simulated$BST
#' so_far_so_good = sss$so_far_so_good
#' N_simulated_nulls_required
#' successes
#' uuu = rejuvenate_successes_fn(
#'     adaptive_conf_level, N_simulated_nulls_limit, N_simulated_nulls_required,
#'     simulated, so_far_so_good, successes, x_observed)
#' uuu
#' 
#' @export
rejuvenate_successes_fn = function(
  adaptive_conf_level,
  N_simulated_nulls_limit,
  N_simulated_nulls_required,
  simulated,
  so_far_so_good,
  successes,
  x_observed
){
  if(so_far_so_good){
    more = unlist(sapply(names(simulated), function(this){
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