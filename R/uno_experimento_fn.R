#' @title The workhorse of \code{gatars_test_size}
#' 
#' @description \code{uno_experimento_fn} begins by calling
#' \code{basic_and_optimized_lu_fn} to obtain the p-values of the
#' basic statistic (\code{bo$p_value_basic}) and \code{x_observed}
#' (\code{bo$x_observed}), and then calls \code{p_value_optimized_fn}
#' to obtain the p-values of the optimized statistics
#' (\code{pvo$p_value_optimized}).
#' 
#' @param adaptive_conf_level \code{params_fn} sets this real number to 0.99.  
#' This if the level of confidence used to to decide if the number of 
#' simulated nulls obtained so far will not be significant.  
#' 
#' @param calculate_optimized A logical.  Set this to \code{TRUE} if genome 
#' resampling is desired (when analyzing real data),   and to \code{FALSE} when
#' calculating the power in the simulations where simulations rather than
#' genome resampling can be used to obtain simulated nulls.
#' 
#' @param g_target A numerical matrix with \code{NNN} rows and \code{MMM} 
#' columns equal to the genotype matrix in the manuscript.
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
#' @param theta_init A vector of length \code{2} that is the initial value of the
#' reparametrization of alpha when I am finding
#' minimum p-value in the full triangle \code{(alpha_B, alpha_S, alpha_T)}.
#' 
#' @param WWW A diagonal (numerical) matrix of dimension \code{MMM} by \code{MMM}
#' with the diagonals equal to the \code{weights}.  (The user will specify
#' \code{weights} in her call to \code{gatars_test_size}.)
#' 
#' @param y_1 A numerical vector of length \code{NNN} equal to what is referred
#' to in the manuscript as \code{y}, the vector of subjects' coded trait 
#' phenotypes.
#' 
#' @param y_2 A numerical vector of length \code{NNN} equal to what is referred
#' to in the manuscript as \code{mu}, the vector of user-specified phenotype
#' predictions. 
#' 
#' @return A list containing
#' \itemize{
#' \item{\code{N_simulated_nulls_required}: } {
#' A positive integer equal to the number of simulated nulls used to calculate the
#' p_value_optimized.
#' }
#' \item{\code{p_value}: } {
#' A named numerical vector containing the p-values of each of the seven statistics.
#' }
#' \item{\code{q}: } {
#' A named numerical vector containing \code{Q_B}, \code{Q_S}, \code{Q_T}, 
#' \code{Q_BS}, \code{Q_BT}, \code{Q_ST}, and \code{Q_BST}, the quadratic form
#' versions of the seven statistics.
#' }
#' \item{\code{x}: } {
#' A named numerical vector containing the x-transformations of the four optimized
#' statistics.
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
#' epsilon = 0.01
#' exclusion_region = NULL
#' sampling_set = gatars_sampling_set(
#'   bim, epsilon, exclusion_region,
#'   genotype, hotspot, target_markers)
#' print(sampling_set)
#' sampling_set = sampling_set$sampling_set
#' adaptive_conf_level = 0.99
#' calculate_optimized = TRUE
#' g_target = genotype[, target_markers]
#' MMM = ncol(g_target)
#' N_simulated_nulls_interval = 20
#' N_simulated_nulls_limit = 100
#' y_1 = yyy = phenotype$y
#' y_2 = mu = phenotype$mu
#' Phi = Phi_fn(Psi, y_1, y_2)
#' theta_init = rep(pi/3, 2)
#' www_num = rep(1, MMM)
#' www = www_num/sum(www_num) * MMM
#' WWW = diag(www)
#' ooo = uno_experimento_fn(
#'   adaptive_conf_level, calculate_optimized, g_target, MMM, 
#'   N_simulated_nulls_interval, N_simulated_nulls_limit, 
#'   Phi, sampling_set, theta_init, WWW, y_1, y_2)
#' ooo
#' 
#' @export
uno_experimento_fn = function(
  adaptive_conf_level, calculate_optimized, g_target, MMM, 
  N_simulated_nulls_interval, N_simulated_nulls_limit, Phi, 
  sampling_set, theta_init, WWW, y_1, y_2, 
  statistics = c("BS", "BT", "ST", "BST")
  ){
  answer = if (rankMatrix(g_target) < MMM) {
    list(message = "error--g_target not full rank")
  } else {
    bo = basic_and_optimized_lu_fn(g_target, Phi, theta_init, WWW, y_1, y_2, statistics)
    p_value_basic = bo$p_value_basic
    p_value_basic = bo$p_value_basic
    q_basic = bo$q_basic
    q_optimized = bo$q_optimized
    qqq = c(q_basic, q_optimized)
    theta = bo$theta
    x_observed = bo$x_observed
    x_observed
    #--------------------------------------------------- p_value_optimized
    ooo = p_value_optimized_fn(
      adaptive_conf_level, calculate_optimized, MMM, 
      N_simulated_nulls_interval, N_simulated_nulls_limit, Phi, sampling_set,
      theta, WWW, x_observed, y_1, y_2, statistics
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