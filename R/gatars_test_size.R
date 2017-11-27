#' @title Genetic Association Tests for Arbitrarily Related Subjects
#' 
#' @description 
#' Test the association between a specified set of genetic markers, called target markers, 
#' and a binary or quantitative trait using subjects with any genealogical relationship
#' by calculating the p-values of the following seven statistics:
#' the three basic statistics given by  
#' the squared burden statistic \eqn{Q_B}, the linear kernel (SKAT) statistic
#' \eqn{Q_S}, and the trait-based statistic \eqn{Q_T}, and their four optimized linear 
#' combinations \eqn{Q_BS}, \eqn{Q_BT}, \eqn{Q_ST}, \eqn{Q_BST}.  
#' See \url{http://stanford.edu/~ggong/gatars/} for more details.
#' 
#' @param phenotype A data.frame containing \code{NNN} rows and at least the two columns
#' \code{y} and \code{mu}.  \code{y[nnn]} is the binary or quantitative trait of the 
#' \code{nnn}-th person.  If the trait is binary, possible values are \code{0} for not affected
#' and \code{1} for affected; if the trait is quantitative, possible values are real numbers.
#' \code{mu[nnn]}, a real number, is  a trait prediction of 
#' \code{y[nnn]}.  (For example, \code{mu[nnn]} can be the mean of the \code{NNN} \code{y[nnn]}
#' values or the value predicted from a logistic or linear regression of \code{y[nnn]}
#' on nongenetic covariates and/or on principal components of ancestry.)
#' 
#' @param Psi A real-valued matrix with \code{NNN} rows and \code{NNN} columns.  
#' \code{Psi[n1, n2]} is the 
#' correlation of the genotype at one marker between the \code{n1}-th and \code{n2}-th 
#' subjects. \code{Psi} can be estimated from known family pedigree structures and/or from
#' the subjects' genotypes at markers independent of those in the target set
#' in software packages such as KING (Manichaikul et al 2010).  The diagonal 
#' elements are all unity and two people known to be non-identical-twin full sibs and who
#' have parents known to be completely unrelated have correlation 0.5.
#' 
#' @param sampling_set An object of class \code{gatars_sampling_set} produced by 
#' \code{gatars_sampling_set}. (See \code{\link{gatars_sampling_set}}.)  
#' The object \code{sampling_set} is a list containing the objects
#' \code{params_sampling_set}, \code{report}, and \code{sampling_set}.  The important object
#' is \code{sampling_set$sampling_set}, a list containing \code{MMM} matrices,
#' the \code{mmm}-th matrix containing the sampling set for the \code{mmm}-th target marker.
#' To run genome resampling, \code{gatars_test_size}
#' creates "simulated genotype matrices" by sampling
#' and binding together one column from each of the \code{MMM} matrices.  
#' 
#' @param N_simulated_nulls An integer equal to the number of simulated genotype matrices to generate
#' when estimating the test sizes of the optimized statistics.
#' 
#' @param weights A vector of length \code{MMM} equal to the length of \code{target_markers} which
#' is an input of \code{gatars_sampling_set}. The entries in \code{weights} are non-negative real
#' numbers.  The size of \code{weights[mmm]} reflects the importance of the \code{mmm}-th target
#' marker.  If \code{weights} is not specified or set equal to \code{NULL}, \code{gatars_test_size}
#' assumes you would like equal weights among the \code{MMM} target markers.
#' 
#' @return A list containing the following named numerical vectors:
#' \itemize{
#' \item{\code{p_value}: } {
#' The p-values of the seven statistics.  The basic statics labeled \code{B}, \code{S}, and \code{T} 
#' are the squared burden, the SKAT, and the trait-based statistics, and their p-values are given by 
#' the \code{davies} function from the package \code{CompQuadForm}. THe p-values of the
#' optimized statistics, labeled \code{BS}, \code{BT}, \code{ST}, and \code{BST}, are gotten by
#' genome resampling. See \url{http://stanford.edu/~ggong/gatars/} for more details.
#' }
#' \item{\code{q}: }{
#' The observed seven statistics. This is \code{Q(alpha)} with alpha evaluated at its optimizing value.
#' }
#' \item{\code{x}: }{
#' The observed seven statistics transformed. This is \code{X_alpha = - log_10 (P_alpha)} 
#' where \code{P_alpha} is the minimal nominal p-value of \code{Q_alpha}.
#' }
#' } 
#' 
#' @template gatars_test_size_examples
#' @examples 
#' 
#' @export
gatars_test_size = function(phenotype, Psi, sampling_set, N_simulated_nulls, weights = NULL){
  g_target = sampling_set$g_target
  MMM = sampling_set$MMM
  sampling_set = sampling_set$sampling_set
  adaptive_conf_level = 0.99
  calculate_optimized = TRUE
  N_simulated_nulls_interval = N_simulated_nulls
  N_simulated_nulls_limit = N_simulated_nulls
  y_1 = yyy = phenotype$y
  y_2 = mu = phenotype$mu
  Phi = Phi_fn(Psi, y_1, y_2)
  theta_init = rep(pi/3, 2)
  www_num = if(!is.null(weights)){
    weights
  } else {
    rep(1, MMM)
  }
  www = www_num/sum(www_num) * MMM
  WWW = diag(www)
  www = t(t(www))
  ooo = uno_experimento_fn(
    adaptive_conf_level, calculate_optimized, g_target, MMM, 
    N_simulated_nulls_interval, N_simulated_nulls_limit, 
    Phi, sampling_set, theta_init, WWW, y_1, y_2)
  ooo
}