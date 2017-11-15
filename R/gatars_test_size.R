#' @title Genetic Association Tests for Arbitrarily Related Subjects
#' 
#' @description 
#' Test the association between a prespecified set of genetic markers, called target markers, 
#' and a binary or quantitative trait  
#' Calculate the p-values of a family of statistics to test the association between
#' a specified set of \code{MMM} markers and a binary or quantitative trait using subjects
#' with any genealogical relationship.
#' These statistics include three basic 
#' statistics--the squared burden statistic \eqn{Q_B}, the linear kernel (SKAT) statistic
#' \eqn{Q_S}, and a new trait-based kernel statistic \eqn{Q_T}--as well as an ensemble of four
#' statistics \eqn{Q_BS}, \eqn{Q_BT}, \eqn{Q_ST}, \eqn{Q_BST} that optimize linear combinations
#' of the three basic statistics.  See \url{http://stanford.edu/~ggong/gatars/} for more details.
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
#' correlation for the genotype at one marker between the \code{n1}-th and \code{n2}-th 
#' subjects. \code{Psi} can be estimated from known family pedigree structures and/or from
#' the subjects' genotypes at markers independent of those in the target set. The diagonal 
#' elements are all unity and two people known to be non-identical-twin full sibs and who
#' have parents known to be completely unrelated have correlation 0.5.
#' 
#' @param sampling_set An object of class \code{gatars_sampling_set} produced by 
#' \code{gatars_sampling_set}. (See \code{\link{gatars_sampling_set}}.)  
#' The object \code{sampling_set} is list containing the objects
#' \code{params_sampling_set}, \code{report}, and \code{sampling_set}.  The important object
#' is \code{sampling_set$sampling_set}, a list containing \code{MMM} matrices. 
#' To run genome resampling, \code{gatars_test_size}
#' creates "simulated genotype matrices" by sampling
#' and binding together one column from each of the \code{MMM} matrices.  
#' 
#' @param N_simulated_nulls An integer equal to the number of simulated genotype matrices to generate
#' when calculating the test sizes of the optimized statistics.


#' @template gatars_test_size_examples
#' 
#' @export
gatars_test_size = function(phenotype, Psi, sampling_set, N_simulated_nulls, weights = NULL){
  N_simulated_nulls = N_simulated_nulls
  params_sampling_set = sampling_set$params_sampling_set
  sampling_set = sampling_set$sampling_set
  alpha_uni_N_increments = 10
  params = params_fn(
    alpha_uni_N_increments, params_sampling_set, phenotype, Psi, sampling_set, N_simulated_nulls, weights)
  ooo = uno_experimento_fn(params, calculate_optimized = TRUE)
  ooo
}