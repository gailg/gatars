#' @title Test the association between genetic markers and disease. 
#' 
#' @description (Genetic Association Tests for Arbitrarily Related Subjects) 
#' Calculate the p-values of a family of statistics that test the association between
#' a specified set of \code{MMM} makrders and a binary trait.  These statistics include three 
#' statistics--the squared burden statistic \eqn{Q_B}, the linear kernel (SKAT) statistic
#' \eqn{Q_S}, a new case-based kernel statistic \eqn{Q_C}--as well as an ensemble of four
#' statistics \eqn{Q_BS}, \eqn{Q_BC}, \eqn{Q_SC}, \eqn{Q_BSC} that optimize linear combinations
#' of the three basic statistics.  See \url{http://stanford.edu/~ggong/gatars/} for more details.
#' 
#' @param fam A data.frame containing \code{NNN} rows and at least the two columns
#' \code{y} and \code{e_y}, where \code{y[nnn]} is an indicator for the \code{nnn}-th person
#' having disease and \code{e_y[nnn]} is is a trait probability, the predicted value of 
#' \code{y[nnn]} based on nongenetic covariates and possibly principal components of ancestry. 
#' After we saw that principal conponents of ancestry were insignificant, we performed a 
#' logistic regression of \code{y} on (1) age, (2) membership in family or case/control data, 
#' and (3) their interaction, and used \code{e_y} to be the fitted values of the logistic 
#' regression.
#' 
#' @param Psi A matrix with \code{NNN} rows and \code{NNN} columns.  \code{Psi[n1, n2]} is the 
#' correlation for the genotype at one marker between the \code{n1}-th and \code{n2}-th 
#' subjects. \code{Psi} can be estimated from known family pedigree structures and/or from
#' the subjects' genotypes at markers independent of those in the target set. The diagonal 
#' elements are all unity and two people known to be non-identical-twin full sibs and who
#' have parents known to be completely unrelated have correlation 0.5.
#' 
#' @param sampling_set An object of class \code{gatars_sampling_set} produced by 
#' \code{gatars_sampling_set}. (See \code{\link{gatars_sampling_set}}.)  The object \code{sampling_set} is list containing the objects
#' \code{params_sampling_set}, \code{report}, and \code{sampling_set}.  The important object
#' is \code{sampling_set$sampling_set}, a list containing \code{MMM} matrices. When it is time 
#' for genome resampling, \code{gatars} creates "simulated genotype matrices" by sampling
#' and cbinding together one column from each of the \code{MMM} matrices.  
#' 
#' @param N_sim_reps An integer equal to the number of "simulated genotype matrices to generate
#' when calculating the test sizes of the optimal statistics.


#' @template gatars_examples
#' 
#' @export
gatars = function(fam, Psi, sampling_set, N_sim_reps, weights = NULL){
  params_sampling_set = sampling_set$params_sampling_set
  sampling_set = sampling_set$sampling_set
  params = params_fn(params_sampling_set, fam, Psi, sampling_set, N_sim_reps, weights)
  ooo = one_experiment_fn(params, calculate_fancy = TRUE)
  ooo$p_value
}