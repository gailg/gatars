#' Calculate \code{A_alpha} in equation (7) of the manuscript
#' 
#' @param alpha_B A real number in the closed interval \code{[0, 1]}.
#' alpha_B + alpha_S + alpha_T should sum to 1.
#' @param alpha_S A real number in the closed interval \code{[0, 1]}.
#' @param alpha_T A real number in the closed interval \code{[0, 1]}.
#' @param MMM A positive integer equal to the length of \code{target_markers}
#' in \code{gatars_sampling_set}
#' 
#' @return A matrix of dimension \code{(2 * MMM)} by \code{(2 * MMM)}
#' 
#' @examples 
#' AAA_fn(1, 0, 0, 3)
#' AAA_fn(0, 1, 0, 3)
#' AAA_fn(0, 0, 1, 3)
#' @export
AAA_fn = function(alpha_B, alpha_S, alpha_T, MMM){
  III = diag(rep(1, MMM))
  JJJ = matrix(rep(1, MMM * MMM), nrow = MMM)
  A_1 = alpha_B * JJJ + (alpha_S + alpha_T) * III
  A_2 = -(alpha_B * JJJ + alpha_S * III)
  A_3 = alpha_B * JJJ + (alpha_S - alpha_T) * III
  AAA = rbind(cbind(A_1, A_2),
              cbind(A_2, A_3))
  AAA
}