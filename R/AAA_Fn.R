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