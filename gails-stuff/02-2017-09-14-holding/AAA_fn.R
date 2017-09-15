#' @export
AAA_fn = function(alpha, JJJ, III){
  A_1 = alpha[1] * JJJ + alpha[2] * III
  A_2 = alpha[3] * JJJ + alpha[4] * III
  AAA =  rbind(cbind(A_1 + A_2, -A_1),
               cbind(-A_1, A_1 - A_2))
  AAA
}