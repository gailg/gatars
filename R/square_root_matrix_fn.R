#' @export
square_root_matrix_fn = function(AAA){
  eee = eigen(AAA)
  PPP = eee$vectors
  Lambda = diag(eee$values)
  square_root = PPP %*% Lambda^(1/2) %*% t(PPP)
  square_root_inverse = PPP %*% diag(eee$values^(-1/2)) %*% t(PPP)
  list(square_root = square_root,
       square_root_inverse = square_root_inverse,
       A_one_half = square_root,
       A_minus_one_half = square_root_inverse)
}