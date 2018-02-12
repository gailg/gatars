#' @export
square_root_matrix_fn = function(AAA){
  eee = eigen(AAA, symmetric = TRUE)
  PPPT = t(eee$vectors)
  Lambda = eee$values
  square_root = crossprod(Lambda^(1/4) * PPPT)
  square_root_inverse = crossprod(Lambda^(-1/4) * PPPT)
  list(square_root = square_root,
       square_root_inverse = square_root_inverse,
       A_one_half = square_root,
       A_minus_one_half = square_root_inverse)
}