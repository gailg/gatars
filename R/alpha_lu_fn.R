#' @export
alpha_lu_fn = function(theta){
  a1 = cos(theta[1])^2
  a2 = ( sin(theta[1]) * cos(theta[2]) )^2
  a3 = ( sin(theta[1]) * sin(theta[2]) )^2
  c(B = a1, S = a2, T = a3)
}