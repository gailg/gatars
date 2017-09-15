#' @export
alpha_lu_with_theta_3_zeroed_fn = function(theta){
  a1 = cos(theta[1])^2
  a2 = ( sin(theta[1]) * cos(theta[2]) )^2
  a3 = 0
  a4 = ( sin(theta[1]) * sin(theta[2]) )^2
  c(a1 = a1, a2 = a2, a3 = a3, a4 = a4)
}