#' @export
Phi_fn = function(Psi, y_1, y_2){
  Phi_11 = t(y_1) %*% Psi %*% y_1
  Phi_12 = t(y_1) %*% Psi %*% y_2
  Phi_22 = t(y_2) %*% Psi %*% y_2
  Phi = matrix(c(Phi_11, Phi_12, Phi_12, Phi_22), nrow = 2)
  Phi
}