#' @import CompQuadForm
#' @export
davies_fn = function(zzz, mu_z, V_z, AAA){
  qqq = as.vector(t(zzz) %*% AAA %*% zzz)
  sss = square_root_matrix_fn(V_z)
  V_z_one_half = sss$A_one_half
  V_z_minus_one_half = sss$A_minus_one_half
  z_tilde = V_z_minus_one_half %*% zzz
  mu_z_tilde =  V_z_minus_one_half %*% mu_z
  V_z_tilde = V_z_minus_one_half %*% V_z %*% V_z_minus_one_half
  A_tilde = V_z_one_half %*% AAA %*% V_z_one_half
  # ------------------------------------------------- breve, lambda, delta
  eigen_A_tilde = eigen(A_tilde)
  UUU = t(eigen_A_tilde$vectors)
  lambda = eigen_A_tilde$values
  Lambda = diag(lambda)
  z_breve = UUU %*% z_tilde
  mu_z_breve = UUU %*% mu_z_tilde
  delta = mu_z_breve^2
  q_spectral_decomp = sum(lambda * z_breve^2)
  # ----------------------------------------------- p_value
  p_value = davies(q = qqq, lambda = lambda, delta = delta,
                   lim = 50000, acc = 0.00005)$Qq
  data.frame(q = qqq, p_value)
}