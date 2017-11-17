#' @export
zzz_and_first_two_moments_fn = function(g_target, Phi, WWW, y_1, y_2){
  z_1 = as.vector(WWW %*% t(g_target) %*% y_1)
  z_2 = as.vector(WWW %*% t(g_target) %*% y_2)
  zzz = c(z_1, z_2)
  
  e_g_target_1 = colMeans(g_target)
  e_g_target = matrix(rep(e_g_target_1, nrow(g_target)), nrow = nrow(g_target), byrow = TRUE)
  e_z_1 = matrix(as.vector(WWW %*% t(e_g_target) %*% y_1), ncol = 1)
  e_z_2 = matrix(as.vector(WWW %*% t(e_g_target) %*% y_2), ncol = 1)
  mu_z = rbind(e_z_1, e_z_2)
  
  V_G = cov(g_target)
  W_VG_W = WWW %*% V_G %*% WWW
  V_z = kronecker(Phi, W_VG_W)
  
  list(zzz = zzz,
       mu_z = mu_z,
       V_z = V_z)
}