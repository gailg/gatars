#' @title Call \code{davies}
#' 
#' @description Given the random vector (of functions of genotypes)
#' \code{Z} with mean vector \code{mu_z} and covariance matrix \code{V_z},
#' and given the matrix \code{AAA} calculate the p-value for the
#' statistic \code{q = t(zzz) AAA  zzz}. 
#' \code{davies_fn} calculates the spectral decomposition of \code{q}
#' to obtain the eigenvalues \code{lambda} and also the 
#' noncentrality parameter \code{delta} that are required by \code{davies}
#' of the \code{ComquadForm} package.
#' 
#' @param zzz A numerical vector of length \code{(2 * MMM)}, one of the objects returned by
#' \code{zzz_and_first_two_moments_fn}.
#' 
#' @param mu_z A numerical vector of length \code{(2 * MMM)}, one of the objects
#' returned by \code{zzz_and_first_two_moments_fn}.
#' 
#' @param V_z A numerical matrix of dimension \code{(2 * MMM)} by \code{(2 * MMM)},
#' one of the objects returned by \code{zzz_and_first_two_moments_fn}.
#' 
#' @param AAA A numerical matrix of dimension \code{(2 * MMM)} by \code{(2 * MMM)}, 
#' the object returned by \code{AAA_fn}
#' @param one_sided A logical equal to TRUE unless the statistic involves the 
#' trait statistic.
#' 
#' @return A one-row data.frame containing the columns \code{q} and \code{p_value}.
#' 
#' @examples 
#' bim = gatars_example$bim
#' genotype = gatars_example$genotype
#' phenotype = gatars_example$phenotype
#' Psi = gatars_example$Psi
#' target_markers = gatars_example$target_markers[3:5]
#' g_target = genotype[, target_markers]
#' MMM = ncol(g_target)
#' NNN = nrow(g_target)
#' e_g_target_1 = colMeans(g_target)
#' p_target = e_g_target_1/2
#' e_g_target = matrix(rep(e_g_target_1, nrow(g_target)), nrow = nrow(g_target), byrow = TRUE)
#' y_1 = yyy = phenotype$y
#' y_2 = mu = phenotype$mu
#' Phi = Phi_fn(Psi, y_1, y_2)
#' www_num = rep(1, MMM)
#' www = www_num/sum(www_num) * MMM
#' WWW = diag(www)
#' zzz_etc = zzz_and_first_two_moments_fn(g_target, Phi, WWW, y_1, y_2)
#' zzz = zzz_etc$zzz
#' mu_z = zzz_etc$mu_z
#' V_z = zzz_etc$V_z
#' AAA = AAA_fn(1, 0, 0, MMM)
#' davies_fn(zzz, mu_z, V_z, AAA, one_sided = TRUE)
#' 
#' @import CompQuadForm
#' 
#' @export
davies_fn = function(zzz, mu_z, V_z, AAA, one_sided){
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
  UUU = Re(t(eigen_A_tilde$vectors))
  lambda = Re(eigen_A_tilde$values)
  Lambda = diag(lambda)
  z_breve = UUU %*% z_tilde
  mu_z_breve = UUU %*% mu_z_tilde
  delta = mu_z_breve^2
  q_spectral_decomp = sum(lambda * z_breve^2)
  # ----------------------------------------------- p_value
  absolute_qqq = abs(qqq)
  right_tail = davies(absolute_qqq, lambda = lambda, delta = delta,
                      lim = 50000, acc = 0.00005)$Qq
  p_value = if(one_sided){
    right_tail
  } else {
    left_tail = davies(- absolute_qqq, lambda = lambda, delta = delta,
                       lim = 50000, acc = 0.00005)$Qq
    right_tail + 1 - left_tail
  }
  data.frame(q = qqq, p_value)
}