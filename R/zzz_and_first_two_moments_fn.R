#' @title Calculate zzz and its first two moments
#' 
#' @description  In the manuscript the first two moments are denoted \code{E_0(Z)}
#' and \code{Cov_0(Z)}, and are given by equation (8).  
#' In \code{gatars}, the first two moements are denoted
#' \code{mu_z} and \code{V_z}.  
#'
#' @param g_target A numerical matrix of dimension \code{NNN} by \code{MMM}
#' equal to what is referred to in the manuscript as \code{G}, the genotype matrix.
#' 
#' @param Phi A numerical matrix of dimension \code{2} by \code{2}.
#' \code{Phi_{k_1, k_2} = y_{k_1} Psi y_{k_2}}.  
#' This  matrix is a useful intermediate calculation for getting
#' \code{V_z}: \code{V_z = kronecker(Phi, W_VG_W)}. 
#' It is of dimension \code{2} by \code{2} because there are two entities 
#' \code{y_1} and \code{y_2}.
#' 
#' @param WWW A diagonal (numerical) matrix of dimension \code{MMM} by \code{MMM}
#' with the diagonals equal to the \code{weights}.  (The user will specify
#' \code{weights} in her call to \code{gatars_test_size}.)
#' 
#' @param y_1 A numerical vector of length \code{NNN} equal to what is referred
#' to in the manuscript as \code{y}, the vector of subjects' coded trait 
#' phenotypes.
#' 
#' @param y_2 A numerical vector of length \code{NNN} equal to what is referred
#' to in the manuscript as \code{mu}, the vector of user-specified phenotype
#' predictions.
#' 
#' @return A list containing 
#' \itemize{
#' \item{\code{zzz}: }{
#' A numerical vector of length \code{(2 * MMM)} equal to 
#' \code{c(WWW t(g_target) y_1, WWW t(g_target) y_2)}
#' and denoted \code{z} in equation (5) in the manuscript.
#' }
#' \item{\code{mu_z}: }{
#' A numerical vector of length \code{(2 * MMM)}  equal to \code{E_0(z)}
#' in equation (8) of the manuscript.
#' }
#' 
#' \item{\code{V_z}: }{
#' A numerical matrix of dimension \code{(2 * MMM)} by \code{(2 * MMM)}
#' equal to \code{Cov_0(z)} in equation (8) of the manuscript.
#' } 
#' }
#' 
#' @examples 
#' bim = gatars_example$bim
#' genotype = gatars_example$genotype
#' phenotype = gatars_example$phenotype
#' Psi = gatars_example$Psi
#' target_markers = gatars_example$target_markers[3:5]
#' 
#' g_target = genotype[, target_markers]
#' 
#' y_1 = yyy = phenotype$y
#' y_2 = mu = phenotype$mu
#' 
#' Phi = Phi_fn(Psi, y_1, y_2)
#' Phi
#' MMM = length(target_markers)
#' WWW =  diag(rep(1, MMM))
#' 
#' zzz_etc = zzz_and_first_two_moments_fn(g_target, Phi, WWW, y_1, y_2)
#' zzz_etc
#' str(zzz_etc)
#' 
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