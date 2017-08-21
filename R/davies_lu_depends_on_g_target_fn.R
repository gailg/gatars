#' @export
davies_lu_depends_on_g_target_fn = function(
  g_target,
  MMM,
  rho_uni,
  Phi,
  theta,
  WWW,
  www,
  y_1,
  y_2
 ){
  JJJ = matrix(rep(1, MMM * MMM), nrow = MMM)
  III = diag(rep(1, MMM))
  zero_M = matrix(rep(0, MMM * MMM), nrow = MMM)
  zero_M_2M = matrix(rep(0, MMM * 2 * MMM), nrow = MMM)
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
  # ----------------------------------------------- straight, fancy_grid
  alpha_uni_zero = expand.grid(a1 = rho_uni, a2 = rho_uni, a3 = rho_uni)
  alpha_uni = alpha_uni_zero[apply(alpha_uni_zero, 1, sum) <= 1, ]
  answer_0 = do.call(rbind, lapply(1:nrow(alpha_uni), function(jjj){
    alpha = alpha_uni[jjj, ]
    a4 = (1 - sum(as.vector(alpha)))
    A_1 = alpha[, 1] * JJJ + alpha[,2] * III
    A_2 = alpha[, 3] * JJJ + a4 * III
    AAA =  rbind(cbind(A_1 + A_2, -A_1),
                 cbind(-A_1, A_1 - A_2))
    davies_answer = davies_fn(zzz, mu_z, V_z, AAA)
    cbind(alpha, a4, davies_answer)
  }))
  a1 = answer_0$a1
  a2 = answer_0$a2
  a3 = answer_0$a3
  a4 = answer_0$a4
  standard_lee_0 = answer_0[a1 + a2 == 1, ]
  optim_skat_0   = answer_0[a2 + a4 == 1, ]
  odd_0          = answer_0[a1 + a4 == 1, ]
  #--------------------------------------------------------------------- pause for straight
  straight = rbind(standard_burden = answer_0[a1 == 1, ],
                   standard_skat   = answer_0[a2 == 1, ],
                   altern_skat     = answer_0[a4 == 1, ])
  q_straight = straight$q
  names(q_straight) = rownames(straight)
  p_value = data.frame(t(straight$p_value))
  names(p_value) = rownames(straight)
  #--------------------------------------------------------------------- pause for straight
  fancy_grid = rbind(
    standard_lee    = standard_lee_0[which.min(standard_lee_0$p_value), ],
    optim_skat      = optim_skat_0[which.min(optim_skat_0$p_value), ],
    odd             = odd_0[which.min(odd_0$p_value), ],
    optim           = answer_0[which.min(answer_0$p_value), ])
  fancy_grid$q = unlist(lapply(fancy_grid$p_value, q_alice_fn))
  # ---------------------------------------------------------------------- lu
  davies_alpha_fn = function(alpha){
    AAA = AAA_fn(alpha, JJJ, III)
    aaa = data.frame(t(alpha))
    names(aaa) = paste0("a", 1:4)
    cbind(aaa, davies_fn(zzz, mu_z, V_z, AAA))
  }
  # --------------------------------------------------- standard_lee
  p_value_a1_a2_fn = function(a1){
    alpha = c(a1, 1 - a1, 0, 0)
    AAA = AAA_fn(alpha, JJJ, III)
    davies_fn(zzz, mu_z, V_z, AAA)$p_value
  }
  ooo = optimize(p_value_a1_a2_fn, interval = c(0, 1))
  a1 = ooo$minimum
  q = q_alice_fn(ooo$objective)
  standard_lee = data.frame(a1, a2 = 1 - a1, a3 = 0, a4 = 0, q, p_value = ooo$objective)
  # ----------------------------------------------------- optim_skat
  p_value_a2_a4_fn = function(a2){
    alpha = c(0, a2, 0, 1 - a2)
    AAA = AAA_fn(alpha, JJJ, III)
    davies_fn(zzz, mu_z, V_z, AAA)$p_value
  }
  ooo = optimize(p_value_a2_a4_fn, interval = c(0, 1))
  a2 = ooo$minimum
  q = q_alice_fn(ooo$objective)
  optim_skat = data.frame(a1 = 0, a2, a3 = 0, a4 = 1 - a2, q, p_value = ooo$objective)
  # --------------------------------------------------- odd
  p_value_a1_a4_fn = function(a1){
    alpha = c(a1, 0, 0, 1 - a1)
    AAA = AAA_fn(alpha, JJJ, III)
    davies_fn(zzz, mu_z, V_z, AAA)$p_value
  }
  ooo = optimize(p_value_a1_a4_fn, interval = c(0, 1))
  a1 = ooo$minimum
  q = q_alice_fn(ooo$objective)
  odd = data.frame(a1, a2 = 0, a3 = 0, a4 = 1 - a1, q, p_value = ooo$objective)
  # ---------------------------------------------------------- lu
  p_value_fn = function(theta){
    alpha = alpha_lu_with_theta_3_zeroed_fn(theta)
    AAA = AAA_fn(alpha, JJJ, III)
    davies_fn(zzz, mu_z, V_z, AAA)$p_value
  }
  lu = optim(
    theta, p_value_fn,
    method = "L-BFGS-B", lower = rep(0, 3), upper = rep(pi/2, 2))
  counts_lu = lu$counts
  lu_theta = lu$par
  fancy_lu = rbind(standard_lee = standard_lee,
                   optim_skat = optim_skat,
                   odd = odd,
                   optim = data.frame(t(alpha_lu_with_theta_3_zeroed_fn(lu_theta)), q = q_alice_fn(lu$value), p_value = lu$value))
  #---------------------------------------------------------------------------------- fancy_best
  fancy_best = do.call(rbind, lapply(rownames(fancy_grid), function(row_name){
    larry = rbind(grid = fancy_grid[row_name, ],
                  lu = fancy_lu[row_name, ])
    winner = which.min(larry$p_value)
    answer = larry[winner, ]
    answer$winner = rownames(larry)[winner]
    answer$difference = - diff(larry$p_value)
    answer
  }))
  rownames(fancy_best) =rownames(fancy_grid)
  q_fancy = fancy_best$q
  names(q_fancy) = rownames(fancy_best)

  qqq = data.frame(t(q_straight), t(q_fancy))

  list(fancy_best = fancy_best,
       fancy_grid = fancy_grid,
       fancy_lu = fancy_lu,
       counts_lu = counts_lu,
       lu_theta = lu_theta,
       p_value = p_value,
       qqq = qqq,
       straight = straight)
}
