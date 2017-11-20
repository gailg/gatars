#` @title 


#' @export
basic_and_optimized_lu_fn = function(g_target, Phi, theta, WWW, y_1, y_2){
  # ---------------------------------------------------- zzz etc
  zzz_etc = zzz_and_first_two_moments_fn(g_target, Phi, WWW, y_1, y_2)
  zzz = zzz_etc$zzz
  mu_z = zzz_etc$mu_z
  V_z = zzz_etc$V_z
  # ------------------------------------------------------ basic
  basic_ones = data.frame(B = c(1, 0, 0),
                          S = c(0, 1, 0),
                          T = c(0, 0, 1))
  rownames(basic_ones) = c("B", "S", "T")
  MMM = ncol(g_target)
  AAA_basic = lapply(1:nrow(basic_ones), function(jjj){
    alpha = basic_ones[jjj, ]
    AAA = AAA_fn(alpha$B, alpha$S, alpha$T, MMM)
    list(alpha = alpha, AAA = AAA)
  })
  basic = do.call(rbind, lapply(1:nrow(basic_ones), function(jjj){
    alpha = basic_ones[jjj, ]
    AAA = AAA_fn(alpha$B, alpha$S, alpha$T, MMM)
    davies_answer = davies_fn(zzz, mu_z, V_z, AAA)
    cbind(alpha, davies_answer)
  }))
  # ------------------------- basic is a data.frame that contains 
  # ------------------------ the quadratic form q and the p-value 
  # ----------------------- of the three basic statistics B, S, T
  # ----------------------------------------------------- BS
  p_value_BS_fn = function(alpha_B){
    AAA = AAA_fn(alpha_B, 1 - alpha_B, 0, MMM)
    davies_fn(zzz, mu_z, V_z, AAA)$p_value
  }
  ooo = optimize(p_value_BS_fn, interval = c(0, 1))
  alpha_B = ooo$minimum
  BS_0 = data.frame(B = alpha_B, S = 1 - alpha_B, T = 0, p_value = ooo$objective)
  possibles = rbind(BS_0, basic[c("B", "S"), names(basic) != "q"])
  BS = possibles[which.min(possibles$p_value), ]
  # ----------------------------------------------------- BT
  p_value_BT_fn = function(alpha_B){
    AAA = AAA_fn(alpha_B, 0, 1 - alpha_B, MMM)
    davies_fn(zzz, mu_z, V_z, AAA)$p_value
  }
  ooo = optimize(p_value_BT_fn, interval = c(0, 1))
  alpha_B = ooo$minimum
  BT_0 = c(B = alpha_B, S = 0, T = 1 - alpha_B, p_value = ooo$objective)
  possibles = rbind(BT_0, basic[c("B", "T"), names(basic) != "q"])
  BT =  possibles[which.min(possibles$p_value), ]
  #------------------------------------------------------ ST
  p_value_ST_fn = function(alpha_S){
    AAA = AAA_fn(0, alpha_S, 1 - alpha_S, MMM)
    davies_fn(zzz, mu_z, V_z, AAA)$p_value
  }
  ooo = optimize(p_value_ST_fn, interval = c(0, 1))
  alpha_S = ooo$min
  ST_0 = c(0, S = alpha_S, T = 1 - alpha_S, p_value = ooo$objective)
  possibles = rbind(ST_0, basic[c("S", "T"), names(basic) != "q"])
  ST =  possibles[which.min(possibles$p_value), ]
  # ----------------------------------------------------- BST
  p_value_fn = function(theta){
    alpha = alpha_lu_fn(theta)
    AAA = AAA_fn(alpha[1], alpha[2], alpha[3], MMM)
    davies_fn(zzz, mu_z, V_z, AAA)$p_value
  }
  lu = optim(
    theta, p_value_fn,
    method = "L-BFGS-B", lower = rep(0, 3), upper = rep(pi/2, 2))
  counts_lu = lu$counts
  theta = lu$par
  BST_0 = c(alpha_lu_fn(theta), p_value = lu$value)
  possibles = rbind(BST_0, BS_0, BT_0, ST_0, basic[, names(basic) != "q"])
  BST =  possibles[which.min(possibles$p_value), ]
  #----------------------------------------------------- optimized
  optimized_0 = as.data.frame(rbind(BS, BT, ST, BST))
  q = unlist(lapply(1:nrow(optimized_0), function(kkk){
    AAA = AAA_fn(optimized_0[kkk, "B"], optimized_0[kkk, "S"], optimized_0[kkk, "T"], MMM)
    as.vector(t(zzz) %*% AAA %*% zzz)
  }))
  xxx = x = unlist(lapply(optimized_0$p_value, q_alice_fn))
  optimized = cbind(optimized_0, q, x)[, c("B", "S", "T", "q", "p_value", "x")]
  rownames(optimized) =  c("BS", "BT", "ST", "BST")
  optimized
  #------------------------------------------------------- output
  p_value_basic = basic$p_value
  names(p_value_basic) = rownames(basic)
  p_value_basic
  
  q_basic = basic$q
  names(q_basic) = rownames(basic)
  q_basic
  
  q_optimized = optimized$q
  names(q_optimized) = rownames(optimized)
  q_optimized
  
  names(xxx) = rownames(optimized)
  list(
    #AAA_basic = AAA_basic,
    basic = basic,
    optimized = optimized,
    p_value_basic = p_value_basic,
    q_basic = q_basic,
    q_optimized = q_optimized,
    theta = theta,
    xxx = xxx)
}