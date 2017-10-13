#' @export
basic_and_optimized_fn = function(alpha_uni, g_target, Phi, WWW, y_1, y_2){
  # ---------------------------------------------------- zzz etc
  zzz_etc = zzz_etc_fn(g_target, Phi, WWW, y_1, y_2)
  zzz = zzz_etc$zzz
  mu_z = zzz_etc$mu_z
  V_z = zzz_etc$V_z
  # -------------------------------------------------- basic_AAA
  basic_ones = alpha_uni[apply(alpha_uni, 1, function(row){
    sum(row == 0) > 1
  }), ]
  MMM = ncol(g_target)
  AAA_basic = lapply(1:nrow(basic_ones), function(jjj){
    alpha = basic_ones[jjj, ]
    AAA = AAA_fn(alpha$B, alpha$S, alpha$T, MMM)
    list(alpha = alpha, AAA = AAA)
  })
  # ------------------------------------------------ df contains 
  # ----------------------B, S, T (the three components of alpha
  # ----------------------- as well as  q, and (nominal) p_value
  df = do.call(rbind, lapply(1:nrow(alpha_uni), function(jjj){
    alpha = alpha_uni[jjj, ]
    AAA = AAA_fn(alpha$B, alpha$S, alpha$T, MMM)
    davies_answer = davies_fn(zzz, mu_z, V_z, AAA)
    cbind(alpha, davies_answer)
  }))
  # ------------------------- basic is a data.frame that contains 
  # ------------------------ the quadratic form q and the p-value 
  # ----------------------- of the three basic statistics B, S, T
  basic = rbind(df[df$B == 1, ],
                df[df$S == 1, ],
                df[df$T == 1, ])
  rownames(basic) =  c("B", "S", "T")
  # --------------------- optimized is a data.frame that contains
  # ------------ the quadratic form q, the nominal p-value, and x
  # ----------------- of the optimized statistics BS, BT, ST, BST
  bs_0 = df[df$B + df$S == 1, ]
  bt_0 = df[df$B + df$T == 1, ]
  st_0 = df[df$S + df$T == 1, ]
  bs = bs_0[which.min(bs_0$p_value), ]
  bt = bt_0[which.min(bt_0$p_value), ]
  st = st_0[which.min(st_0$p_value), ]
  bst = df[which.min(df$p_value), ] 
  optimized = rbind(bs, bt, st, bst)
  xxx = unlist(lapply(optimized$p_value, q_alice_fn))
  optimized$x = xxx
  rownames(optimized) =  c("BS", "BT", "ST", "BST")
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
    AAA_basic = AAA_basic,
    basic = basic,
    optimized = optimized,
    p_value_basic = p_value_basic,
    q_basic = q_basic,
    q_optimized = q_optimized,
    xxx = xxx)
}