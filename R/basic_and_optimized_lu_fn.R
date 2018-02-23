#' @title Calculate the p-values of the basic statistics and the nominal p-values of the optimized statistics
#' 
#' @description For the basic statistics, calculate \code{Q} and its p-value.
#' For the optimized statistics, calculate the optimizing \code{Q(alpha)}
#' and its nominal p-value. The statistics \code{BS}, \code{BT},
#' and \code{ST} require optimizing on one variable so I use 
#' \code{optimize}.  The statistic \code{BST}, I use the lu
#' reparametrization given in the function \code{alpha_lu_fn}.
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
#' @param theta A vector of length \code{2} that is the initial value of the
#' reparametrization of alpha when I am finding
#' minimum p-value in the full triangle \code{(alpha_B, alpha_S, alpha_T)}
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
#' @return A list that contains the following objects
#' \itemize{
#' \item{\code{basic}: }{
#' A data frame containing the columns 
#' \code{B, S, T} (for \code{alpha_B, alpha_S, alpha_T}), \code{q} and \code{p-value}
#' and a row for each of the basic statistics \code{B}, \code{S}, and \code{T}.
#' }
#' \item{\code{optimized}: }{
#' A data frame containing the same columns as \code{basic} and the additional column
#' \code{x} and a row for each of the optimized statistics
#' \code{BS}, \code{BT}, \code{ST}, and \code{BST}.
#' }
#' \item{\code{p_value_basic}: }{
#' A named numerical vector containing the p-values of the three basic statistics.
#' }
#' \item{\code{q_basic}: }{
#' A named numerical vector containing the \code{Q_B}, \code{Q_S}, and \code{Q_T}.
#' }
#' \item{\code{q_optimized}: }{
#' A named numerical vector containing the optimized statistics \code{Q_BS},
#' \code{Q_BT}, \code{Q_ST}, and \code{Q_BST}.
#' }
#' \item{\code{theta}: }{
#' A vector of length \code{2} containing the value of \code{theta} that optimizing
#' \code{BST}.
#' }
#' \item{\code{x_observed}: }{
#' A named numerical vector containing the transformed values of the nominal p-values
#' of the optimized statistics \code{BS}, \code{BT}, \code{ST}, and \code{BST}. It is this
#' vector that will be compared with the simulated results to determine the true 
#' p-value.
#' }
#' }
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
#' theta_init = rep(pi/3, 2)
#' statistics = c("BS", "BT", "ST", "BST")
#' bo = basic_and_optimized_lu_fn(g_target, Phi, theta_init, WWW, y_1, y_2, statistics)
#' bo$optimized


#' @export
basic_and_optimized_lu_fn = function(g_target, Phi, theta, WWW, y_1, y_2, statistics){
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
    one_sided = if(jjj <= 2){
      TRUE
    } else {
      FALSE
    }
    davies_answer = davies_fn(zzz, mu_z, V_z, AAA, one_sided)
    cbind(alpha, davies_answer)
  }))
  # ------------------------- basic is a data.frame that contains 
  # ------------------------ the quadratic form q and the p-value 
  # ----------------------- of the three basic statistics B, S, T
  # ----------------------------------------------------- BS
  BS = if("BS" %in% statistics){
    p_value_BS_fn = function(alpha_B){
      AAA = AAA_fn(alpha_B, 1 - alpha_B, 0, MMM)
      davies_fn(zzz, mu_z, V_z, AAA, one_sided = TRUE)$p_value
    }
    ooo = optimize(p_value_BS_fn, interval = c(0, 1))
    alpha_B = ooo$minimum
    BS_0 = data.frame(B = alpha_B, S = 1 - alpha_B, T = 0, p_value = ooo$objective)
    possibles = rbind(BS_0, basic[c("B", "S"), names(basic) != "q"])
    BS = possibles[which.min(possibles$p_value), ]
    row.names(BS) = "BS"
    BS
  } else {
    NULL
  }
  # ----------------------------------------------------- BT
  BT = if("BT" %in% statistics){
    p_value_BT_fn = function(alpha_B){
      AAA = AAA_fn(alpha_B, 0, 1 - alpha_B, MMM)
      davies_fn(zzz, mu_z, V_z, AAA, one_sided = FALSE)$p_value
    }
    ooo = optimize(p_value_BT_fn, interval = c(0, 1))
    alpha_B = ooo$minimum
    BT_0 = c(B = alpha_B, S = 0, T = 1 - alpha_B, p_value = ooo$objective)
    possibles = rbind(BT_0, basic[c("B", "T"), names(basic) != "q"])
    BT =  possibles[which.min(possibles$p_value), ]
    row.names(BT) = "BT"
    BT
  } else {
    NULL
  }
  #------------------------------------------------------ ST
  ST = if("ST" %in% statistics){
    p_value_ST_fn = function(alpha_S){
      AAA = AAA_fn(0, alpha_S, 1 - alpha_S, MMM)
      davies_fn(zzz, mu_z, V_z, AAA, one_sided = FALSE)$p_value
    }
    ooo = optimize(p_value_ST_fn, interval = c(0, 1))
    alpha_S = ooo$min
    ST_0 = c(0, S = alpha_S, T = 1 - alpha_S, p_value = ooo$objective)
    possibles = rbind(ST_0, basic[c("S", "T"), names(basic) != "q"])
    ST =  possibles[which.min(possibles$p_value), ]
    row.names(ST) = "ST"
    ST
  } else {
    NULL
  }
  # ----------------------------------------------------- BST
  BST = if("BST" %in% statistics){
    p_value_BST_fn = function(theta){
      alpha = alpha_lu_fn(theta)
      AAA = AAA_fn(alpha[1], alpha[2], alpha[3], MMM)
      davies_fn(zzz, mu_z, V_z, AAA, one_sided = FALSE)$p_value
    }
    lu = optim(
      theta, p_value_BST_fn,
      method = "L-BFGS-B", lower = rep(0, 2), upper = rep(pi/2, 2))
    counts_lu = lu$counts
    theta = lu$par
    BST_0 = c(alpha_lu_fn(theta), p_value = lu$value)
    possibles = rbind(BST_0, basic[, names(basic) != "q"])
    BST =  possibles[which.min(possibles$p_value), ]
    row.names(BST) = "BST"
    BST
  } else {
    NULL
  }
  #----------------------------------------------------- optimized
  gather = list(BS, BT, ST, BST)
  optimized_0 = as.data.frame(do.call(rbind, gather[!sapply(gather, is.null)]))

  q = unlist(lapply(1:nrow(optimized_0), function(kkk){
    AAA = AAA_fn(optimized_0[kkk, "B"], optimized_0[kkk, "S"], optimized_0[kkk, "T"], MMM)
    as.vector(t(zzz) %*% AAA %*% zzz)
  }))
  x_observed = x = unlist(lapply(optimized_0$p_value, q_alice_fn))
  optimized = cbind(optimized_0, q, x)[, c("B", "S", "T", "q", "p_value", "x")]
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
  
  names(x_observed) = rownames(optimized)
  list(
    #AAA_basic = AAA_basic,
    basic = basic,
    optimized = optimized,
    p_value_basic = p_value_basic,
    q_basic = q_basic,
    q_optimized = q_optimized,
    theta = theta,
    x_observed = x_observed)
}