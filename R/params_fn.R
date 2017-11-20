#' @title The parameters that run the innards of \code{uno_experimento_fn}
#' 
#' @description This is my attempt to make the variables that run
#' \code{uno_experimento_fn} global but safer. Massage the inputs
#' and output the variables that will be neccessary for the innards
#' of \code{uno_experimento_fn}.
#' 
#' @param bim A \code{data.frame} containing the three columns
#' \code{chromosome}, \code{snp}, \code{bp}, and containing \code{LLL} rows 
#' corresponding to the \code{MMM} target markers and those used to build
#' the sampling sets. The \code{lll}-th row of \code{bim} summarizes the 
#' \code{lll}-th marker and corresponds to the \code{lll}-th column of 
#' \code{genotype} (described below).  The column labeled \code{chromosome} 
#' contains integers between \code{1} and \code{22} (other integers may be 
#' included, but only chromosomes \code{1} through \code{22} will be used),
#' the \code{snp} column contains character strings that name the marker,
#' and the \code{bp} column contains integers giving the markers' bp
#' positions.  Because the \code{gatars} dat set \code{hotspot} (described
#' below) is given in Build hg38/GRCh38, \code{bp} must also be expressed in
#' Build hg38/GRCh38.) The object \code{bim} could be the 
#' \code{.map} or \code{.bim} file from plink. 
#' @export
params_fn = function(
  params_sampling_set,
  phenotype,
  Psi,
  sampling_set,
  N_simulated_nulls,
  weights = NULL
){
  adaptive_conf_level = 0.99
  N_simulated_nulls_interval = N_simulated_nulls
  N_simulated_nulls_limit = N_simulated_nulls
  rho_uni = seq(0, 1, by = 0.1)
  test_size = 0.01
  y_1 = yyy = phenotype$y
  y_2 = mu = phenotype$mu
  Phi = Phi_fn(Psi, y_1, y_2)
  MMM = params_sampling_set$MMM
  theta_init = rep(pi/3, 2)
  www_num = if(!is.null(weights)){
    weights
  } else {
    rep(1, MMM)
  }
  www = www_num/sum(www_num) * MMM
  WWW = diag(www)
  www = t(t(www))
  list(
    adaptive_conf_level = adaptive_conf_level,
    mu = mu,
    g_target = params_sampling_set$g_target,
    MMM = params_sampling_set$MMM,
    NNN = params_sampling_set$NNN,
    N_simulated_nulls_interval = N_simulated_nulls_interval,
    N_simulated_nulls_limit = N_simulated_nulls_limit,
    Phi = Phi,
    sampling_set = sampling_set,
    test_size = test_size,
    theta_init = theta_init,
    WWW = WWW,
    www = www,
    yyy = yyy)
}
