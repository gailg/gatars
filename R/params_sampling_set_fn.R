#' @export
params_sampling_set_fn = function(
  bim,
  genotype,
  target_markers,
  exclusion_region,
  hotspot,
  epsilon_on_log_scale,
  columns = NULL
){
  g_target_0 = genotype[, target_markers]
  columns = if(is.null(columns)) 1:ncol(g_target_0) else columns
  MMM = length(columns)
  g_target = g_target_0[, columns]
  NNN = nrow(g_target)
  e_g_target_1 = colMeans(g_target)
  p_target = e_g_target_1/2
  e_g_target = matrix(rep(e_g_target_1, nrow(g_target)), nrow = nrow(g_target), byrow = TRUE)
  list(
    bim = bim,
    genotype = genotype,
    e_g_target = e_g_target,
    epsilon_on_log_scale = epsilon_on_log_scale,
    exclusion_region = exclusion_region,
    g_target = g_target,
    hotspot = hotspot,
    MMM = MMM,
    NNN = NNN,
    p_target = p_target,
    target_markers = target_markers)
}
