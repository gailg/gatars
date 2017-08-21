#' @export
q_alice_fn = function(p_value){
  p_value = if(p_value < 2.220446e-16) 2.220446e-16 else p_value
  - log10(p_value)
}