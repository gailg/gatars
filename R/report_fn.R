#' @export
report_fn = function(chosen_sim){
  do.call(rbind, lapply(chosen_sim, function(this){
    ppp = apply(this, 2, mean)/2
    min = min(ppp)
    max = max(ppp)
    N = length(ppp)
    data.frame(min, max, N)
  }))
}
#  "2015-07-16 09:22:34 PDT"