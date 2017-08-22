#' @export
print.gatrs_sampling_set = function(x, ...){
  params_sampling_set = x$params_sampling_set
  report = x$report
  minimum_sampling_set_size = min(report$set_size)
  list(   
    sampling_set_report = report,
    minimum_sampling_set_size = minimum_sampling_set_size)
}
# S3method(print, gatrs_sampling_set)