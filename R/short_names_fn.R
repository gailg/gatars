#' @export
short_names_fn = function(names){
  ifelse(names =="standard_burden", "B",
  ifelse(names =="standard_skat",   "S",
  ifelse(names == "altern_skat",    "T",
  ifelse(names == "standard_lee",   "BS",
  ifelse(names == "optim_skat",     "ST",
  ifelse(names == "odd",            "BT",
                                    "BST"))))))
}
