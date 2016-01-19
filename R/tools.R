
#' Function to create startvalues for nlsur models
#' @param model nlsur model
#' @param data the data frame used for evaluation
#' @param val value
#'
#' @export
getstartvals <- function(model, data, val) {
  # automatic creation of start values
  modelparameters <- unlist(lapply(model, all.vars))
  svals <- modelparameters[which(!modelparameters %in% names(data))]
  svals <- unique(svals[order(svals)])
  # svals
  strtvls <- rep(val, length(svals))
  names(strtvls) <- svals

  strtvls
}
