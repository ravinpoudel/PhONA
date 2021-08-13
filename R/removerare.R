#' Remove rare observation based on prevalence percentage and/or minimum count
#'
#' This function removes rare feature based on revalence percentage and/or minimum count
#'
#' @param x Vector of observations for a feature or a column of a feature from a feature matrix.
#' @return Boolean. If a feature is within the defined criteria, return True, else False
#' @export
#'
remove_rare <- function(x, minOTUcount=10, prevalence_percent= 5){
  nsample <- length(x)
  px <- sum(x > 0) ## prevalence
  pxp <- round(px/nsample * 100, 1)
  dx <- sum(x)
  ifelse ((pxp > prevalence_percent & dx > minOTUcount),TRUE, FALSE)
}
