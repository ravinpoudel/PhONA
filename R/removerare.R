remove_rare <- function(x, minOTUcount=10, prevalence_percent= 5){
  nsample <- length(x)
  px <- sum(x > 0) ## prevalence
  pxp <- round(px/nsample * 100, 1)
  dx <- sum(x)
  ifelse ((pxp > prevalence_percent & dx > minOTUcount),TRUE, FALSE)
}
