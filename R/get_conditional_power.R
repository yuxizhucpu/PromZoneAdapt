#' Get the conditional power based on the interim result
#'
#' This function computes the conditional power 
#' based on the interim results
#' 
#' @param zt: test statistic at the first stage
#' @param cb: critical boundary at the second stage
#' @param x2: information vector
#' @return conditional power


#' @export

get_conditional_power <- function(zt, cb, x2){
  
  conditional_power <- 1-pnorm((cb-zt*sqrt(x2[1])-(zt*(1-x2[1]))/sqrt(x2[1]))/(sqrt(1-x2[1])))
  return(conditional_power)
  
}