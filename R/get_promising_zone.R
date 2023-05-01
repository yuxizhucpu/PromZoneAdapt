#' Get the promising zone and sample size adaptation
#'
#' This function computes the promising zone and 
#' provide newly adapted sample size for second stage 
#' if the interim result is located in the promising zone
#' 
#' @param beta: type II error
#' @param cb: critical boundary at the second stage
#' @param x2: information vector
#' @param s.size: sample size designed at the beginning
#' @param s.max: maximum sample size based on the budget
#' @return a list of promising zone and a table for sample size adaptation based on interim test statistics



#' @export

get_promising_zone <- function(beta,cb,x2,s.size, s.max){
  
  d2new=rep(0, length(seq(0.01,3,0.01)))
  cpnew=rep(0, length(seq(0.01,3,0.01)))
  b1adjust=rep(0, length(seq(0.01,3,0.01)))
  
  niter=1
  for (zt in seq(0.01,3,0.01)) { 
    
    res=rep(0,length((s.size*(1-x2[1])):(s.max-s.size*x2[1])))
    inter1=1
    
    for (dnew in (s.size*(1-x2[1])):(s.max-s.size*x2[1])) {
      
      res[inter1]=qnorm(beta)-(cb*sqrt(s.size*x2[1]+dnew)-zt*sqrt(s.size*x2[1])-zt*(dnew/sqrt(s.size*x2[1])))/(sqrt(dnew))
      inter1=inter1+1
    }
    
    if (sum(res>=0)==0) {
      d2new[niter]=s.max-s.size * x2[1]
    } else if (sum(res<=0)==0) {
      d2new[niter]=s.size * x2[1]
    } else {
      d2new[niter]= min(((s.size*x2[1]):s.max)[which(res>0)])
    }
    
    cpnew[niter]=1-pnorm((cb*sqrt(s.size*x2[1]+d2new[niter])-zt*sqrt(s.size*x2[1])-zt*(d2new[niter]/sqrt(s.size*x2[1])))/(sqrt(d2new[niter])))
    
    b1adjust[niter]=((cb-zt*sqrt(x2[1]))/sqrt((1-x2[1])))*(sqrt((d2new[niter]/(d2new[niter]+s.size*x2[1]))))+ zt*sqrt(1-(d2new[niter]/(d2new[niter]+s.size*x2[1])))
    
    
    niter=niter+1
    
  }
  
  newd2_power=cbind(d2new, cpnew,b1adjust,cb, zi=seq(0.01,3,0.01))
  
  idx=intersect(which(cpnew>=0.9), which((b1adjust-cb)<0))
  
  promisingspace=c(seq(0.01,3,0.01)[idx][1],seq(0.01,3,0.01)[idx][length(seq(0.01,3,0.01)[idx])])
  conditional_power_range=round(c(get_conditional_power(promisingspace[1], cb, x2),get_conditional_power(promisingspace[2], cb, x2)),2)
  
  restore=list(promisingarea=promisingspace,conditionalpower= conditional_power_range, d2_power= newd2_power[idx,] )
  
  return( restore )
}