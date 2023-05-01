hermite <- function (points, z) {
  p1 <- 1/pi^0.4
  p2 <- 0
  for (j in 1:points) {
    p3 <- p2
    p2 <- p1
    p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
  }
  pp <- sqrt(2 * points) * p2
  c(p1, pp)
}

mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")
  
  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }
  
  ## rotate, scale, translate points
  eig <- eigen(sigma) 
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
}


#' Generate the critical boundaries of two-stage design
#'
#' This function creates the critical boundaries of two-stage design
#' based on the alpha spending function
#' user can specify the alpha spending function by themself
#' 
#' @param x1: type 1 error (alpha)
#' @param x2: information vector
#' @param x3: alpha spending function. Pb and OBFb
#' @return critical boundaries of two-stage design



#' @export
get_two_stage_boundary <- function(x1, x2,x3){
  
  if (x3=="Pb") {
    asf=function(information) x1*log(1+(exp(1)-1)*information) 
  } else if (x3=="OBFb"){
    asf=function(information) 2-2*pnorm(qnorm(1-(x1/2))/sqrt(information))
  }
  
  
  a1=asf(x2[1])
  
  a2=asf(x2[2])-asf(x2[1])
  
  b1=round(qnorm(1-(a1/2)),3)
  
  sig <- matrix(c(1,sqrt(x2[1]),sqrt(x2[1]),1),2,2)
  pts <- mgauss.hermite(120, mu=c(0,0), sigma=sig, prune=0.1)
  
  gfun <- function(x,b2) prod( abs(x[1]) <= b1 , abs(x[2])> b2)
  
  res_store=rep(0,400)
  
  j=1
  for (i in seq(0.01,4,0.01)) {
    res_store[j]=sum(apply(pts$points, 1, gfun,b2=i) * pts$weights)
    j=j+1
  }
  
  b2= max(round( seq(0.01,4,0.01)[ which(abs(res_store-a2)<0.001)], 3))
  
  return(c(b1,b2))
}
