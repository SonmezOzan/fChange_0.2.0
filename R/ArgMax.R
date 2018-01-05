
#' Exact quantile calcuations for Confidence inervals
#'
#' @name Argmax
#' @keywords internal



int_x = function(x, t1, t2){
  t2/(t1+t2) + 2*t1*(2*pi)^(-1/2)*sqrt(x)*exp(-t1^2*x/2)+
    ((t1*(t1+2*t2))/(t2*(t1+t2)))*exp(2*t2*(t1+t2)*x)*pnorm(-(t1+2*t2)*sqrt(x))-
    (2*t1^2*x+(t1^2+2*t2^2+2*t1*t2)/(t2*(t1+t2)))*pnorm(-t1*sqrt(x))
}

CdF = function(x, c1, c2, s1, s2){
  t1 = c1/s1
  t2 = c2*s1/s2^2
  t1p = c2/s2
  t2p = c1*s2/s1^2
  if(x>0 & !(is.na(x))){
    cdf = t2/(t1+t2) + int_x(x, t1p, t2p)
    if (is.na(cdf)){out = 1}else(out=cdf)
  } else{
    cdf = -int_x(-x, t1, t2) + t2/(t1+t2)
    if (is.na(cdf)){out = 0}else(out=cdf)
  }
  return(out)
}

newton <- function(f, tol=1e-7, x0=0, N=1000){
  h <- 1e-7
  i <- 1
  x1 <- x0
  p <- numeric(N)
  while(i <= N){
    df.dx <- (f(x0+h)-f(x0))/h
    x1 <- (x0 - (f(x0)/df.dx))
    p[i] <- x1
    i = i+1
    if(abs(x1-x0)<tol & !(is.na(abs(x1-x0)))) break
    x0=x1
  }
  return(p[(i-1)])
}

quant <- function(c.1, c.2, s.1, s.2, p){
  ftn = function(x){
    CdF(x, c1=c.1, c2=c.2, s1=s.1, s2=s.2) - p
  }
  newton(ftn)
}






