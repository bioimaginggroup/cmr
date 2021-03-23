#' Draw random vectors from multivariate Gaussian in canonical form
#'
#' @param n Number of draws
#' @param b b parameter
#' @param P Precision matrix
#'
#' @return matrix with n columns, vector if n=1
#'
#' @export
#' @examples
#' P<-matrix(c(1,.5,.5,1),ncol=2)
#' b=c(2,0)
#' # expected value and covariance matrix
#' Sigma = solve(P)
#' mu = b%*%Sigma
#' # sample
#' x<-rmvnormcanon(1000,b,P)
#' mu.hat=apply(x,1,mean)
#' print(mu.hat-mu)
#' Sigma.hat=var(t(x))
#' print(Sigma.hat-Sigma)
#' 
rmvnormcanon <- function(n,b, P){
  L <- chol(P) 
  u <- solve(t(L), b) 
  v <- solve(L, u) 
  if (n==1)
  {
    z <- rnorm(n = length(b), mean = 0, sd = 1)
    m <- solve(L, z) 
    Gamma <- v + m
  }
  if (n>1)
  {
    z<-parallel::mclapply(1:n,function(i,n)rnorm(n = n, mean = 0, sd = 1),n=length(b))
    m<-parallel::mclapply(z,function(z,L)return(solve(L,z)),L)
    Gamma<-parallel::mclapply(m,function(m,v)return(v+m),v)
    Gamma<-matrix(unlist(Gamma),ncol=n)
  }
  
  return(Gamma)
}
