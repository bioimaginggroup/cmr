rmvnormp <- function(n,m,Q)
{
require(Matrix)
Q=chol(Q)
w=solve(t(Q),m)
mu=solve(Q,w)
z=rnorm(length(m))
v=solve(Q,z)
beta=mu+v
return(beta)
}
