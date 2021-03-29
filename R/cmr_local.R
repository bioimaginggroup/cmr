#' Spline analysis of cardiovascular magnetic resonance imaging
#'
#' @param data 3d array of CMR signal 
#' @param mask 2d array of mask. Voxel with 0 or FALSE will be omitted from analysis
#' @param input input function
#' @param quantiles quantiles used for credible interval, default: c(0.25, 0.75)
#' @param cores number of cores to use in parallel computing
#'
#' @return list of mbf (point estimation) and ci (credible interval)
#' @export
#' @import splines Matrix parallel 
#' @importFrom stats rnorm rgamma median quantile
#' 
#' @examples
#' \donttest{
#' oldpar <- par(no.readonly = TRUE)
#'  library(cmR)
#'  data(cmrsim)
#'  local.mbf=local.ci=array(NA,c(30,30,3))
#'  for (i in 1:3){
#'   mask=array(NA,c(30,30))
#'   mask[cmrdata_sim[,,i,1]!=0]=1
#'   temp=cmr.local(cmrdata_sim[,,i,], mask, input_sim, cores=2)
#'   local.mbf[,,i]=t(as.matrix(temp$mbf))
#'   local.ci[,,i]=t(as.matrix(temp$ci))
#'  }
#'  par(mfrow=c(2,1))
#'  imageMBF(maxresp_sim, zlim=c(0,5))
#'  imageMBF(local.mbf, zlim=c(0,5))
#'  imageMBF(local.ci, zlim=c(0,0.8))
#' par(oldpar)
#' }
#' 
cmr.local<-function(data,mask,input,quantiles=c(.25,.75), cores=1)
{
  if(.Platform$OS.type == "windows") { cores=1 }
XX<-dim(data)[1]
YY<-dim(data)[2]
N<-sum(!is.na(mask))
T=dim(data)[3]

# compute B

zeit<-((1:T)-1)/60
# knots<-seq(-6.5,44,length=25)/60
# knots<-knots[-seq(13,25,by=2)]
# knots<-knots[-seq(10,18,by=2)]

#Eventuell bessere Wahl?
#knots=c(seq(-5,1,by=2),seq(3,13,by=1),seq(14,36,by=2))

# knots=seq(-2,33,length=24)
# knots<-knots[-seq(14,30,by=2)]
# 
knots<-seq(-1,T+1,length=T-2)
knots=knots/60#-1/60
k<-4
p<-length(knots)-k
B<-splines::splineDesign(knots,zeit,k,outer.ok=TRUE)

# compute A

A<-array(0,c(length(zeit),length(zeit)))
D<-array(0,c(dim(B)))
ni<-zeit
for (i in 1:length(zeit))
for (j in 1:length(zeit))
if (zeit[j]<=zeit[i])ni[i]=j

for (i in 1:length(zeit))
for (j in 1:ni[i])
{
A[i,j]<-input[1+ni[i]-j]
}
A<-A*mean(diff(zeit))

# compute D, DD

D<-A%*%B

i<-j<-x<-c()
for (ii in 1:T)
for (jj in 1:p)
if (D[ii,jj]!=0)
{
i=c(i,ii)
j=c(j,jj)
x=c(x,D[ii,jj])
}

D.sparse=Matrix::sparseMatrix(i, j, x=x, dims=c(T,p))
DD<-Matrix::t(D.sparse)%*%D.sparse

DC <- Ct <- c()
for (i in 1:XX)
for (j in 1:YY)
if (!is.na(mask[i,j]))
{
Ct=c(Ct,data[i,j,])
DC=c(DC,as.vector(Matrix::t(D.sparse)%*%data[i,j,]))
}

# prepare Q (Q is temporal prior)

Q.x<-c(1,2,3,1,2,3,4)
Q.y<-c(1,1,1,2,2,2,2)
for (i in 3:(p-2))
{
Q.x<-c(Q.x,i-2,i-1,i,i+1,i+2)
Q.y<-c(Q.y,i,i,i,i,i)
}
Q.x<-c(Q.x,p-3,p-2,p-1,p,p-2,p-1,p)
Q.y<-c(Q.y,p-1,p-1,p-1,p-1,p,p,p)


tauq2Q<-array(0,c(5*(p-2)+6,p-2))
tauq2Q[1:10,1]<-c(1,-2,1,-2,4,-2,0,1,-2,1)
for (i in 2:(p-2))
{
tauq2Q[0:10+(i-1)*5,i]<-c(1,-2,1,0,-2,4,-2,0,1,-2,1)
}
tauq2Q<-tauq2Q[-(p*5-8+c(0,4)),]

tauq<-rep(1,p-2)

Q.sparse=Matrix::sparseMatrix(Q.x,Q.y,x=as.vector(tauq2Q%*%tauq),dims=c(p,p))

coord<-c()

for (i in 2:XX)
for (j in 2:YY)
if (!is.na(mask[i,j]))
{
coord<-cbind(coord,c(i,j))
}

taueps=rep(1/10,N)

Q.klein<-matrix(c(1,-2,1,-2,4,-2,1,-2,1),nrow=3)

tauq.local<-rep(1,p-2)
Q.sparse=Matrix::sparseMatrix(Q.x,Q.y,x=as.vector(tauq2Q%*%tauq.local),dims=c(p,p))
taueps.local=1/10

system.time(temp<-parallel::mclapply(1:N,cmr.voxel,data,coord,Q.sparse, D.sparse, taueps.local, tauq.local, DD, T, p, B, Q.klein, Q.x, Q.y, tauq2Q, mc.cores=cores))
  
NRI<-300
response=array(NA,c(T,N,NRI))

for (voxel in 1:N)
{
response[,voxel,]=temp[[voxel]]
}

resp.max=apply(response,2:3,max)
q4<-function(x)stats::quantile(x,c(.25,.75),na.rm=TRUE)
resp.q4=apply(resp.max,1,q4)
resp.q4=resp.q4[2,]-resp.q4[1,]

resp.i<-Matrix::sparseMatrix(coord[1,],coord[2,],x=apply(resp.max,1,stats::median),dims=c(XX,YY))
resp.i4<-Matrix::sparseMatrix(coord[1,],coord[2,],x=resp.q4,dims=c(XX,YY))

resp.i[resp.i==0]<-NA
resp.i4[resp.i4==0]<-NA

return(list("mbf"=resp.i,"ci"=resp.i4))
}


cmr.voxel<-function(voxel,data,coord,Q.sparse, D.sparse, taueps.local, tauq.local, DD,T, p, B, Q.klein, Q.x, Q.y, tauq2Q){
  tauq.l.s<-taueps.l.s<-beta.l.s<-c()
  Ct <- data[coord[1,voxel],coord[2,voxel],]
  DC <- Matrix::t(D.sparse)%*%Ct
  
  for (iter in 1:1000)
  {
    # update beta
    
    L = Q.sparse + taueps.local*DD
    b = (taueps.local*DC)[,1]
    beta.local=rmvnormcanon(1,b,L)
    
    a=1+T/2
    bb=1e-3+0.5*sum((D.sparse%*%beta.local-Ct)^2)
    taueps.local=stats::rgamma(1,a,bb)
    
    for (i in 1:(p-2))
    {
      which=rep(FALSE,p)
      which[i:(i+2)]=TRUE
      aa=1+1/2
      bb=1e-3+0.5*beta.local[which]%*%Q.klein%*%beta.local[which]
      tauq.local[i]=rgamma(1,aa,bb[1,1])
    }
    
    Q.sparse=Matrix::sparseMatrix(Q.x,Q.y,x=as.vector(tauq2Q%*%tauq.local),dims=c(p,p))
    
    if (((iter%/%3)==iter/3)&iter>100)
    {
      tauq.l.s=cbind(tauq.l.s,tauq.local)
      taueps.l.s=c(taueps.l.s,taueps.local)
      beta.l.s=cbind(beta.l.s,beta.local)
    }
  }
  
  resp.local<-c()
  for (i in 1:300)
    resp.local<-cbind(resp.local,B%*%beta.l.s[,i])
  return(resp.local)
}

