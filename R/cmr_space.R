#' Spatial spline analysis of cardiovascular magnetic resonance imaging
#'
#' @param data 3d array of CMR signal 
#' @param mask 2d array of mask. Voxel with 0 or FALSE will be omitted from analysis
#' @param input input function
#' @param quantiles quantiles used for credible interval, default: c(0.25, 0.75)
#'
#' @return list of mbf (point estimation) and ci (credible interval)
#' @export
#' @import splines Matrix 
#' @importFrom stats rnorm rgamma median quantile
#' @examples 
#' \donttest{
#' oldpar <- par(no.readonly = TRUE)
#'  library(cmR)
#'  data(cmrsim)
#'  mask=array(NA,c(30,30))
#'  space.mbf=space.ci=array(NA,c(30,30,3))
#'  for (i in 1:3){
#'   mask=array(NA,c(30,30))
#'   mask[cmrdata_sim[,,i,1]!=0]=1
#'   temp=cmr.space(cmrdata_sim[,,i,], mask, input_sim, cores=2)
#'   space.mbf[,,i]=t(as.matrix(temp$mbf))
#'   space.ci[,,i]=t(as.matrix(temp$ci))
#'   }
#'  par(mfrow=c(2,1))
#'  imageMBF(maxresp_sim, zlim=c(0,5))
#'  imageMBF(space.mbf, zlim=c(0,5))
#'  imageMBF(space.ci, zlim=c(0,0.8))
#'  par(oldpar)
#' }
#' 
cmr.space<-function(data,mask,input,quantiles=c(.25,.75))
{
  if(.Platform$OS.type == "windows") { cores=1 }
XX<-dim(data)[1]
YY<-dim(data)[2]
N<-sum(!is.na(mask))
T=dim(data)[3]

# compute B

zeit<-((1:T)-1)/60
#knots<-seq(-6.5,44,length=25)/60
#knots<-knots[-seq(13,25,by=2)]
#knots<-knots[-seq(10,18,by=2)]
# knots=c(seq(-5,1,by=2),seq(3,13,by=1),seq(14,36,by=2))
# knots=seq(-2,33,length=24)
# knots<-knots[-seq(14,30,by=2)]
knots=seq(-1,T+1,length=T-2)
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

rw=2

if (rw==1)
{
Q.x<-c(1,2,1)
Q.y<-c(1,1,2)
for (i in 2:(p-1))
{
Q.x<-c(Q.x,i,i+1,i)
Q.y<-c(Q.y,i,i,i+1)
}
Q.x<-c(Q.x,p)
Q.y<-c(Q.y,p)


tauq2Q<-array(0,c(length(Q.y),p-1))
for (i in 1:(p-1))
{
tauq2Q[1:4+(i-1)*3,i]<-c(1,-1,-1,1)
}


tauq<-rep(1,p-1)

Q.sparse=Matrix::sparseMatrix(Q.x,Q.y,x=as.vector(tauq2Q%*%tauq),dims=c(p,p))

}

if (rw==2)
{
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
}

#prepare R

coord<-c()

for (i in 1:XX)
for (j in 1:YY)
if (!is.na(mask[i,j]))
{
coord<-cbind(coord,c(i,j))
}

nei<-c()
for (i in 1:(N-1))
for (j in (i+1):N)
if(sum((coord[,i]-coord[,j])^2)<2)
nei<-cbind(nei,c(i,j))

NEI <- dim(nei)[2]

R.x <- c(1:N,nei[1,],nei[2,])
R.y <- c(1:N,nei[2,],nei[1,])

taur2R<-array(0,c(2*NEI+N,NEI))
for (j in 1:NEI)
{
taur2R[nei[1,j],j]<-taur2R[nei[1,j],j]+1
taur2R[nei[2,j],j]<-taur2R[nei[2,j],j]+1
taur2R[N+j,j]<- -1
taur2R[N+NEI+j,j]<- -1
}

taur<-rep(1,NEI)
R.sparse=Matrix::sparseMatrix(R.x,R.y,x=as.vector(taur2R%*%taur),dims=c(N,N))

QR=Matrix::kronecker(R.sparse,Q.sparse)
taueps=rep(1/10,N)

IDD = Matrix::kronecker(Matrix::Diagonal(N,taueps),DD)
IDC = rep(taueps,each=p)*DC

if (rw==2)Q.klein<-matrix(c(1,-2,1,-2,4,-2,1,-2,1),nrow=3)
if (rw==1)Q.klein<-matrix(c(1,-1,-1,1),nrow=2)
R.klein<-matrix(c(1,-1,-1,1),nrow=2)

tauq.s<-taueps.s<-beta.s<-taur.s<-c()

do.taueps<-function(i, D.sparse,beta,p,Ct,T)
{
  bb=1e-3+sum(0.5*(D.sparse%*%beta[(1:p)+(i-1)*p]-Ct[(1:T)+(i-1)*T])^2)
  return(stats::rgamma(1,1+T/2,bb))
}

for (iter in 1:200)
{

# update beta
QR=Matrix::kronecker(Matrix::Diagonal(N),Q.sparse)+Matrix::kronecker(R.sparse,Matrix::Diagonal(p))

L = QR + IDD
b = IDC
beta=rmvnormcanon(1,b,L)

taueps=unlist(parallel::mclapply(1:N, do.taueps,D.sparse,beta,p,Ct,T))

IDD = Matrix::kronecker(Diagonal(N,taueps),DD)
IDC = rep(taueps,each=p)*DC

#QkR<-kronecker(Q.klein,R.sparse)
#QkR<-kronecker(R.sparse,Q.klein)
for (i in sample(1:(p-rw)))
{
#which=rep(FALSE,p)
#which[i:(i+rw)]=TRUE
#which=rep(which,N)
beta1=beta[(0:(N-1))*p+i]
beta2=beta[(0:(N-1))*p+i+1]
if(rw==2)beta3=beta[(0:(N-1))*p+i+2]
aa=1+(N-1)/2
beta1=beta1-beta2
if(rw==2)beta1=beta1-beta2+beta3
#bb=1+0.5*beta[which]%*%QkR%*%beta[which]
bb=1+0.5*beta1%*%R.sparse%*%beta1
tauq[i]=rgamma(1,aa,bb[1,1])
}

Q.sparse=Matrix::sparseMatrix(Q.x,Q.y,x=as.vector(tauq2Q%*%tauq),dims=c(p,p))

if(iter>2)
{
#QRk<-kronecker(Q.sparse,R.klein)
#QRk<-kronecker(R.klein,Q.sparse)
for (i in sample(1:NEI))
{
#which=rep(FALSE,N)
#which[nei[,i]]=TRUE
#which=rep(which,each=p)
beta1=beta[(nei[1,i]-1)*p+(1:p)]
beta2=beta[(nei[2,i]-1)*p+(1:p)]
aa=1+(p-rw)/2
#bb=1+0.5*beta[which]%*%QRk%*%beta[which]
#bb=1 for under stress!
#bb=2 for at rest!
bb=1+0.5*(beta1-beta2)%*%Q.sparse%*%(beta1-beta2)
taur[i]=rgamma(1,aa,bb[1,1])
}

R.sparse=Matrix::sparseMatrix(R.x,R.y,x=as.vector(taur2R%*%taur),dims=c(N,N))
}

if (iter>10)
{
tauq.s=cbind(tauq.s,tauq)
taur.s=cbind(taur.s,taur)
taueps.s=c(taueps.s,taueps)
beta.s=cbind(beta.s,beta)
}
}


NRI=dim(beta.s)[2]
response<-array(NA,c(T,N,NRI))
for (i in 1:NRI)
for (j in 1:N)
response[,j,i]<-B%*%beta.s[(1:p)+(j-1)*p,i]

resp.max=apply(response,2:3,max)
q4<-function(x)stats::quantile(x,c(.25,.75),na.rm=TRUE)
resp.q4=apply(resp.max,1,q4)
resp.q4=resp.q4[2,]-resp.q4[1,]

resp.i<-Matrix::sparseMatrix(coord[1,],coord[2,],x=apply(resp.max,1,stats::median),dims=c(XX,YY))
resp.i4<-Matrix::sparseMatrix(coord[1,],coord[2,],x=resp.q4,dims=c(XX,YY))
resp.i[resp.i==0]<-NA
resp.i4[resp.i4==0]<-NA

return(list("mbf"=resp.i,"ci"=resp.i4))
#,"beta.s"=beta.s,"coord"=coord,"D"=D,"B"=B,"taur"=taur.s))
}
