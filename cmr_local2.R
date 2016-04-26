cmr.local.single<-function(Ct,D.sparse,tauq.local,Q.sparse,taueps.local,DD)
{
if (is.na(Ct[1]))return(rep(NA,300))
DC <- t(D.sparse)%*%Ct
tauq.l.s<-taueps.l.s<-beta.l.s<-c()

for (iter in 1:1000)
{
# update beta

L = Q.sparse + taueps.local*DD
b = taueps.local*DC
L=chol(as.matrix(L))
w=solve(t(L),b)
mu=solve(L,w)
z=rnorm(p)
v=solve(L,z)
beta.local=mu+v

a=1+T/2
bb=1e-3+0.5*sum((D.sparse%*%beta.local-Ct)^2)
taueps.local=rgamma(1,a,bb)

for (i in 1:(p-2))
{
which=rep(FALSE,p)
which[i:(i+2)]=TRUE
aa=1+1/2
bb=1e-3+0.5*beta.local[which]%*%Q.klein%*%beta.local[which]
tauq.local[i]=rgamma(1,aa,bb[1,1])
}

Q.sparse=sparseMatrix(Q.x,Q.y,x=as.vector(tauq2Q%*%tauq.local),dims=c(p,p))

if (((iter%/%3)==iter/3)&iter>100)
{
tauq.l.s=cbind(tauq.l.s,tauq.local)
taueps.l.s=c(taueps.l.s,taueps.local)
beta.l.s=cbind(beta.l.s,beta.local[,1])
}
}

resp.local<-c()
for (i in 1:300)
resp.local<-cbind(resp.local,B%*%beta.l.s[,i])

mbf.local=apply(resp.local,2,max)
return(mbf.local)
}


cmr.local<-function(data,mask,aif,quantiles=c(.25,.75))
{

require(Matrix)

XX<-dim(data)[1]
YY<-dim(data)[2]
N<-sum(!is.na(mask))

T=30

# compute B

library(splines)
zeit<-((1:T)-1)/60
#knots<-seq(-6.5,44,length=25)/60
#knots<-knots[-seq(13,25,by=2)]
#knots<-knots[-seq(10,18,by=2)]
knots=c(seq(-5,1,by=2),seq(3,13,by=1),seq(14,36,by=2))
knots=seq(-2,33,length=24)
knots<-knots[-seq(14,30,by=2)]
knots<-seq(-1,T+1,length=T-2)
knots=knots/60#-1/60
k<-4
p<-length(knots)-k
B<-splineDesign(knots,zeit,k,outer.ok=TRUE)

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
A[i,j]<-aif[1+ni[i]-j]
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

D.sparse=sparseMatrix(i, j, x=x, dims=c(T,p))
DD<-t(D.sparse)%*%D.sparse


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

Q.sparse=sparseMatrix(Q.x,Q.y,x=as.vector(tauq2Q%*%tauq),dims=c(p,p))




## lokale Schätzung

tauq.local<-rep(1,p-2)
Q.sparse=sparseMatrix(Q.x,Q.y,x=as.vector(tauq2Q%*%tauq.local),dims=c(p,p))
taueps.local=1/10
tauq.l.s<-taueps.l.s<-beta.l.s<-c()
print(dim(data))
for (i in 1:dim(data)[1])
for (j in 1:dim(data)[2])
if(is.na(mask[i,j]))data[i,j,1]=NA
mbf.local=apply(data,c(1,2),cmr.local.single,D.sparse,tauq.local,Q.sparse,taueps.local,DD)
print(dim(mbf.local))
qr<-function(x)
{
  if (sum(is.na(x)>0))
    return(NA)
  else
    return(quantile(x,.75)-quantile(x,.25))
}

  
return(list("mbf"=apply(mbf.local,2:3,median),"ci"=apply(mbf.local,2:3,qr)))
}



