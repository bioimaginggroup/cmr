library(Matrix)
library(cmr)

data(sim)
resp[resp==0]<-NA
data[data==0]=NA

for (i in 1:dim(data)[1])
for (j in 1:dim(data)[2])
for (k in 1:3)
if (!is.na(data[i,j,k,1]))
data[i,j,k,]=data[i,j,k,]+rnorm(30,0,sqrt(30))

space.mbf=array(NA,c(30,30,3))
space.ci=array(NA,c(30,30,3))
local.mbf=array(NA,c(30,30,3))
local.ci=array(NA,c(30,30,3))

for (i in 1:3)
{
mask=array(NA,c(30,30))
mask[data[,,i,1]!=0]=1
temp=cmr.local(data[,,i,],mask,aif)
local.mbf[,,i]=t(as.matrix(temp$mbf))
local.ci[,,i]=t(as.matrix(temp$ci))
}

for (i in 1:3)
{
mask=array(NA,c(30,30))
mask[data[,,i,1]!=0]=1
temp=cmr.space(data[,,i,],mask,aif)
space.mbf[,,i]=t(as.matrix(temp$mbf))
space.ci[,,i]=t(as.matrix(temp$ci))
}

space.mbf[space.mbf==0]=NA
space.ci[space.ci==0]=NA
local.mbf[local.mbf==0]=NA
local.ci[local.ci==0]=NA
resp.resp[resp.resp==0]=NA

for (i in 1:3)
{
resp.resp[,,i]=t(resp.resp[,,i])
local.mbf[,,i]=t(local.mbf[,,i])
local.ci[,,i]=t(local.ci[,,i])
}


image.mbf(resp.resp,zlim=c(0,5))

image.mbf(space.mbf,zlim=c(0,5))

image.mbf(local.mbf,zlim=c(0,5))

image.mbf(space.ci,zlim=c(0,.8))

image.mbf(local.ci,zlim=c(0,.8))

image.mbf(resp.resp,zlim=c(0,5))

par(fin=10*c(dim(resp[,,1])/max(dim(resp[,,1]))))
par(mai=c(0,0,0,0))
par(cex=2)
image.plot(resp.resp[,,1],zlim=c(0,5),legend.only=TRUE,legend.width=3.5)

image.mbf(space.mbf,zlim=c(0,5))

image.mbf(local.mbf,zlim=c(0,5))

image.mbf(space.ci,zlim=c(0,.8))

image.mbf(local.ci,zlim=c(0,.8))

