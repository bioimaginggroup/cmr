image.mbf<-function(img,zlim=NULL)
{
require(fields)
Z=dim(img)[3]
if (is.null(zlim))
zlim=c(0,max(img,na.rm=TRUE))

img[img>zlim[2]]=zlim[2]
img[img<zlim[1]]=zlim[1]

sum.na<-function(x)return(sum(x,na.rm=TRUE))
yrange=range(which(apply(img[,,],2,sum.na)!=0))
xrange=range(which(apply(img[,,1],1,sum.na)!=0))
drange=max(diff(xrange)-diff(yrange),0)/2
yrange=yrange+floor(c(-1,1)*drange)

xrange=range(which(apply(img[,,1],1,sum.na)!=0))
fullimg=img[xrange[1]:xrange[2],yrange[1]:yrange[2],1]
for (i in 2:Z)
{
fullimg=rbind(fullimg,rep(NA,diff(yrange)+1))
xrange=range(which(apply(img[,,i],1,sum.na)!=0))
fullimg=rbind(fullimg,img[xrange[1]:xrange[2],yrange[1]:yrange[2],i])
}
par(fin=10*c(dim(fullimg)/max(dim(fullimg))))

image(fullimg,zlim=zlim,axes=FALSE,col = tim.colors(64))
image.plot(fullimg,zlim=zlim,axes=FALSE,col = tim.colors(64),add=TRUE)
}
