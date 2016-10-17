#' Title Plotting of (voxelwise) cardiac MBF 
#'
#' @param img 3d array ob MBF values
#' @param zlim limits of MBF, default: NULL means zlim=c(0,max(img,na.rm=TRUE))
#' @param reverse reverse color scheme
#'
#' @return plots
#' @export
#' @import fields graphics
#'
#' @examples
#' library(cmr)
#' data(sim)
#' image.mbf(resp)
imageMBF<-function(img,zlim=NULL,reverse=TRUE)
{
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
fullimg=rbind(fullimg,rep(NA,diff(yrange)+1))
xrange=range(which(apply(img[,,i],1,sum.na)!=0))
fullimg=rbind(fullimg,img[xrange[1]:xrange[2],yrange[1]:yrange[2],i])
}
fullimg=rbind(fullimg,rep(NA,diff(yrange)+1))
fullimg=rbind(fullimg,rep(NA,diff(yrange)+1))
par(pin=5*c(dim(fullimg)/max(dim(fullimg))))
#par(mai=c(0,0,0,0))
farbe=tim.colors(64)
if(reverse)farbe=rev(farbe)
#add  space for legend
for (i in 1:floor(dim(fullimg)[1]/6))fullimg <- rbind(fullimg,rep(NA,diff(yrange)+1))
image(fullimg,zlim=zlim,axes=FALSE,col = farbe)
image.plot(fullimg,zlim=zlim,legend.only=TRUE,legend.width=1.8,add=TRUE)
}

