#' Bullseye plot for cmr
#'
#' @param werte 
#' @param zlim 
#' @param reverse 
#' @param text 
#'
#' @return plot
#' @export
#' @import fields
#'
bullseye<-function(werte,zlim=NULL,reverse=TRUE,text=TRUE,cex=2)
{
if(is.null(zlim))zlim=c(0,max(werte,na.rm=TRUE))

ziel0<-werte
ziel<-ziel0-zlim[1]
ziel<-ziel/(zlim[2]-zlim[1])
ziel[ziel>=1]<-1

colo<-rev(tim.colors(265)[-265])
if(reverse)colo<-tim.colors(265)[-265]
par(mai=c(0,0,0,0))
plot(c(0,1),c(0,1),col="white",axes=FALSE,asp=1,ylab="",xlab="")
 
for (i in 1:6)
polygon(0.4+0.4*c(cos(-pi/2+pi/6+seq(2*i,2*i+2,by=0.1)*pi/6),0.75*cos(-pi/2+pi/6+seq(2*i+2,2*i,by=-0.1)*pi/6)),
        0.5+0.4*c(sin(-pi/2+pi/6+seq(2*i,2*i+2,by=0.1)*pi/6),0.75*sin(-pi/2+pi/6+seq(2*i+2,2*i,by=-0.1)*pi/6)),
col=colo[263*ziel[i]+1],border="black")
for (i in 1:6)
polygon(0.4+0.3*c(cos(-pi/2+pi/6+seq(2*i,2*i+2,by=0.1)*pi/6),0.667*cos(-pi/2+pi/6+seq(2*i+2,2*i,by=-0.1)*pi/6)),
        0.5+0.3*c(sin(-pi/2+pi/6+seq(2*i,2*i+2,by=0.1)*pi/6),0.667*sin(-pi/2+pi/6+seq(2*i+2,2*i,by=-0.1)*pi/6)),
col=colo[263*ziel[i+6]+1],border="black")
for (i in 1:4)
polygon(0.4+0.2*c(cos(seq(2*i-1,2*i+1,by=0.1)*pi/4),0.5*cos(seq(2*i+1,2*i-1,by=-0.1)*pi/4)),
        0.5+0.2*c(sin(seq(2*i-1,2*i+1,by=0.1)*pi/4),0.5*sin(seq(2*i+1,2*i-1,by=-0.1)*pi/4)),
col=colo[263*ziel[i+12]+1],border="black")


#for (i in 0:262)
#polygon(c(0.88,0.93,0.93,0.88),c(0.2+i*.6/263,0.2+i*.6/263,0.2+(i+1)*.6/263,0.2+(i+1)*.6/263),col=colo[i+1],border="transparent")
#polygon(c(0.88,0.93,0.93,0.88),c(0.2,0.2,0.8,0.8))

#bla<-seq(zlim[2],zlim[1],length=5)
#for (i in 0:4)
#text(0.93,0.8-i*0.15,round(bla[i+1],1),cex=1,pos=4)
#dev.off()

if(text){
  par(cex=cex)
text(0.88,0.5,"LCX")
text(0.1,0.85,"LAD")
text(0.1,0.15,"RCA")
}
image.plot(fullimg,zlim=zlim,legend.only=TRUE,legend.width=1.8,add=TRUE)
}

bullseye.1<-function(werte,zlim=NULL,reverse=TRUE)
{

if(is.null(zlim))zlim=c(0,max(werte,na.rm=TRUE))

ziel0<-werte
ziel<-ziel0-zlim[1]
ziel<-ziel/(zlim[2]-zlim[1])
ziel[ziel>=1]<-1

colo<-tim.colors(64)
if(reverse)colo<-rev(tim.colors(64))
par(mai=c(0,0,0,0))
plot(c(0,1),c(0,1),col="white",axes=FALSE,asp=1,ylab="",xlab="")
 
for (i in 1:6)
{
polygon(0.45+0.4*c(cos(-pi/2+pi/6+seq(2*i,2*i+2,by=0.1)*pi/6),0.66*cos(-pi/2+pi/6+seq(2*i+2,2*i,by=-0.1)*pi/6)),
        0.55+0.4*c(sin(-pi/2+pi/6+seq(2*i,2*i+2,by=0.1)*pi/6),0.66*sin(-pi/2+pi/6+seq(2*i+2,2*i,by=-0.1)*pi/6)),
col=colo[63*ziel[i]+1],border="black")
}
#dev.off()
}

pseudobullseye<-function(mbf,zlim=c(0,6),text=TRUE,reverse=FALSE,center=FALSE){
if (center)
  {
    sum.na<-function(x)sum(x,na.rm=TRUE)
    x<-apply(mbf,1,sum.na)!=0
    xx<-sum(x)
    y<-apply(mbf,2,sum.na)!=0
    yy<-sum(y)
    xy<-36
    mbf.neu<-array(NA,c(xy,xy,3))
  for (i in 1:3)
  {
    x<-apply(mbf[,,i],1,sum.na)!=0
    xi<-sum(x)
    y<-apply(mbf[,,i],2,sum.na)!=0
    yi<-sum(y)
    x0<-floor((xy-xi)/2)
    y0<-floor((xy-yi)/2)
    mbf.neu[x0+1:xi,y0+1:yi,i]=mbf[x,y,i]
  }
    mbf<-mbf.neu
  }
par(mar=c(0,0,0,0))
par(pty="s",bty="n")
co=tim.colors(64)
if(reverse)co=rev(tim.colors(64))
image(mbf[,,3],zlim=zlim,col=co,axes=FALSE)
if(text){
  par(cex=2)
text(0.97,0.48,"LCX")
text(0.15,0.83,"LAD")
text(0.15,0.13,"RCA")
}
#par(plt=c(.16,.84,.13,.81),new=TRUE)
par(plt=c(.15,.85,.15,.85),new=TRUE)
par(pty="s",bty="n")
image(mbf[,,2],zlim=zlim,axes=FALSE,col=co)
#par(plt=c(.30,.74,.24,.68),new=TRUE)
par(plt=c(.3,.7,.3,.7),new=TRUE)
par(pty="s",bty="n")
image(mbf[,,1],zlim=zlim,axes=FALSE,col=co)

}
