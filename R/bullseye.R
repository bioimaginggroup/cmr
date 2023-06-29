#' Bullseye plot
#'
#' @param x vector of length 16 or 17
#' @param lim limits of x values
#' @param reverse boolean, reverse colors?
#' @param legend boolean, add legend?
#' @param text boolean, should text legend be added?
#' @param cex cex for text legend
#'
#' @return plot
#' @export
#' @import fields 
#' @importFrom graphics polygon par plot 
#' @importFrom plotrix draw.circle
#' @examples
#'   bullseye(1:16)
bullseye<-function(x, lim=NULL,reverse=TRUE,legend=TRUE,text=TRUE,cex=1)
{
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))

if (length(x)<16|length(x)>17)stop("x has to be of length 16 or 17")
werte=x
zlim=lim

if(is.null(zlim))zlim=c(0,max(werte,na.rm=TRUE))

ziel0<-werte
ziel<-ziel0-zlim[1]
ziel<-ziel/(zlim[2]-zlim[1])
ziel[ziel>=1]<-1

colo<-rev(fields::tim.colors(265)[-265])
if(reverse)colo<-fields::tim.colors(265)[-265]
par(mai=c(0,0,0,0))
plot(c(0,1),c(0,1),col="white",axes=FALSE,asp=1,ylab="",xlab="")
 
for (i in 1:6)
  graphics::polygon(0.4+0.4*c(cos(-pi/2+pi/6+seq(2*i+2,2*i+4,by=0.1)*pi/6),
                              0.75*cos(-pi/2+pi/6+seq(2*i+4,2*i+2,by=-0.1)*pi/6)),
        0.5+0.4*c(sin(-pi/2+pi/6+seq(2*i+2,2*i+4,by=0.1)*pi/6),
                  0.75*sin(-pi/2+pi/6+seq(2*i+4,2*i+2,by=-0.1)*pi/6)),
col=colo[263*ziel[i]+1],border="black")
for (i in 1:6)
  graphics::polygon(0.4+0.3*c(cos(-pi/2+pi/6+seq(2*i+2,2*i+4,by=0.1)*pi/6),0.667*cos(-pi/2+pi/6+seq(2*i+4,2*i+2,by=-0.1)*pi/6)),
        0.5+0.3*c(sin(-pi/2+pi/6+seq(2*i+2,2*i+4,by=0.1)*pi/6),0.667*sin(-pi/2+pi/6+seq(2*i+4,2*i+2,by=-0.1)*pi/6)),
col=colo[263*ziel[i+6]+1],border="black")
for (i in 1:4)
graphics::polygon(0.4+0.2*c(cos(seq(2*i-1,2*i+1,by=0.1)*pi/4),0.5*cos(seq(2*i+1,2*i-1,by=-0.1)*pi/4)),
        0.5+0.2*c(sin(seq(2*i-1,2*i+1,by=0.1)*pi/4),0.5*sin(seq(2*i+1,2*i-1,by=-0.1)*pi/4)),
col=colo[263*ziel[i+12]+1],border="black")
if (length(x)==17)
  plotrix::draw.circle(0.4,0.5,0.1,col=colo[263*ziel[17]+1],border="black")


#for (i in 0:262)
#polygon(c(0.88,0.93,0.93,0.88),c(0.2+i*.6/263,0.2+i*.6/263,0.2+(i+1)*.6/263,0.2+(i+1)*.6/263),col=colo[i+1],border="transparent")
#polygon(c(0.88,0.93,0.93,0.88),c(0.2,0.2,0.8,0.8))

#bla<-seq(zlim[2],zlim[1],length=5)
#for (i in 0:4)
#text(0.93,0.8-i*0.15,round(bla[i+1],1),cex=1,pos=4)

if(text){
  par(cex=cex)
text(0.88,0.5,"LCX")
text(0.1,0.85,"LAD")
text(0.1,0.15,"RCA")
}
if(legend)fields::image.plot(werte,zlim=zlim,legend.only=TRUE,legend.width=1.8,add=TRUE,cex=cex)
}
