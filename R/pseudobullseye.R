#' Pseudo bullseye plot
#'
#' @param x 3D array
#' @param lim limits of x values
#' @param reverse boolean, reverse colors?
#' @param legend boolean, add legend?
#' @param text boolean, should text legend be added?
#' @param cex cex for text legend
#' @param legend.width Width in characters of the legend strip.
#' @param center boolean, should input x be centered before plotting
#'
#' @return plots
#' @export
#' @importFrom graphics par
#' @importFrom fields tim.colors image.plot
#' @examples
#' data(cmrsim)
#' pseudobullseye(maxresp_sim)
#' 
pseudobullseye<-function(x, lim=range(x,na.rm=TRUE), legend=FALSE, text=TRUE, reverse=FALSE, center=TRUE, cex=1, legend.width=1){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  mbf=x
  zlim=lim
  if (center)
  {
    x1<-apply(mbf,1,sum,na.rm=TRUE)!=0
    xx<-sum(x1)
    y<-apply(mbf,2,sum,na.rm=TRUE)!=0
    yy<-sum(y)
    mbf.neu<-array(NA,dim(x))
    for (i in 1:dim(x)[3])
    {
      x1<-apply(mbf[,,i],1,sum,na.rm=TRUE)!=0
      xi<-sum(x1)
      y<-apply(mbf[,,i],2,sum,na.rm=TRUE)!=0
      yi<-sum(y)
      x0<-floor((dim(x)[1]-xi)/2)
      y0<-floor((dim(x)[2]-yi)/2)
      mbf.neu[x0+1:xi,y0+1:yi,i]=mbf[x1,y,i]
    }
    mbf<-mbf.neu
  }
  par(pty="s",bty="n")
  par(plt=c(.2,.9,.15,.85))
  co=tim.colors(64)
  if(reverse)co=rev(tim.colors(64))
  image(mbf[,,3],zlim=zlim,col=co,axes=FALSE)
  if(text){
    par(cex=cex)
    text(0.97,0.48,"LCX")
    text(0.15,0.83,"LAD")
    text(0.15,0.13,"RCA")
  }
  fields::image.plot(mbf,zlim=zlim,legend.only=TRUE,legend.width=legend.width,add=FALSE, horizontal = TRUE)
  par(plt=c(.31,.79,.26,.74),new=TRUE)
  par(pty="s",bty="n")
  image(mbf[,,2],zlim=zlim,axes=FALSE,col=co)
  #par(plt=c(.30,.74,.24,.68),new=TRUE)
  par(plt=c(.38,.72,.33,.67),new=TRUE)
  par(pty="s",bty="n")
  image(mbf[,,1],zlim=zlim,axes=FALSE,col=co)
}
