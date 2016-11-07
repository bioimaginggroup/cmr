#' Title
#'
#' @param mbf 
#' @param zlim 
#' @param text 
#' @param reverse 
#' @param center 
#'
#' @return
#' @export
#'
#' @examples
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
