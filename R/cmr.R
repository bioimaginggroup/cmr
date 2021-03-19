#' Bayesian analysis of cardiovascular magnetic resonance imaging
#'
#' @param data 3D or 4D array of CMR signal 
#' @param input input function
#' @param mask 2D or 3D array of mask. Voxel with 0 or FALSE will be ommited from anaysis. Default NULL: use NA values in data as mask
#' @param method "spatial" or "local"
#' @param quantiles quantiles used for credible intervall, default: c(0.25, 0.75)
#' @param cores number of cores for parallel computation. Spatial model only computes slices parallel, local can be parallelized on voxel level
#' @return list of mbf (point estimation) and ci (credible interval)
#' @importFrom parallel mclapply
#' @export

cmr<-function(data, input, mask=NULL, method="spatial", quantiles=c(.25,.75), cores=getOption("mc.cores", 1))
{
  if (length(dim(data))==4)
  {
    mbf <- ci <- array(NA,dim(data)[1:3])
    I <- 1:dim(mbf)[3]
    if (cores>1)temp<-parallel::mclapply(I,do.cmr,data,input,mask,method,cores=floor(cores/dim(mbf)[3]),mc.cores=cores)
    if (cores==1)temp<-lapply(I,do.cmr,data,input,mask,method,cores=1)
    
    for (i in I)
      {
        mbf[,,i]=t(as.matrix(temp[[i]]$mbf))
        ci[,,i]=t(as.matrix(temp[[i]]$ci))
    }
  }
  if (length(dim(data))==3)
  {
    mbf <- ci <- array(NA,dim(data)[1:2])
    if (method=="local")
    {
      if (is.null(mask))
      {
        mask=array(NA,c(30,30))
        mask[data[,,i]!=0]=1
      }
      temp=cmr.local(data,mask,input)
      mbf=t(as.matrix(temp$mbf))
      ci=t(as.matrix(temp$ci))
    }
    if (method=="spatial")
    {
      if (is.null(mask))
      {
        mask=array(NA,c(30,30))
        mask[data[,,i]!=0]=1
      }
      temp=cmr.space(data,mask,input)
      mbf=t(as.matrix(temp$mbf))
      ci=t(as.matrix(temp$ci))
    }
  }
  if (length(dim(data))>4){return(NULL)}
  if (length(dim(data))<2){return(NULL)}
  
  return(list("mbf"=mbf,"ci"=ci))
}

do.cmr<-function(i,data,input,mask,method,cores)
{
  if (is.null(mask))
  {
    mask0=array(NA,c(30,30))
    mask0[data[,,i,1]!=0]=1
  }
  else
  {
    mask0<-mask[,,i,]
  }
  if (method=="local")temp=cmr.local(data[,,i,],mask0,input,cores=cores)
  if (method=="spatial")temp=cmr.space(data[,,i,],mask0,input)
  return(temp)
}
