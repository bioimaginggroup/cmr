% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_mbf.R
\name{imageMBF}
\alias{imageMBF}
\title{Plotting of (voxelwise) cardiac MBF}
\usage{
imageMBF(img, zlim = NULL, reverse = TRUE)
}
\arguments{
\item{img}{3d array ob MBF values}

\item{zlim}{limits of MBF, default: NULL means zlim=c(0,max(img,na.rm=TRUE))}

\item{reverse}{reverse color scheme}
}
\value{
plots
}
\description{
Plotting of (voxelwise) cardiac MBF
}
\examples{
data(cmrsim)
imageMBF(maxresp_sim)
}
