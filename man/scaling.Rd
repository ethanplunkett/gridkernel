\name{scaling}
\alias{downscale}
\alias{upscale}
\title{
functions to change the resolution of a grid
}
\description{
These functions change the resolution of a grid. They can downscale and upscale based on integer scaling of the resolution. 
}
\usage{
downscale(x, factor, match.extent = TRUE)
upscale(x, factor, max.p.na = 0.2, use.old = FALSE)
}
\arguments{
  \item{x}{
A grid object
}
  \item{factor}{
The factor by which the horizontal and vertical resolution will be either increased or decreased. 
}
  \item{match.extent}{
   Downscaling results in a loss of pixels around the edge.  If match extent is TRUE then the downscaled image will have the same extent as x and will include a band of NA pixels around the edge. If it is FALSE the extent will be smaller and will not include the NA cells.
}
\item{max.p.na}{
	If the proportion of NA values (from the original resolution) in an upscaled cell exceeds \code{max.p.na} than it will receive an \code{NA} value during	 upscaling; otherwise the upscaled cell is given the mean value of all non-NA cells that it covers. 
}

 \item{use.old}{if \code{FALSE} then the newer \pkg{Rcpp} based implementation is used. 
  Otherwise the old, pure R, version is used.  This argument is here to 
  facilitate testing.}

}

\value{
A grid of similar extent containing the upscaled or downscaled version of x.
}
\author{
Ethan Plunkett
}
\note{

Downscaling with the \code{match.extent} argument set to TRUE can gaurantee that the new grid is the same extent as the lower resolution grid. However, upscaling often expands the extent of the grid slightly because cells that are partially filled will be created in the upscaled grid. So when used in combination you will usually end up with mismatched extents.

Calling the function \code{\link{matchextent}} can fix this.  For an example see the code for the \code{\link{gaussiansmooth.grid}} function.
} 


\seealso{
\code{\link{match.extent}}
}
\examples{
#
}
\keyword{ spatial }

