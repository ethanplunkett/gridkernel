\name{kernelscalingdetails}
\alias{kernelscalingdetails}
\title{
Determines the upscaling factor, the realized tilesize, and the required buffer used by gaussiansmooth and kernelsmooth when coarsening grids to increase efficiency
}
\description{
This is used internally in \code{\link{kernelsmooth}} and \code{\link{gaussiansmooth}} to determine the upscaling factor, realized tilesize, and buffer.
It is provided to users as it is necessary to anticipate the realized tilesize  correctly to use the tiles arguments to the character methods of the smoothing functions (in \pkg{gridio}).
}
\usage{
kernelscalingdetails(max.r, cellsize, tilesize = 2000, kernel.dim = 21, method = 3)
}

\arguments{
  \item{max.r}{ radius of the kernel in the units of the projection (mapunits)
  it should be set to the same value as in the call to the smoothing function (\code{gaussiansmooth}, \code{kernelsmooth}.
  Note for gaussiansmooth code{max.r} defaults to 3 times the sd.
}
  \item{cellsize}{ cellsize of the grid being smoothed
}
  \item{tilesize}{ the target tilesize as specified in calls to the smoothing functions.
}
  \item{kernel.dim}{ the target kernel dimension (in the coarse grid) as specified in the calls to the smoothing funcitons.
}
  \item{method}{ The method used when calling gaussiansmooth.character (in \pkg{gridio}) in all other cases just leave at the default.
}
}
\details{
If you want to run on a subset of tiles and are using gaussiansmooth.character or kernelsmooth.character in the gridio package then use this function to precalculate the tilesize so you can figure out which tiles to run on.  The tilesize will depend on the target tilesize, the cellsize, and max.r.  Which in practice means you'll have a different tilesize for each kernel bandwidth.
}
\value{
 \item{factor }{The factor used when changing the resolution of a grid.  This is the ratio of the dimension of a coarse grid cell to that of a fine grid cell. Eg with 30 metercells and factor of 4 the coarse representation will have 120 m cells}
  \item{tilesize}{This is the realized tilesize used by the smoothing functions (when reading from disk in tiles).  Since the realized tilesize is a multiple of the factor the user requested tilesize can't usually be met and instead the closes possible multiple of factor is used.}
}

\author{
Ethan Plunkett
}

\seealso{
\code{\link{gaussiansmooth}},\code{\link{kernelsmooth}}
}
