\name{gaussiansmooth}
\alias{gaussiansmooth}
\alias{gaussiansmooth.grid}
\title{
function to apply a gaussian smoothing to a grid.
}
\description{
This function returns a grid of the same dimension as the input in which each
cell has a gaussian weighted mean of cells around it in the  input grid. To
speed up processing with large bandwidths (SD) it may upscale the data prior to
applying a Gaussian kernel.  The amount of upscaling is determined by the
\code{kernel.dim} argument which sets the target dimension (in cells) of the
kernel at the courser resolution. The function can use a weights table to run
directly on a landcover grid substituting the cover values with weights as
defined in the weights table and with the \code{no.match} and \code{na.value}
arguments.  Edge correction is used; \code{NA} cells and cells outside of the
region are not included in the calculation.
}
\usage{
gaussiansmooth(x, sd, max.r, kernel.dim, ...)
\method{gaussiansmooth}{grid}(x, sd, max.r=3*sd, kernel.dim=21, max.p.na=1,
        sample.points, weights, no.match=0, use.old = FALSE,
        na.value = NA, ...)
}

\arguments{
  \item{x}{
	\code{x} should be a grid object in memory.
}
  \item{sd}{
	The standard deviation of the kernel used to do the smoothing (in map units).
}
  \item{max.r}{
	The maximum radius that weights are calculated for (in map units).
	It defaults to three times \code{sd}.
}
  \item{kernel.dim}{
      The dimension of the kernel in the upscaled grid. It must be an odd number.
      This controls the amount of upscaling, which is adjusted such that a
      \code{max.r} radius kernel has this cell dimension.  Smaller values result
      in more upscaling and faster processing at the expense of precision.
      Use \code{Inf} for no upscaling.  The default value of this parameter
      is 21.
}

\item{max.p.na}{
	If the proportion of \code{NA} values (from the original resolution) in an upscaled
  cell exceeds \code{max.p.na} than it will receive a \code{NA} value during
  upscaling; otherwise the upscaled cell is given the mean value of all
  non-NA cells that it covers.
}

 \item{sample.points}{
	(optional) a dataframe or matrix with columns "x" and "y" containing the map
  unit coordinates of points.  If specified the function will return both the
  samples and a grid (in a list).
}
\item{weights}{
	(optional) a dataframe or matrix with columns "value" and "weight".  If
  \code{weights} is provided than all the cells in the grid will be replaced
  with their associated weight; any cell whose value is not in the value column
  will be replaced with the values specified by \code{no.match} or \code{na.value}.
  If \code{weights} does not contain columns
  named "value" and "weight" and contains only two columns it will be assumed
  that the first is "value" and the second is "weight".
}

\item{no.match}{
If weights are supplied this determines the value to use when cells in \code{x}
 are not \code{NA} and do not match any of the values in \code{weights}.
 }

 \item{na.value}{
If weights are supplied this determines the value to use when cells in \code{x}
 are \code{NA}.  Note if the \code{weights} argument contains a value for
 \code{NA} as well that will take precedence.
 }

  \item{use.old}{if \code{FALSE} then the newer \pkg{Rcpp} based implementation
  is used. Otherwise the old, pure R, version is used.  This argument is here to
  facilitate testing.}


 \item{\dots}{
	arguments to be passed to other methods.
 }

}
\value{
If \code{sample.points} are not provided then the function will return a grid.
If \code{sample.points} are provided the grid method will return a two item
list: \code{sample} containing a vector of sampled values, and \code{grid}.
}

\author{
Ethan Plunkett
}

\note{
It is left up to the user to make sure that the \code{weights} table is
complete. Any values that occur in the grid but are not in the "value" column
of \code{weights} will be replaced with \code{no.match} prior to smoothing
unless they are \code{NA} in which case they will be replaced with \code{na.value}.
}

\seealso{
\code{\link{kernelsmooth}}, \code{\link{makegaussiankernel}},
\code{\link{upscale}}, and  \code{\link{downscale}} are all called by this function.
\code{\link{kernelsmooth}} does the same thing but allows for user specified
kernels.
\code{\link{calckernel}} allows calculating a kernel weighted average for a
single cell.
\code{\link{kernsource}} is similar but works by summing kernels centered on
each cell, it may be faster if the non-zero values in the grid are sparse and
the SD is small.
\code{swap} is used internally to substitute values when weights are used.
}
\examples{

# Matrix method
m <- matrix(1:60, 6, 10)
k <- matrix(1, 3, 3)
k[2,2] <- 3
k <- k/sum(k)
m2 <- kernelsmooth(m, k)


}
\keyword{ spatial }
