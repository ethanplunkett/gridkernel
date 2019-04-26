\name{makegaussiankernel}
\alias{makegaussiankernel}
\title{
function to make a gaussian kernel
}
\description{
this function makes a kernel where the value in each cell is a gaussian function of the distance between that cell and the center of the kernel.
}
\usage{
makegaussiankernel(sd, max.r = 3 * sd, cellsize = 1, include.center.cell = TRUE)
}
\arguments{
  \item{sd}{
the standard deviation in map units of the gaussian function
}
  \item{max.r}{
the maximum radius of the kernel in map units.  This defaults to three standard deviations. A smaller radius speeds calculation but loses information at the tail.
}
  \item{cellsize}{
the size of each cell in the grid in map units
}
  \item{include.center.cell}{
if FALSE the value of the center cell will be forced to zero.  This was helpful in my original use of the function but probably not to you.
}
}
\value{
a square, numeric matrix which sums to 1 and in which the value of each cell is a gaussian function of the distance from the cell to the center of the matrix.
}
\author{
Ethan Plunkett
}
\seealso{
\code{\link{calckernel}}, \code{\link{kernelsmooth}}
}
\examples{
k <- makegaussiankernel(sd=2)
m <- matrix(1:60, 6, 10)
sm <- kernelsmooth(m, k)

}
\keyword{ spatial }

