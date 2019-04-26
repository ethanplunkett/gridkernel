\name{calckernel}
\alias{calckernel}
\alias{calc.kernel}
\alias{calckernel.grid}
\alias{calckernel.matrix}
\title{
 functions to calculate the sum of the product of a kernel and a patch of a grid or matrix around a focal cell.
}
\description{
  This function returns the sum of the product of each cell in the kernel and the corresponding cell in the matrix around the focal cell.
}
\usage{
\method{calckernel}{grid}(obj, kernel, x, y, row=y2r(y, list=obj), col=x2c(x, list=obj), ...)
\method{calckernel}{matrix}(obj, kernel, row, col, use.old=FALSE, ...)
calckernel(obj, ...)
calc.kernel(obj, ...) # depreciated
}
\arguments{
  \item{obj}{
either a matrix or an object of class "grid" containing the data
}
  \item{kernel}{
the kernel. It should be a square, odd dimensioned, numeric matrix that sums to 1.
}
  \item{row, col}{
the row and column of the focal cell.  The center of the kernel will be aligned with this cell.
}
  \item{x, y}{
the x and y coordinates in map units of the focal cell.  
}
  \item{use.old}{if \code{FALSE} then \pkg{Rcpp} based implementation is used. Otherwise an older, pure R, version is used.  This is here to facilitate testing.}
 \item{...}{
Arguments to be passed to methods
}

}
\details{
this function returns the sum of the product of the kernel and data around the focal cell. If (as required) 
the kernel sums to one this is the weighted average of the cells.
If the kernel overlaps the edge or \code{obj} contains \code{NA} values than only the available data will be used and the 
kernel will be renormalized such that the sum of the used cells is one; the result will be the weighted average of the available data.
}
\value{
A weighted average of the cells in the matrix surrounding the focal cell. 
}
\author{
Ethan Plunkett
}
\note{
Please make sure the kernel sums to 1.  This could be done in the function but 
that would be inefficient in cases where the kernel is repeatedly applied to 
the matrix.  Use \code{kernel <- kernel/sum(kernel, na.rm = TRUE)} to properly
normalize a kernel. If the kernel doesn't sum to one than the results near
edges and \code{NA} values will be inconsistent.

}
\seealso{
\code{\link{kernelsmooth}}, \code{\link{makegaussiankernel}}
}
\examples{
m <- matrix(1:60, 6, 10)
k <- matrix(1, 3, 3)
k[2,2] <- 4
k <- k/sum(k)

kernelsmooth(m, k)

calckernel(m, k, 2, 2)

}
\keyword{ spatial }

