\name{kernelsmooth}
\alias{kernelsmooth}
\alias{kernelsmooth.grid}
\alias{kernelsmooth.matrix}
\title{
function to apply a kernel to each cell in a grid or matrix with edge and NA correction
}
\description{
These functions apply a kernel to a every cell in a grid or matrix. The returned
object will be of the same type and dimension  of the input but every cell will
contain the sum of the product of the kernel and the corresponding cells of
the input. If the \code{kernel.fun} argument is used than depending on the value
of \code{max.r} and \code{kernel.dim} the grid method may perform upscaling
before applying the kernel, this makes the algorythm much faster for large
bandwidth smoothing operations.
}
\usage{
kernelsmooth(x, ...)

\method{kernelsmooth}{grid}(x, kernel, ...)
\method{kernelsmooth}{matrix}(x, kernel, use.old=FALSE, ...)
\method{kernelsmooth}{grid}(x, kernelfun, max.r, kernel.dim = 21, max.p.na = 1,
sample.points, weights, no.match = 0, use.old = FALSE, na.value = NA, ...)

}
\arguments{
  \item{x}{
   An object of class "grid" or "matrix" containing the raster data to be
   smoothed.
}
  \item{kernel}{
  A kernel. It should be a square, odd dimensioned, numeric matrix that sums to
  1.Use either this argument or \code{kernelfun}.
}
 \item{kernelfun}{
  A function that generates a square, odd dimensioned, matrix that sums to 1.
  It must have \code{max.r} and \code{cellsize} arguments.  Other arguments to
  \code{kernelfun} can be passed via \code{...}.
  Use either this argument or \code{kernel},  if you use both \code{kernel}
  will be ignored.   See \code{\link{makegaussiankernel}} for one such function.
}
  \item{max.r}{
  The maximum radius that weights are calculated for (in map units) - will be
  passed to \code{kernelfun} and used along with \code{kernel.dim}
  to determine how much upscaling is appropriate.
}
  \item{kernel.dim}{
  The cell dimension of the kernel in the upscaled grid. It must be an odd
  number.  This controls the amount of upscaling, which is adjusted such that a
  \code{max.r} radius kernel has this cell dimension.  Smaller values result in
  more upscaling and faster processing at the expense of precision.  The default
  value 21 is a good trade off between speed and precision for guassian kernels
  but may not be appropriate  for other kernels. See
}
  \item{max.p.na}{
	If the proportion of NA values (from the original resolution) in an upscaled
	cell exceeds \code{max.p.na} than it will receive an NA value during
	upscaling; otherwise the upscaled cell is given the mean value of all
	non-NA cells that it covers.
}

 \item{sample.points}{
	(optional) a dataframe or matrix with columns "x" and "y" containing the map
	unit coordinates of points.  If specified the character method will return
	the a vector of sample values (from the smoothed surface), and the grid
	method will return both the samples and a grid (in a list).
}
\item{weights}{
	(optional) a dataframe or matrix with columns "value" and "weight".
	If \code{weights} is provided than all the cells in the grid will be replaced
	with their associated weight; any cell whose value is not in the value column
	will be replaced with \code{no.match}, NA cells in \code{x} will be replaced
	with \code{na.value}.  If \code{weights} does not contain columns named
	"value" and "weight" and contains only two columns it will be assumed that
	the first is "value" and the second is "weight".
}
\item{no.match}{
  If weights are supplied Non-NA values in x that don't match any value in
  original will be replaced with this value.
}
\item{na.value}{
 If weights are supplied this determines the value to use when cells in \code{x}
 are \code{NA}.  Note if the \code{weights} argument contains a value for
 \code{NA} that will take precedence.
 }

 \item{use.old}{if \code{FALSE} then the newer \pkg{Rcpp} based implementation is used.
  Otherwise the old, pure R, version is used. This argument is here to
  facilitate testing;  leave it at FALSE.}

\item{\dots}{
Arguments to be passed to other methods. If x is a character (path) these are
passed on to \code{kernelfun}.
}

}
\details{
If a gaussian kernel is used than the result is a smoothed version of the
original grid or matrix. Other kernels could be used for other types of
smoothing, or potentially be used to enhance edges or increase contrast.

This function adjusts for missing values in the matrix and at the edges.

The matrix version is equivalent to calling \code{\link{calckernel}} on each
cell of \code{x} to calculate the corresponding value in  the output.

If the \code{kernel} argument is  used the grid version is the same but returns
a grid object. If the \code{kernelfun} argument is used than the grid version
may upscale \code{x} to a coarser resolution, make a kernel appropriate for the
new cellsize, apply the kernel to the upscaled grid and then downscale
(with interpolation) back to the original resolution.


}

\value{
an object of the same class as \code{x} in which each cell is the result of
applying a kernel to the correspononding focal cell in \code{x}, or if
\code{sample.points} are provided a list with a \code{grid} and the sampled
values.

}
\author{
Ethan Plunkett
}


\seealso{
\code{\link{calckernel}}, \code{\link{makegaussiankernel}}

}
\examples{


# Matrix example
m <- matrix(1:110, 10, 11)
m[5, ] <- 80
k1 <- matrix(1, 3, 3)  # 3 x 3 kernel with more weight in center
k1[2,2] <- 3
k1 <- k1/sum(k1)
sm1 <- kernelsmooth(m, k1)

# An alternative to gaussiansmooth
k2 <- makegaussiankernel(sd = 1.5,  )
sm2 <- kernelsmooth(m, k2)



\dontrun{
# Visualize the results
# Note image rotates the matrix 90 CCW when plotting it.
cols <- heat.colors(n = 110)
image(m, col = cols)  # original
image(sm1, col = cols) # 3 x 3 with extra weight in center
image(sm2, col = cols)  # 1.5 cell
}

# Grid example
g  <- as.grid(m, xll = 0, yll = 0, cellsize = 10)
gs1  <- kernelsmooth(g, k1)
gs2 <- kernelsmooth(g, k2)


\dontrun{
# Visualize
plot(g)
plot(gs1)
plot(gs2)
}

# With upscaling
# Upscaling only comes into play with large kernel bandwidth and a
# kernel function. Here will make  a big example landscape and
#  use make.gaussian.kernel as the kernel function the result will be the same
# as if we had used gaussiansmooth()

dim  <- 1000  # cells
m <- matrix(runif(dim^2,  min = 0, max  = 100), dim, dim)
m[500:520, ] <- 0  # Low value horizontal band
m[ , 500:520] <- 100   # High value vertical band
m[1:100, 1:100] <- NA  # Patch of NA in upper left
g <- as.grid(m, xll = 0, yll = 0, cellsize = 30)

# Do smoothing  with  upscaling
sd <- 700  # meters not cells
a <- system.time({
  gs <- kernelsmooth(g, kernelfun = makegaussiankernel,  max.r = 3 * sd, sd = sd  )
})


\dontrun{
# Compare to results without upscaling
#  SLOW!!
b <- system.time({
  gs2 <- kernelsmooth(g, kernelfun = makegaussiankernel, kernel.dim = Inf,
                      max.r = 3 * sd, sd = sd  )
})
(a/b)[3]   #  0.007
cor(as.numeric(gs$m), as.numeric(gs2$m), use = "complete")  # 0.99969
mean( as.numeric(gs$m  - gs2$m), na.rm = TRUE)  # 0.0013
# Upscaling took 0.007 of the time and is almost  perfectly correlated to the
# non-upscaled version.
# However, upscaling did produce some artifacts in the interior of the NA region.
}

#  Pareto Kernel example
pk <- makeparetokernel(scale=10, shape=.8, max.r=250, cellsize=30)

# Plot a cross  section of the kernel
plot(k[ , ceiling(nrow(k)/2)])

\dontrun{
# By passing the kernel directly:
ps <- kernelsmooth(g, pk)

# With a kernel function
# Note, pareto  kernels are pointy, you might want to increase kernel.dim
ps2 <- kernelsmooth(g, kernelfun = makeparetokernel, max.r = 250, scale = 10, shape = 0.8 )
}

}
\keyword{ spatial }

