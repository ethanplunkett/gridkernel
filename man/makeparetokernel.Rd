\name{makeparetokernel}
\alias{makeparetokernel}

\title{
function to make a Pareto kernel
}
\description{
This function creates a Pareto kernel for use with \code{\link{kernelsmooth}}.
}
\usage{
makeparetokernel(scale, shape, max.r = max.r, cellsize)
}

\arguments{
  \item{scale, shape}{
These control the shape of the kernel.  See \code{\link[texmex]{dgpd}} in the
\pkg{texmex} package. }

\item{max.r}{
The maximum radius for which which weights are calculated.
}
  \item{cellsize}{
The cellsize of the data to which the kernel will be applied.
}
}
\value{
this function returns a square matrix where each cell contains a weight that
is a Pareto function of the distance to the center of the kernel.

}
\author{
Ethan Plunkett
}

\seealso{
\code{\link{makegaussiankernel}},
A working example is in \code{\link{kernelsmooth}}
}
\keyword{ spatial}

