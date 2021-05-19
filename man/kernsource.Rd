\name{kernsource}
\alias{kernsource}

\title{
 	Source-based standard (non-resistant) Gaussian kernel estimator
}
\description{
	This function uses Gaussian kernels to calculate a smooth surface by summing kernels centered on all cells that have non-zero values. 

}
\usage{
kernsource(x, bandwidth, cellsize = 30, search = 3, buffer = 0)
}

\arguments{
  \item{x}{
source grid (matrix) with >=0 values for source points.  Kernels are weighted by values in x; use 1's for unweighted kernels.
}
  \item{bandwidth}{
bandwidth (h), in map units (the standard deviation of the Gaussian function)
}
  \item{cellsize}{
	cellsize in map units
}
  \item{search}{
search distance, in s.d. search = 3 gives 99% of kernel; use values as low as 2 for quick and sloppy results.
}
  \item{buffer}{
edge buffer, in meters.  Source cells in buffer are ignored (presumably these are dealt with in adjacent blocks).  Buffer must be either 0 or >= bandwidth * search / cellsize.
}
}
\details{
This version does not do edge corrections.   Results will be biased downward at hard data edges.  Source data should be buffered by bandwidth * search / cellsize; write results in summation mode.

B. Compton, 28 Apr 2011.  From kern.r, 16 Jan 2008; translated from APL function KERN, 25 Sep. 2003.
}
\value{
A matrix with the a smoothed surface.
}

\author{
Brad Compton
}
\note{
From Brad Compton (5/28/2011)
Ethan and I spent a while on the theory and practice of non-resistant kernels yesterday.  Here's what we came up with:

There are two approaches to estimating kernels: source-based and recipient-based.  Both give the same results, but there are significant differences in efficiency depending on the situation.  I'll use an example: we're building kernels around developed cells with weights for each of the four development classes. The source-based approach, familiar from home-range analysis, centers a kernel on each developed cell, multiplying the kernel by the cell's weight.  It then adds the kernel to the result, cell-by-cell.  The recipient-based approach (commonly used in image processing), on the other hand, visits each cell in the landscape for which you want a result, centers a kernel on that cell, multiplies the kernel value by the weights it overlaps, and then sums the volume of result and places it in the focal cell.  Both give you the same results, but when you have a small number of source cells (as in a home range analysis or the development example), a source-based approach is far more efficient, because run-time is proportional to the number of source cells.  If, on the other hand, you have a restricted set of cells that you want a result for and lots of source cells (say you want values for each vernal pool based on how much forest surrounds it), then the recipient-based approach is faster.

The other issue is edge correction.  This needs to be done if you have hard data edges in your landscape so you don't bias results downward near the edges.  I think this will not be an issue for the LCC project, except on the Maine-Quebec/New Brunswick border, and this won't really matter for development modeling, because there's nothing much on either side.  Edge correction is fairly costly, so we don't want to do it if it's not necessary.

I've attached kernsource.r, which is a source-based version without edge correction.  This is a version of my old kern.r with the arguments changed to be understandable by a human-being, instead of matching bkde2D's incomprehensible arguments.  Someday, I'll add edge-correction to this, but I don't think we need it now.  Liz, I presume you've got data beyond the Hampden County edge and will then clip the results to Hampden County.  I hope this version is acceptably fast for large kernels.  I'd try bandwidths up to 10 or 20 km or so.  Let me know how long this ends up taking.  I just did a quick test on a 1000x2000 matrix with 200 1's and a bandwidth of 10 km; it took a few seconds.

Ethan has a version that's recipient-based, with edge correction.  For applications where you only want to get results for a restricted set of cells and there are lots of source cells, this will be faster.  And where we need edge correction, this is the only version to use until I add edge correction.

If you're calling either version in blocks for a large landscape, you'll need to read a buffer around each block of at least bandwidth (in cells) * search distance.  For kernsource, you'll write  the entire block+buffer in transparent summation mode (results are added to the grid), for kern.recipient, you'll throw away the buffer when writing in normal replacement mode.

At some point, Ethan and I will integrate these two functions into one function with options for the various approaches.

}


\seealso{
\code{\link{calckernel}}, \code{\link{kernelsmooth}}
}
\examples{
#
}
\keyword{ spatial }

