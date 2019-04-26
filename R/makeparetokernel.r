makeparetokernel<-function(scale, shape, max.r = max.r,
                             cellsize) {

  max.r.cells <- max.r/cellsize
  size = ceiling(max.r.cells) * 2 + 1
  center = ceiling(max.r.cells) + 1
  kernel <- new("matrix", 0, size, size)

  rows <- row(kernel)
  cols <- col(kernel)
  dists <- sqrt( (rows - center)^2 + (cols - center)^2 )

  kernel[ , ] <- texmex::dgpd(dists, sigma = scale, xi = shape, log = FALSE)
  kernel[ dists > max.r] <- 0

  kernel <- kernel/sum(kernel)
  # This last part deletes the cells at the edge if they are all zero
  if (all(kernel[1, ] == 0, kernel[, 1] == 0,
          kernel[nrow(kernel),] == 0, kernel[, ncol(kernel)] == 0))
    kernel <- kernel[2:(nrow(kernel) - 1), 2:(ncol(kernel) - 1)]
  return(kernel)
}
