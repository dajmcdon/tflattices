#' Trend filter penalty matrix
#'
#' @param side_lengths vector of side lengths for the grid
#' @param orders vector of orders of the trend filter penalty (denoted `k`). 
#' @param wrap vector indicating any dimensions to "wrap". That is whether
#'   the left side wraps around to tough the right side (like a cylinder).
#'
#' @return a sparse matrix. See [Matrix::sparseMatrix()].
#' @export
#'
#' @examples
#' penalty_matrix(c(4,5), c(1,0))
penalty_matrix <- function(side_lengths, 
                           orders = rep(0L, length(side_lengths)),
                           wrap = rep(FALSE, length(side_lengths))) {
  d <- length(side_lengths)
  if ((no <- length(orders)) == 1) orders <- rep(orders, d)
  if ((nw <- length(wrap)) == 1) wrap <- rep(wrap, d)
  if (no != 1L && d != no) 
    stop("side_lengths and orders are not of compatible lengths.")
  if (nw != 1L && d != nw) 
    stop("side_lengths and wrap are not of compatible lengths.")
  side_lengths <- as.integer(side_lengths)
  orders <- as.integer(orders)
  wrap <- as.integer(wrap)
  stopifnot(all(side_lengths > 0), all(orders >= 0), all(wrap %in% c(0L, 1L)))
  D <- builddmat(side_lengths, orders, wrap)
  return(D)
}