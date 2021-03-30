#' sv_cbind
#'
#' Helper function called by `densify` used to bind
#' a list of sparse vectors into a sparse matrix columnwise
#' @param ... A series of sparse vectors
#' @return A sparse matrix, created columnwise from the vectors
#' @importFrom Matrix sparseMatrix
#' @importClassesFrom Matrix dsparseVector
#' @importFrom methods as slot
#' @source https://stackoverflow.com/questions/8843700/creating-sparse-matrix-from-a-list-of-sparse-vectors/8979207#8979207

sv_cbind <- function (...) {
  input <- lapply(list(...), as, "dsparseVector")
  thelength <- unique(sapply(input,length))
  stopifnot(length(thelength) == 1)
  return(sparseMatrix(
    x = unlist(lapply(input, slot, "x")),
    i = unlist(lapply(input, slot, "i")),
    p = c(0, cumsum(sapply(input, function(x) {length(x@x)}))),
    dims = c(thelength, length(input))
  ))
}
