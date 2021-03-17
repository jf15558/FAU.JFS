#' as_tvertseq
#'
#' Function to coerce a tgraph object to a tvertseq object.
#' Currently, only the highest ranked or lowest ranked
#' elements can be coerced to focal elements in a tvertseq,
#' depending on the argument specification, although
#' future versions may be adapted to allow any single rank
#' to be made focal, with parent or child relationships
#' retained in the tvertseq
#' @param x A tgraph object
#' @param mode A character vector of length one specifying
#' whether parent or child relationships will be recorded
#' in the tvertseq. As only highest or lowest rank
#' conversion is currently supported, mode 'parent' will
#' result in conversion of the lowest rank, while mode
#' 'child' will result in conversion of the highest rank
#' @return A tvertseq object
#' @import igraph

as_tvertseq <- function(x, mode = c("parent", "child")) {

  if(!exists("x")) {
    stop("Please provide an object to convert")
  }
  if(class(x) != "tgraph") {
    stop("x must be a tgraph object")
  }
  # default behaviour is to retrieve parents
  if(length(mode) > 1 | !mode %in% c("parent", "child")) {
    stop("Mode should be specified as 'parent' or 'child'")
  }
  if(mode == "parent") {
    pos <- max(x$ranks)
  } else {
    pos <- min(x$ranks)
  }
  out <- list()
  out[[1]] <- igraph::induced_subgraph(x$taxa, which(igraph::V(x$taxa)$rank == pos))
  out[[1]] <- igraph::delete.edges(out[[1]], edges = igraph::E(out[[1]]))
  out[[2]] <- x$ranks
  out[[3]] <- x$taxa
  out[[4]] <- mode
  names(out) <- c("taxa", "ranks", "seq", "mode")
  class(out) <- "tvertseq"
  return(out)
}
