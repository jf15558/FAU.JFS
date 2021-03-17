#' get_duplicates
#'
#' Function to detect elements with multiple higher
#' classifications in a tgraph object. Detection is
#' performed simply by calculating the per node degree
#' of each node in the graph representation of the
#' elements and their relationships. Degree > 1 indicates
#' duplicate classifications
#' @param x A tgraph object
#' @param ranks An optional numeric or character vector
#' indicating which ranks in x are to be checked for
#' duplicates
#' @return A tvertseq with the duplicate classifications
#' as focal elements
#' @import igraph

get_duplicates <- function(x, ranks = NULL) {

  # check for valid arguments
  if(!exists("x") | class(x) != "tgraph") {
    stop("Please supply a tgraph object")
  }
  if(is.null(ranks)) {
    ranks <- as.numeric(x$ranks)
  }
  if(is.numeric(ranks)) {
    if(any(isFALSE(ranks %in% x$ranks))) {
      stop("One or more of the supplied ranks is invalid")
    }
  } else {
    if(any(isFALSE(ranks %in% names(x$ranks)))) {
      stop("One or more of the supplied ranks is invalid")
    }
    ranks <- as.numeric(x$ranks[which(names(x$ranks) %in% ranks)])
  }

  # ensure that the tgraph object is in graph format
  if(x$type == "edgelist") {
    x <- tgraph_tedgelist(x)
  }
  # test for duplicate vertices and return tvertseq
  inv <- igraph::V(x$taxa)[which(V(x$taxa)$rank %in% ranks & igraph::V(x$taxa)$degree > 1)]
  out <- get_taxa(x, inv, mode = "parent")
  return(out)
}
