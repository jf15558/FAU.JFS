#' as_tgraph
#'
#' Function to convert a tvertseq object to a tgraph object
#' This function simply takes the graph portion of the
#' tvertseq as the graph for the tgraph, nullifying the
#' relevance of parent or child relationships
#' @param x A tvertseq object
#' @return A tgraph object
#' @import igraph

as_tgraph <- function(x) {

  if(!exists("x")) {
    stop("Please provide an object to convert")
  }
  if(class(x) != "tvertseq") {
    stop("x must be a tvertseq object")
  }
  out <- list()
  out[[1]] <- x$seq
  out[[2]] <- x$ranks
  out[[3]] <- "graph"
  out[[4]] <- "Resolved"
  names(out) <- c("taxa", "ranks", "type", "status")
  igraph::V(out$taxa)$degree <- igraph::degree(graph = out$taxa, mode = "out")
  if(max(unique(igraph::V(out[[1]])$degree)) > 1) {
    out$status <- "Unresolved"
  }
  class(out) <- "tgraph"
  return(out)
}
