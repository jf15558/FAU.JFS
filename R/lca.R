#' lca
#'
#' Function to retrieve the last common ancestor of a set
#' of elements in a tgraph object
#' @param graph A tgraph object
#' @param ... A vector of vertex names or numbers for which
#' the common ancestor will be derived
#' @return The name of the common element
#' @import igraph

lca <- function(graph, ...) {
  dots = c(...)
  path = igraph::ego(graph, order = length(igraph::V(graph)), nodes = dots, mode = "out")
  max(Reduce(intersect, path))
}
