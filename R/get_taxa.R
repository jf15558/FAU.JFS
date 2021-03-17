#' get_taxa
#'
#' Function to retrive the parent or child relationships
#' for a set of elements in a tgraph object. The length
#' of edges beyond which relationships will not be returned
#' can be set in the function arguments
#' @param x A tgraph object
#' @param taxa The names of the elements whose relationships
#' will be retrived, or the numeric ids of nodes in the
#' graph portion of the tgraph object (useful for
#' retrieving relationships where the name name of a node
#' may not be known)
#' @param mode The directionality of the relationships
#' to retrieve
#' @param step A positive integer specifinyg the
#' neighbourhood of the relationships to retrive. Specifying
#' a number greater than the number of ranks in the tgraph
#' will not cause a failure, and will instead retrieve
#' all relationships in the direction specified in mode
#' @return a tvertseq object containing the focal elment
#' and its relationships, as specified by mode
#' @import igraph

get_taxa <- function(x, taxa, mode = c("parent", "child"), step = NULL) {

  # check for valid arguments
  if(!exists("x") | class(x) != "tgraph") {
    stop("Please supply a tgraph object")
  }
  if(!exists("taxa")) {
    stop("Please supply a taxon or taxa to retrive")
  }
  if(length(mode) > 1 | !mode %in% c("parent", "child")) {
    stop("Mode should be specified as 'parent' or 'child'")
  }
  if(mode == "parent") {
    mode2 <- "out"
  } else {
    mode2 <- "in"
  }
  if(is.null(step)) {
    n <- length(x$ranks)
  }
  if(length(taxa) == 0) {
    out <- list()
    out[[1]] <- taxa
    out[[2]] <- x$ranks
    out[[3]] <- NA
    out[[4]] <- NA
  } else {
    if(class(taxa) == "igraph.vs") {
      taxa <- names(taxa)
    }
    if(is.character(taxa)) {
      taxa <- which(igraph::V(x$taxa)$name %in% taxa)
    }

    # get paths for each taxon
    try <- igraph::ego(x$taxa, order = n, nodes = taxa, mode = mode2)
    # combine into single graph and get the unique set of vertices
    try <- lapply(try, function(y) {igraph::induced_subgraph(x$taxa, y)})
    try <- igraph::V(Reduce(igraph::union, try))$name
    # get the subgraph for those vertices
    try <- igraph::induced_subgraph(x$taxa, try)

    # get subgraph for focal taxa
    g2 <- igraph::make_graph(rep(as.character(igraph::V(x$taxa)[taxa]$name), each = 2))
    g2 <- igraph::simplify(g2)

    # transfer vertex attributes
    ca <- do.call(cbind.data.frame, igraph::get.vertex.attribute(x$taxa))
    ca <- ca[!colnames(ca) %in% "name"]
    if(length(ca) != 0) {
      for(i in 1:ncol(ca)) {
        g2 <- igraph::set.vertex.attribute(g2, name = colnames(ca)[i], index = igraph::V(g2), value = ca[match(igraph::V(g2)$name, rownames(ca)),i])
      }
    }

    # create output object
    out <- list()
    out[[1]] <- g2
    out[[2]] <- x$ranks
    out[[3]] <- try
    out[[4]] <- mode
  }
  names(out) <- c("taxa", "ranks", "seq", "mode")
  class(out) <- "tvertseq"
  return(out)
}
