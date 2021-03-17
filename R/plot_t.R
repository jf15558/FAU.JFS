#' plot_t
#'
#' Function to plot the parent or child relationships of an
#' element in a tgraph or tvertseq object. The length
#' of edges beyond which relationships will not be plotted
#' can be set in the function arguments
#' @param x a tgraph or tvertex object
#' @param taxa A character vector of element names whose
#' relationships will be plotted
#' @param taxa2 An optional numeric vector of element names
#' whose relationships will be plotted
#' @param mode The direction of the relationships to be
#' plotted
#' @param step A positive integer specifinyg the
#' neighbourhood of the relationships to plot. Specifying
#' a number greater than the number of ranks in the t* object
#' will not cause a failure, and will instead plot all
#' relationships in the direction specified in mode
#' @return A plot of the relationships of the specified
#' elements
#' @import igraph

plot_t <- function(x, taxa = NULL, taxa2 = NULL, mode = c("parent", "child"), step = NULL) {

  # check arguments
  if(!(class(x) %in% c("tgraph", "tvertseq"))) {
    stop("Please supply either a tgraph object and taxon name, or a tvertseq object")
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

  # retrieve taxa specified if a tvertex and taxon is supplied
  if(class(x) == "tvertseq") {
    if(!is.null(taxa)) {
      if(!is.null(taxa2)) {
        taxa2 <- igraph::V(x$seq)[taxa2]$name
        taxa <- unique(c(taxa, taxa2))
      }
      x <- as_tgraph(x)
      x <- get_taxa(x, taxa = taxa, mode = mode, step = step)
    } else {
      x <- x
    }
  }

  # retrieve tvertseq object if a tgraph and taxon is supplied
  if(class(x) == "tgraph") {
    if(is.null(taxa)) {
      stop("If a tgraph is supplied, then taxon must also be specified")
    }
    # if taxa2 is additionally specified, add the vertex names into taxa
    if(!is.null(taxa2)) {
      taxa2 <- igraph::V(x$taxa)[taxa2]$name
      taxa <- unique(c(taxa, taxa2))
    }
    x <- get_taxa(x, taxa = taxa, mode = mode, step = step)
  }

  # plot the retrieved tvertseq object (works for single or multiple taxa)
  lay <- igraph::layout_with_sugiyama(x$seq, layers = as.numeric(igraph::V(x$seq)$rank), hgap = 1)
  rseq <- seq(from = -1, to = 1, length.out = length(x$ranks))
  par(mar = c(0, 5, 0, 5))
  par(las = 1)
  plot(x$seq, layout = lay$layout, edge.arrow.mode = 0)
  axis(4, labels = rev(names(x$ranks)), at = rseq, col = NA)
}
