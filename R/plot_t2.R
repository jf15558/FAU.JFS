#' plot_t2
#'
#' Function to plot the parent or child relationships of an
#' element in a hierarchically organised dataframe. An
#' accompaniment to @seealso plot_t which internally coerces
#' the input into a format that would be plotted as in
#' plot_t. The length of edges beyond which relationships
#' will not be plotted can be set in the function arguments
#' @param x a dataframe containing hierarchically organised
#' data in columns
#' @param taxon A character vector of element names whose
#' relationships will be plotted (these must be of the same
#' rank)
#' @param trank A character vector of length one corresponding
#' to the column name in x in which taxa is located
#' @param ranks A character vector corresponding to the column
#' names in x, given in hierarchical order
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
#' @export

plot_t2 <- function(x, taxon, trank, ranks, mode = c("parent", "child", "all"), step = NULL) {

  if(length(mode) == 3) {
    stop("Please specify the plotting mode")
  }

  tx <- unique(x[which(x[,trank] == taxon), ranks])
  tg <- tgraph(tx)
  if(mode != "all") {
    tv <- get_taxa(tg, taxa = taxon, mode = mode)
    lay <- igraph::layout_with_sugiyama(tv$seq, layers = as.numeric(igraph::V(tv$seq)$rank), hgap = 1)
    rseq <- seq(from = -1, to = 1, length.out = length(unique(igraph::V(tv$seq)$rank)))
    par(mar = c(0, 5, 0, 5))
    par(las = 1)
    plot(tv$seq, layout = lay$layout, edge.arrow.mode = 0)
    axis(4, labels = rev(names(tv$ranks)[unique(igraph::V(tv$seq)$rank[order(igraph::V(tv$seq)$rank)])]), at = rseq, col = NA)
  } else {
    tv <- tg
    lay <- igraph::layout_with_sugiyama(tv$taxa, layers = as.numeric(igraph::V(tv$taxa)$rank), hgap = 1)
    rseq <- seq(from = -1, to = 1, length.out = length(tv$ranks))
    par(mar = c(0, 5, 0, 5))
    par(las = 1)
    plot(tv$taxa, layout = lay$layout, edge.arrow.mode = 0)
    axis(4, labels = rev(names(tv$ranks)[unique(igraph::V(tv$taxa)$rank[order(igraph::V(tv$taxa)$rank)])]), at = rseq, col = NA)
  }
}
