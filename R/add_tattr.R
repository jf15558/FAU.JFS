#' add_tattr
#'
#' Function to add attributes to the vertices of a
#' t* object
#' @param tob A t* object to add attributes to
#' @param x a vector of element names in tob. Alternatively
#' a dataframe with element names and attribute values as its
#' columns
#' @param xattr a vector comprising the attributes
#' of the elements in x. This must be the same
#' length as x
#' @param attr If x is a dataframe, a vector specifying the
#' attribute value columns in x
#' @param taxon If x is a dataframe, a vector specifying the
#' taxon name columns in x
#' @param tadd If x is not specified, tadd can be supplied as a t*
#' object, whose attributes will be added to the corresponding
#' element names in tob
#' @return An object of the same class as tob, with updated
#' attributes
#' @import igraph

add_tattr <- function(tob, x = NULL, xattr = NULL, attr = "attr", taxon = "taxon", tadd = NULL) {

  # check arguments
  if(!exists("tob")) {
    stop("Please supply a tgraph, tvertseq or tvertex to add attributes to")
  }
  if(!is.null(x) & !is.null(tadd)) {
    stop("x and tadd cannot be specified simultaneously")
  }
  if(!exists("x") & is.null(tadd)) {
    stop("Please supply vectors of names and their attributes or an equivalent dataframe. Alternatively,
    specify a t* object or list of the same to extract attributes from (tadd)")
  }
  if(any(attr %in% c("rank", "inf_rank", "degree"))) {
    warning("The attribute names rank, inf_rank and degree are reserved.
            Other attribute names should be used")
  }
  if(is.null(x) & !is.null(tadd)) {
    x <- get_tattr(tadd)
    taxon <- "name"
    attr <- colnames(x)[-1]
  }
  if(all(!is.vector(x), !is.data.frame(x))) {
    stop("Please supply vectors of names and their attributes or an equivalent dataframe. Alternatively,
    specify a t* object or list of the same to extract attributes from (tadd)")
  }

  # coerce to dataframe, ensuring both columns are of mode character
  if(is.data.frame(x)) {
    if(is.null(attr)) {
      attr <- as.character(colnames(x))
    }
    if(!(any(c(taxon, attr) %in% colnames(x)))) {
      stop("The supplied column names must match those in the dataframe")
    }
    attr <- attr[!attr %in% taxon]
    df <- x[,c(taxon, attr)]
    colnames(df) <- c("name", attr)
  }
  if(is.vector(x)) {
    if(is.null(xattr)) {
      stop("If x is a vector, then their attributes must also be provided")
    }
    if(length(x) != length(xattr)) {
      stop("x and xattr should be of the same length")
    }
    df <- cbind.data.frame(x, xattr)
    colnames(df) <- c("name", attr)
  }
  attr <- colnames(df)

  # just in case, drop any NA taxon names if present (NA attributes are fine)
  df <- df[!is.na(df[,"name"]),]
  # trim to unique combinations
  df <- unique(df)
  # check that multiple states are not given
  inv <- table(df[,"name"])
  if(max(inv) > 1) {
    stop("There there can only be one value of each attribute per taxon")
  }

  # done from second column as the first column is the name column and should not be overwritten
  if(class(tob) == "tgraph") {
    for(i in 2:length(attr)) {
      tob$taxa <- igraph::set.vertex.attribute(tob$taxa, name = attr[i], index = igraph::V(tob$taxa),
                                               value = df[match(igraph::V(tob$taxa)$name, df[,"name"]),attr[i]])
    }
  }
  if(class(tob) == "tvertseq") {
    for(i in 2:length(attr)) {
      tob$seq <- igraph::set.vertex.attribute(tob$seq, name = attr[i], index = igraph::V(tob$seq),
                                              value = df[match(igraph::V(tob$seq)$name, df[,"name"]),attr[i]])
    }
    tob$taxa <- igraph::delete.vertices(tob$seq, v = igraph::V(tob$seq)$name[!igraph::V(tob$seq)$name %in% igraph::V(tob$taxa)$name])
    tob$taxa <- igraph::delete.edges(tob$taxa, edges = igraph::E(tob$taxa))
  }
  if(class(tob) == "tvertex") {
    for(i in 2:length(attr)) {
      tob$taxa <- igraph::set.vertex.attribute(tob$taxa, name = attr[i], index = igraph::V(tob$taxa),
                                               value = df[match(igraph::V(tob$taxa)$name, df[,"name"]),attr[i]])
    }
  }
  return(tob)
}
