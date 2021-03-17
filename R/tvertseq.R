#' tvertseq
#'
#' Function to create a tvertseq object from scratch. The
#' tvertseq object displays the higher or lower relationships
#' for a given set of elements. These relationships may
#' intersect between elements and elements need not be at the
#' same level. The function can be used directly, but
#' requires more inputs compared to tvertex or tgraph and
#' exists primarily to complete the logical hierarchy between
#' these objects, and as a widely used carrier internally
#' within t* functions
#' @param x A four column dataframe consisting of two pairs of
#' columns. The first column specifies the focal elements in
#' the tvertseq object. The second column specifies elements
#' related, row by row, to the elements in column one. This
#' first pair of columns is effectively an edgelist. The third
#' and fourth columns are the respective ranks of the
#' corresponding elements in columns one and two
#' @param hrank A named vector of consective integers starting
#' from one. One is the highest rank in any t* object. The
#' names in hrank are the names of the actual rank system.
#' @param taxon The name of the focal element column in x
#' @param link The name of related elements column in x
#' @param ranks The name of the column denoting the ranks of
#' the focal elements in x
#' @param lranks The name of the column denoting the ranks of
#' the related elements in x
#' @param mode The direction of the relationships in the focal
#' and related elements column. If parent, then the
#' relationship elements exist at a higher rank than their
#' corresponding focal elements, and vice versa if child
#' @return A tvertseq object
#' @import igraph

tvertseq <- function(x, hrank, taxon = "taxon", link = "link", ranks = "rank", lranks = "lrank", mode = c("parent", "child")) {

  # check arguments
  if(!exists("x") | !is.data.frame(x)) {
    stop("Please supply a dataframe containing columns with taxa, their linked taxa, and the respective ranks of both")
  }
  if(!exists("hrank")) {
    stop("Please supply a vector of ranks in hierarchical order")
  }
  if(isFALSE(is.character(hrank))) {
    stop("hrank should be of mode character")
  }
  if(!(any(c(taxon, ranks) %in% colnames(x)))) {
    stop("The supplied column names must match those in the dataframe")
  }
  if(any(isFALSE(c(is.character(x[,taxon]), is.character(x[,lranks]),
                   is.character(x[,link]), is.character(x[,lranks]))))) {
    stop("The taxon and rank columns should be of mode character")
  }
  if(length(mode) > 1 | all(!mode %in% c("parent", "child"))) {
    stop("Mode should be specified as 'parent' or 'child'")
  }
  # invert if child
  if(mode == "child") {
    a <- taxon
    b <- link
    c <- ranks
    d <- lranks
    taxon <- b
    links <- a
    ranks <- d
    lranks <- c
  }

  # get table of edges
  cat(paste0("Filtering data"), "\n")
  # extract to dataframe in case of additional columns
  df <- x[,c(taxon, link, ranks, lranks)]
  # patch NA links to ORPHAN
  nas1 <- which(is.na(df[,link]))
  df[nas1,lranks] <- "ORPHAN"
  df[nas1,link] <- "ORPHAN"

  # ensure that hrank contains unique values
  hrank <- unique(hrank)
  # check that all ranks in the dataframe are represented in hrank
  inv1 <- which(!df[,ranks] %in% c(hrank, "ORPHAN"))
  inv2 <- which(!df[,lranks] %in% c(hrank, "ORPHAN"))
  inv <- c(inv1, inv2)
  if(length(inv) != 0) {
    m <- paste(c(unique(df[inv1,ranks], unique(df[inv2,lranks]))), collapse = ", ")
    paste0("Some taxa have ranks (", m, ") not present in hrank and will be added as informal vertices")
    warning(paste0("Some taxa have ranks (", m, ") not present in hrank and will be added as informal vertices"))
  }
  nas <- which(is.na(df[,taxon]))
  if(length(nas) > 0) {
    df <- df[!is.na(df[,taxon]),]
    warning("Some taxa had NA values and have been dropped")
  }
  # set up hierarchy
  fr <- 0:length(hrank)
  names(fr) <- c("ORPHAN", hrank)

  # get table of edges
  cat(paste0("Generating table of edges"), "\n")
  # add to dataframe
  df$r1 <- fr[match(df[,ranks], names(fr))]
  df$r2 <- fr[match(df[,lranks], names(fr))]
  # get the frequencies of links before trimming
  cat(paste0("Tabulating edge frequencies"), "\n")
  df_n <- paste(df[,taxon], df[,link], sep = "|")
  df_n2 <- table(df_n)
  df$freq <- df_n2[match(df_n, names(df_n2))]
  df <- unique(df)
  d_name <- cbind.data.frame(c(df[,taxon], df[,link]), c(df[,ranks], df[,lranks]),c(df[,"r1"], df[,"r2"]))
  d_name <- unique(d_name)
  # check for name re-use at multiple ranks
  test <- table(d_name[,1])
  if(max(test) > 1) {
    inv2 <- paste0(names(test[which(test > 1)]), collapse = ", ")
    stop(paste0("The following taxa has multiple ranks: ", inv2))
  }

  # make the graph
  cat(paste0("Building tvertseq"), "\n")
  # create graph, simplifying to remove loops from NA-filled rows
  g <- igraph::graph_from_edgelist(as.matrix(df[,c(taxon, link)]))
  igraph::E(g)$freq <- df$freq
  if("ORPHAN" %in% df[,1] | "ORPHAN" %in% df[,2]) {
    g <- igraph::delete.vertices(g, "ORPHAN")
  }
  igraph::V(g)$rank <- d_name[match(igraph::V(g)$name, d_name[,1]),3]
  # account for informal ranks if needed
  igraph::V(g)$inf_rank <- NA
  igraph::V(g)$inf_rank[is.na(igraph::V(g)$rank)] <- d_name[match(igraph::V(g)$name[is.na(igraph::V(g)$rank)], d_name[,1]),2]
  igraph::V(g)$degree <- igraph::degree(g, mode = "out")

  # get subgraph for focal taxa
  g2 <- igraph::make_graph(rep(df[,taxon], each = 2))
  g2 <- igraph::simplify(g2)
  igraph::V(g2)$inf_rank <- igraph::V(g)$inf_rank[match(igraph::V(g2)$name, igraph::V(g)$name)]
  igraph::V(g2)$rank <- igraph::V(g)$rank[match(igraph::V(g2)$name, igraph::V(g)$name)]
  igraph::V(g2)$degree <- igraph::V(g)$degree[match(igraph::V(g2)$name, igraph::V(g)$name)]

  # make tvertseq
  out <- list()
  out[[1]] <- g2
  out[[2]] <- fr[-1]
  out[[3]] <- g
  out[[4]] <- mode
  names(out) <- c("taxa", "ranks", "seq", "mode")
  class(out) <- "tvertseq"
  return(out)
}
