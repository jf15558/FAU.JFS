#' tgraph_tedgelist
#'
#' Function to convert the relationships in a tgraph
#' from a graph format to an edgelist. The current
#' format of the supplied object is detected and
#' automatically converted to the other format
#' @param from a tgraph object
#' @return a tgraph object with the alternative
#' formatting of the element relationships in from
#' @import igraph

tgraph_tedgelist <- function(from) {

  if(!exists("from")) {
    stop("Please supply an tgraph object")
  }
  if(class(from) != "tgraph") {
    stop("Please supply an tgraph object")
  }

  if(from$type == "graph") {
    x <- from$taxa
    # extract vertex attributes
    ats <- get_tattr(x = from)
    # coerce graph to edgelist and add in the ranks
    g_list <- as.data.frame(igraph::as_edgelist(x))
    colnames(g_list) <- c("taxon", "parent")
    # grab orphan vertices (no edges so absent above)
    g_0 <- names(igraph::V(x))[!(names(igraph::V(x)) %in% g_list$taxon)]
    g0_r <- igraph::V(x)[g_0]$rank
    g_list2 <- cbind.data.frame(g_0, g_0)
    colnames(g_list2) <- colnames(g_list)
    # bind in the orphan vertices
    g_list <- rbind.data.frame(g_list, g_list2)
    # bind in vertex attributes
    for(i in 2:ncol(ats)) {
      vec <- ats[match(g_list[,1], ats$name), i]
      g_list <- cbind(g_list, vec)
    }
    colnames(g_list) <- c("taxon", "parent", colnames(ats)[2:ncol(ats)])
    # return
    from$taxa <- g_list
    from$type <- "edgelist"
    return(from)
  }

  if(from$type == "edgelist") {
    x <- from$taxa
    # convert the edgelist to a graph
    t_graph_int <- igraph::graph_from_edgelist(as.matrix(x[,1:2]))
    t_graph_int <- igraph::simplify(t_graph_int, remove.multiple = FALSE)
    from$taxa <- t_graph_int
    from$type <- "graph"
    # add in attributes, redoing degree checking to be safe
    ats <- x[,c(1, 3:ncol(x))]
    ats$degree <- igraph::degree(from$taxa, igraph::V(from$taxa)[x$taxon], mode = "out")
    from <- add_tattr(tob = from, x = ats, attr = NULL)
    return(from)
  }
}
