#' get_tattr
#'
#' Function to retrieve the vertex attributes of a t*
#' object
#' @param x A t* object
#' @return A dataframe of the vertex names in x, along
#' with their attributes
#' @import igraph
#' @importFrom utils head
#' @importFrom stats na.omit
#' @importFrom dplyr bind_rows
#' @importFrom data.table setDF .SD

get_tattr <- function(x) {

  if(!exists("x")) {
    stop("Please supply a tgraph, tvertex or tvertex object, or a list of such objects")
  }
  if(typeof(x) != "list") {
    stop("Please supply a tgraph, tvertex or tvertex object, or a list of such objects")
  }
  if(class(x) %in% c("tgraph", "tvertseq", "tvertex")) {
    x <- list(x)
  } else {
    chk <- unlist(lapply(x, class))
    if(!all(chk %in% c("tgraph", "tvertseq", "tvertex"))) {
      stop("Please supply a tgraph, tvertex or tvertex object, or a list of such objects")
    }
  }

  # extract all attributes and bind into dataframe
  x <- lapply(x, function(x) {do.call(cbind.data.frame, igraph::get.vertex.attribute(x$taxa))})
  df1 <- do.call(dplyr::bind_rows, x)

  # remove conflicting values, then add back any rows which contained all NA values
  dfn <- as.data.frame(data.table::setDT(df1)[, lapply(data.table::.SD, function(x) head(na.omit(x), 1L)), by = "name"])
  df <- rbind(dfn, df1[!df1$name %in% dfn$name,])

  # return
  return(df)
}
