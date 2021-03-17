#' pbdb_kingdoms
#' 
#' Helper list of the kingdom-level assignment of
#' phylum names in the PBDB. Used to build a kingdom
#' column in a PBDB download if taxonomy fields were
#' requested. List of three character vectors, named
#' animals, plants, protists. There is only a single
#' fungal phylum in the database so this is assigned
#' directly, while all other phyla are assigned to
#' bacteria
#' 
#' @docType data
#' @usage data(pbdb_kingdoms)
#' @keywords datasets
"pbdb_kingdoms"