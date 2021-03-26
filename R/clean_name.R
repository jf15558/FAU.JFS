#' roxygen documentation
#'
#' clean_name
#'
#' Function which bundles a series of cleaning routines into a
#' single process. First any words in brackets are removed,
#' followed by a series of user-defined terms if given. Next
#' Roman and Arabic numerical are removed, then abbreviations
#' up to five letters (abbreviations are matched by the
#' following dot e.g ABFS.). By default, characters for removal
#' are replaced by a white space to prevent accidental collapse
#' of strings. However, there may be specific cases where a
#' collapse is required and so terms given in collapse are
#' dealt with next. After collapsing, rogue all rogue punctation
#' is removed, then isolated lowercase letters, then isolated
#' groups of capitals up to 5 characters long. Finally, white
#' spaces greater than 1 are removed, along with trailing white
#' space, any remaining strings longer than 2 words subsetted
#' to the first word, the first letter of each string capitalised
#' and zero length strings converted to NA
#' @param x a character vector of names to clean
#' @param terms a character vector of terms to remove from
#' elements of x. Terms are only removed as whole words, rather
#' than if they also happen to occur as strings within elements
#' of x
#' @param collapse a character vector of strings which should
#' collapsed (i.e. replaced by "", rather than the default " ").
#' If one of the collapse terms is a special regex character, it
#' will need to be escaped, e.g. "\\-"
#' @return a character vector the same length as x. Elements
#' which were reduced to zero characters during cleaning are
#' returned as NA
#' @import stringr
#' @export

clean_name <- function(x, terms = NULL, collapse = NULL) {

  # set up important objects
  id_clean <- x
  cat(paste0("Beginning cleaning"), "\n")

  # remove bracketed terms
  id_clean <- gsub("(\\(+.+\\))", "", id_clean)

  # remove terms if specified
  if(!is.null(terms)) {
    id_clean <- gsub(paste0("(?i)\\b", terms, "+\\b", collapse = "|"), " ", id_clean)
  }
  cat(paste0("Term cleaning done"), "\n")

  # Roman and arabic numerals
  id_clean <- gsub("\\b[XVI]+\\b|[0-9]", " ", id_clean)
  # clear abbreviations (those which are formatted correctly with a period, up to five characters e.g. sp. or indet.) {,5}
  id_clean <- gsub("(?i)[a-z]{,20}\\.", " ", id_clean)
  # collapse any terms
  if(!is.null(collapse)) {
    collapse <- paste0(collapse, collapse = "|")
    id_clean <- gsub(collapse, "", id_clean)
  }
  # any other punctuation
  id_clean <- gsub("[[:punct:]]", " ", id_clean)
  # trim isolated letter groups (up to 5 letters for capitals e.g species CCCDE, 1 for lowercase e.g. sp. a)
  id_clean <- gsub("\\b[A-Z]{1,5}\\b|\\b[a-z]{1}\\b", " ", id_clean)
  cat(paste0("Regular cleaning done"), "\n")

  # resolve whitespaces
  id_clean <- gsub("\\s+", " ", id_clean)
  id_clean <- gsub("^ | $", "", id_clean)
  # for the remaining 3+ word ids, take the first word
  id_clean[sapply(strsplit(id_clean, " "), length) > 2] <- stringr::word(id_clean[sapply(strsplit(id_clean, " "), length) > 2])
  # ensure capitalisation of the first letter of the first word (PBDB standard) and lowercase for all other characters
  id_clean[!is.na(id_clean)] <- tolower(id_clean[!is.na(id_clean)])
  id_clean[!is.na(id_clean)] <- paste0(toupper(substr(id_clean[!is.na(id_clean)], 1, 1)), substr(id_clean[!is.na(id_clean)], 2, nchar(id_clean[!is.na(id_clean)])))
  # zero length strings to NA
  id_clean[which(nchar(id_clean) == 0)] <- NA
  cat(paste0("Final formatting done"), "\n")

  # return cleaned name
  return(id_clean)
}
