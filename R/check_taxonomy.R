#' check_taxonomy
#'
#' Function to implement a multi-step cleaning routine for
#' hierarchically structured taxonomic data. The first
#' part of the routine calls @seealso spell_check to flag
#' spelling errors within a given taxonomic group. The second
#' part of the routine calls @seealso check_ranks to flag
#' name re-use at different levels. Some of these cases may
#' arise when a name has been unfortunately, but validly
#' used to refer to different groups at different taxonomic
#' levels. The third part of the routine calls @seealso
#' find_duplicates to flag variable higher classifications for
#' a given taxon. The function flags errors by default, but can
#' return a fully cleaned taxonomic scheme with conflicting
#' classifications eliminated or split by the addition of suffixes
#' to give unique taxonomic names. In the case of non-rank unique
#' names, the suffixes are of the form _RNK + the column number,
#' e.g _RNK1. For conflicting classifications, cleaning is
#' performed using @seealso resolve_duplicates, and in the case
#' of re-use of a name at the same rank for genuinely different
#' taxa, suffixes are capital letters, e.g. AlloniaA, AlloniaB
#' @param x A dataframe with hierarchically organised
#' taxonomic information. If x only comprises the taxonomic
#' information, @param ranks does not need to be specified, but the
#' columns must be in order of decreasing taxonomic rank
#' @param ranks The column names of the taxonomic data fields
#' in x. These must be provided in order of decreasing taxonomic
#' rank
#' @param routine A character vector determining the flagging
#' and cleaning routines to employ. Valid values are spell_check
#' (flag potential spelling errors), discrete_rank (check that
#' taxonomic names are unique to their rank), duplicate_tax (flag
#' conflicting higher classifications of a given taxon)
#' @param jw a numeric greater than 0 and less than 1. This is
#' the distance threshold below which potential synonyms will be
#' considered. Called by @seealso spell_check
#' @param str A positive integer specifying the
#' number of matching characters at the beginning of synonym
#' pairs. By default 1, i.e. the first letters must match. Called
#' by @seealso spell_check
#' @param str2 If not NULL, a positive integer specifying the
#' number of matching characters at the end of synonym pairs.
#' Called by @seealso spell_check
#' @param method2 A character string of length one corresponding
#' to one of the methods used by @seealso afind. One of "osa",
#' "lv", "dl", "hamming", "lcs", "qgram", "cosine",
#' "running_cosine", "jaccard", or "soundex". Called by @seealso
#' spell_check
#' @param q q-gram size. Only used when method2 is "qgram",
#' "cosine" or "Jaccard". Called by @seealso spell_check
#' @param pref If not NULL, a character vector of prefixes (which
#' will be used at all ranks) or a list of rank-specific prefixes,
#' which may result in erroneously low JW distances. Synonyms will only
#' be considered if both terms share the same prefix. If a list, this
#' should be given in descending rank order. Called by @seealso
#' spell_check
#' @param suff If not NULL, a character vector of prefixes (which
#' will be used at all ranks) or a list of rank-specific prefixes,
#' which may result in erroneously low JW distances. Synonyms will only
#' be considered if both terms share the same suffix. If a list, this
#' should be given in descending rank order. Called by @seealso
#' spell_check
#' @param exclude If not NULL, a character vector of group names (which
#' will be used at all ranks) or a list of rank-specific group names,
#' which should be skipped - useful for groups which are known
#' to contain potentially similar terms. If a list, this
#' should be given in descending rank order. Called by @seealso
#' spell_check
#' @param jump The maximum number of ranks between the rank of
#' classification conflict and the next common classification (if present),
#' below which the divergence will be taken as conflicting.
#' Called by @seealso resolve_duplicates
#' @param plot A logical speciying if conflicting classifications should
#' be plotted. Called by @seealso resolve_duplicates
#' @param spell_clean If TRUE, the function will return a cleaned version
#' of the supplied taxonomic dataframe, using the supplied threshold
#' for the similarity method given by method2, to automatically update
#' any names in pairs of flagged synonyms to the more frequent spelling.
#' This is not recommended, however, so the argument is FALSE by default
#' and the threshold left as NULL
#' @param thresh The threshold for the similarity method given by method2,
#' below which flagged pairs of names will be considered synonyms and
#' resolved automatically. See @seealso spell_check for details on method2
#' @param tax_clean If TRUE, the function will return a cleaned version
#' of the supplied taxonomic dataframe, using @seealso resolve_duplicates
#' to resolve conflicts in the way documented by the function. Both
#' spell_clean and tax_clean can both be TRUE to return a dataset cleaned
#' by both methods
#' @param append If TRUE, any suffixes used during cleaning will
#' be retained in the cleaned version of the data. This is preferable
#' as it ensures that all taxonomic names are rank-discrete and
#' uniquely classified
#' @return A list with elements corresponding to the outputs of the
#' chosen flagging routines (three by default), plus
#' a cleaned taxonomic dataframe if either spell_clean or tax_clean or both
#' are TRUE. See @seealso spell_clean, @seealso check_ranks and @seealso
#' find_duplicates for details of the structure of the outputs
#' @importFrom stats na.omit
#' @export

check_taxonomy <- function(x, ranks = NULL, routine = c("spell_check", "discrete_rank", "duplicate_tax"),
                           jw = 0.1, str = 1, str2 = NULL, method2 = "jaccard", q = 1, pref = NULL, suff = NULL, exclude = NULL,
                           jump = 3, plot = FALSE,
                           spell_clean = FALSE, thresh = NULL, tax_clean = TRUE, append = TRUE) {

  # check that data has minimally been supplied
  if(!exists("x")) {
    stop("Please supply x as a dataframe of taxonomic assignments")
  }
  # coerce to dataframe with column names to be safe
  if(!is.data.frame(x)) {x <- as.data.frame(x)}
  if(is.null(colnames(x))) {colnames(x) <- as.character(1:ncol(x))}

  # check that ranks are column names of x
  if(is.null(ranks)) {ranks <- colnames(x)}
  if(!all(ranks %in% colnames(x))) {
    stop("Not all elements of argument ranks are column names in x")
  }
  # check that ranks are in hierarchical order
  if(length(unique(x[,ranks[length(ranks)]])) < length(unique(x[,ranks[(length(ranks) - 1)]]))) {
    warning("Higher taxonomy is more diverse than lower taxonomy. Are the columns in x
            or the column names specified in 'ranks' supplied in descending hierarchical order?")
  }
  # subset to columns
  x <- x[,ranks]

  # check that the data is character
  if(!all(apply(x, 2, class) == "character")) {
    stop("Not all columns in x are of class character")
  }

  # check cleaning routines have been correctly supplied
  if(!is.character(routine)) {
    stop("Routine should be a character vector containing one or more of the following:
         spell_check, discrete_rank, report_tax, resolve_tax")
  }
  if(!any(routine %in% c("spell_check", "discrete_rank", "duplicate_tax"))) {
    stop("All elements of argument routine are invalid.
         Valid elements are spell_check", "discrete_rank", "duplicate_tax")
  }
  if(!all(routine %in% c("spell_check", "discrete_rank", "duplicate_tax"))) {
    warning("Some elements of argument routine are invalid and will be ignored.
            Valid elements are spell_check", "discrete_rank", "duplicate_tax")
  }
  # ensure routine vector is clean and correctly ordered
  routine <- unique(routine)
  routine <- as.vector(na.omit(routine[match(c("spell_check", "discrete_rank", "duplicate_tax"), routine)]))

  # check additional flags
  if("spell_check" %in% routine) {

    # check any supplied prefixes
    if(is.null(pref)) {
      pref <- as.list(1:ncol(x))
      pref_list <- lapply(pref, function(x) {x <- NULL})
    } else {
      if(is.atomic(pref)) {
        pref_list <- lapply(1:ncol(x), function(x) {x <- pref})
      }
      if(is.list(pref)) {pref_list <- pref}
      if(length(pref_list) != length(ranks)) {
        stop("Prefixes should be supplied either as a vector which will be used at all taxonomic levels, or as a list of
             vectors to be used at each specific level in x (length must equal number of columns in x/number of ranks")
      }
      if(!all(unlist(lapply(pref_list, class))) == "character") {
        stop("Not all elements of pref are of class character")
      }
      pref_list <- lapply(pref_list, function(x) {as.vector(na.omit(x))})
    }

    # check any supplied suffixes
    if(is.null(suff)) {
      suff <- as.list(1:ncol(x))
      suff_list <- lapply(suff, function(x) {x <- NULL})
    } else {
      if(is.atomic(suff)) {
        suff_list <- lapply(1:ncol(x), function(x) {x <- suff})
      }
      if(is.list(suff)) {suff_list <- suff}
      if(length(suff_list) != length(ranks)) {
        stop("Suffixes should be supplied either as a vector which will be used at all taxonomic levels, or as a list of
             vectors to be used at each specific level in x (length must equal number of columns in x/number of ranks")
      }
      if(!all(unlist(lapply(suff_list, class))) == "character") {
        stop("Not all elements of pref are of class character")
      }
      suff_list <- lapply(suff_list, function(x) {as.vector(na.omit(x))})
    }

    # check any supplied exclusions
    if(is.null(exc)) {
      exc <- as.list(1:ncol(x))
      exc_list <- lapply(exc, function(x) {x <- NULL})
    } else {
      if(is.atomic(exc)) {
        exc_list <- lapply(1:ncol(x), function(x) {x <- exc})
      }
      if(is.list(exc)) {exc_list <- exc}
      if(length(exc_list) != length(ranks)) {
        stop("Exclusions should be supplied either as a vector which will be used at all taxonomic levels, or as a list of
             vectors to be used at each specific level in x (length must equal number of columns in x/number of ranks")
      }
      if(!all(unlist(lapply(exc_list, class))) == "character") {
        stop("Not all elements of pref are of class character")
      }
      exc_list <- lapply(exc_list, function(x) {as.vector(na.omit(x))})
    }

    if(spell_clean) {
      warning("As spell checking is approximate, automatic resolution is not recommended. Any flagged names should be properly checked")
      if(is.null(thresh)) {
        stop("If spell_clean has been requested, thresh must be specified. This threshold is specific to method2, see spell_check documentation")
      }
      if(!is.numeric(thresh) | any(thresh < 0)) {
        stop("Thresh must be a positive numeric to be used at all ranks, or vector of values to be used individually at each rank. See spell_check for details on value choice")
      }
      if(length(thresh == 1)) {
        thresh <- rep(thresh, length(ranks))
      }
      if(length(thresh != length(ranks))) {
        stop("Thresh must be a positive numeric to be used at all ranks, or vector of values to be used individually at each rank. See spell_check for details on value choice")
      }
    }
  }

  # set up output list
  out <- list()

  # check spelling
  if("spell_check" %in% routine) {

    spell_list <- list()
    for(i in 1:(length(ranks) - 1)) {
      spell_list[[i]] <- spell_check(x = x, terms = ranks[i + 1], groups = ranks[i],
                                     jw = jw, str = str, str2 = str2, method2 = method2,
                                     q = q, pref = pref_list[[i + 1]], suff = suff_list[[i + 1]], exclude = exc_list[[i + 1]])
      spell_list[[i]] <- cbind.data.frame(level = rep(ranks[i + 1], nrow(spell_list[[i]])), spell_list[[i]])
      message(paste0(nrow(spell_list[[i]]), " potential synonyms flagged at the ", ranks[i + 1], " level"))
    }
    out[[1]] <- do.call(rbind, spell_list)
    if(spell_clean) {
      ob <- out[[1]]
      for(j in 1:nrow(ob)) {
        if(ob$m2[j] < thresh[i]) {
          x[which(x[,ob$level[j]] == c(ob$t1[j], ob$t2[j])[which.min(c(ob$freq1[j], ob$freq2[j]))])] <-
            c(ob$t1[j], ob$t2[j])[which.max(c(ob$freq1[j], ob$freq2[j])), ob$level[j]]
        }
      }
    }
  }

  # check rank discretion (only reporting if requested)
  rank_check <- check_ranks(x, ranks = rev(ranks))
  crossed_all <- rev(ob$crossed_all)
  if("discrete_rank" %in% routine) {
   out[[2]] <- rank_check
   message(paste0(length(unique(unlist(crossed_all)))), " cross-rank names identified")
  }

  # check classification
  if("duplicate_tax" %in% routine) {
    out[[3]] <- find_duplicates(x = x, ranks = ranks)
    message(paste0(nrow(out[[3]]), " conflicting classifications identified"))
    if(tax_clean) {
      rnks <- names(crossed_all)[!names(crossed_all) %in% ranks[length(ranks)]]
      for(i in 1:length(rnks)) {
        x[which(x[,rnks[i]] %in% crossed_all[[i]]), rnks[i]] <- paste0(x[which(x[,rnks[i]] %in% crossed_all[[i]]), rnks[i]], "_RNK", i)
      }
      x <- resolve_duplicates(x = x, ranks = ranks)
    }
  }

  # format output
  if(tax_clean | spell_clean) {
    if(!append) {
      # remove the duplicate classification suffixes
      x <- apply(x, 2, function(y) {gsub("[A-Z]$", "", x = y)})
      # remove the rank discrete suffixes
      x <- apply(x, 2, function(y) {gsub("_RNK*$", "", x = y)})
    }
    routine <- c(routine, "cleaned_data")
    out[[4]] <- x
  }
  # remove null elements if present
  to_remove <- is.null(out)
  if(sum(to_remove) > 0) {out <- out[!to_remove]}
  names(out) <- routine
  return(out)
}
