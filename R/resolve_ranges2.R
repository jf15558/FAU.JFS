#' resolve_ranges2
#'
#' Function to generate a consensus age for assemblages
#' of fossil data in x, given a table of taxonomic
#' ranges. The need for error-checking is informed by
#' the error codes for the individual fossil
#' occurrences within each collection - if there is no
#' error, then the consensus age is unchanged. If
#' errors are present, then a consensus age for a
#' threshold proportion of taxa is searched for using
#' the overlap of the ranges for those taxa, as given
#' in range table y. Taxa whose occurrences lie outside
#' this consensus age are flagged as potential taxonomic
#' errors. If the threshold consensus partially overlaps
#' with the assemblage age, this overlap is returned to
#' present overzealous alteration of the age - otherwise
#' the complete consensus age is returned. If a
#' consensus age cannot be found, the original assemblage
#' age is returned, and each occurrence in the collection
#' flagged as potential taxonomic errors.
#' @param x Fossil occurrence data grouped into
#' spatiotemporally distinct assemblages
#' @param y A stratigraphic range dataset from which
#' consensus assemblage ages will be derived
#' @param assemblage The column name of the assemblage
#' groups in x
#' @param srt The column name of stratigraphic bases for
#' each element in both x and y - i.e. x and y must have
#' this same name for that column
#' @param end The column name of stratigraphic tops for
#' each element in both x and y - i.e. x and y must have
#' this same name for that column
#' @param taxon The column name denoting the taxon names
#' in both x and y - i.e. x and y must have
#' this same name for that column
#' @param err The column name flagging age errors for
#' occurrences in x. Age errors can be derived using
#' @seealso flag_ranges. Otherwise, this can be
#' fudged by adding an error column to x which will
#' result in each collection being assessed:
#' `x$err <- rep("0R0, times = nrow(x))`
#' @param prop A numeric, between 0 and 1, denoting the
#' threshold percentage of taxa in the assemblage for
#' which a consensus age must be found
#' @param allow.zero A logical determining if, in the
#' case of a collection LAD being equal to the consensus
#' age FAD (i.e. a pointwise overlap), that pointwise
#' age will be taken as the revised age. The resultant
#' collection age will have no uncertainty as a result,
#' which may be unrealistic. The default behaviour is
#' FALSE, in which case pointwise overlaps will be
#' ignored and the revised age taken instead
#' @param verbose A logical determining if the progress
#' of the redating procedure should be reported
#' @return A list of two dataframes, the first recording
#' the results of the consensus redating procedure for
#' each assemblage in x, the second recording any flags
#' (if any) for each occurrence in x
#' @import data.table
#' @importFrom stats complete.cases
#' @export

resolve_ranges2 <- function(x, y, assemblage = "collection_no", srt = "max_ma", end = "min_ma", taxon = "genus",
                            err = NULL, prop = 0.75, allow.zero = FALSE, verbose = TRUE) {

  if(!exists("x") | !exists("y")) {
    stop("Both x and y must be supplied")
  }
  if(!is.data.frame(x) | !is.data.frame(y)) {
    stop("Both x and y must be dataframes")
  }
  if(!assemblage %in% colnames(x)) {
    stop("assemblage must be a column name in x")
  }
  if(!all(c(taxon, srt, end) %in% colnames(x))) {
    stop("Arguments genus, srt and end must all be column names in x and y")
  }
  if(!all(c(taxon, srt, end) %in% colnames(y))) {
    stop("Arguments genus, srt and end must all be column names in x and y")
  }
  if(is.null(err)) {
    err <- "age_flag"
    x$age_flag <- rep("0R0", times = nrow(x))
  } else {
    if(!err %in% colnames(x)) {
      stop("err must be a colum name in x")
    }
  }
  if(length(taxon) > 1) {
    stop("taxon must a character vector of length 1")
  }
  if(!all(class(x[,srt]) == "numeric", class(x[,end]) == "numeric",
          class(y[,srt]) == "numeric", class(y[,end]) == "numeric")) {
    stop("srt and end columns in x and y must all be numeric")
  }
  if(any(x[,srt] < x[,end])) {
    stop("One or more maximum ages in x are smaller than their corresponding minimum ages")
  }
  if(any(y[,srt] < y[,end])) {
    stop("One or more maximum ages in y are smaller than their corresponding minimum ages")
  }
  if(length(prop) != 1 | class(prop) != "numeric") {
    stop("prop must be a numeric of length 1")
  }
  if(prop > 1 | prop < 0) {
    stop("prop must be greater than 0 and less than or equal to 1")
  }
  if(prop < 0.5) {
    warning("Prop is quite a low value - a minimum of 0.6 is desirable")
  }

  # global variable workaround
  . <- lb <- ub <- b <- N <- .N <- .SD <- .EACHI <- NULL
  # tabulate error codes in each collection (ignoring uncoded and correct ages R1R, 000)
  codes <- c("000", "R1R", "0R0", "00R", "R00", "0R1", "1R0")
  tabs <- data.frame(unique(x[,assemblage]), NA, NA, NA, NA, NA, NA, NA)
  for(i in 1:length(codes)) {
    vals <- tapply(x[,err], x[,assemblage], function(x) {length(which(x == codes[i]))})
    tabs[,i + 1] <- vals[match(tabs$unique.x...assemblage.., as.numeric(names(vals)))]
  }
  colnames(tabs) <- c("assemblage", codes)

  # objects to store output. This is initialised with the assemblage ages as given in the PBDB
  # and nodes for unrevised (revision = 0) and 100% congruency (prop = 1).
  # 1 revision = "revised", 2 revision = "unsolved"
  z <- data.frame(tabs$assemblage, NA, NA, "valid", 1)
  colnames(z) <- c("assemblage", "FAD", "LAD", "revision", "prop")
  z$FAD <- x[match(tabs$assemblage, x[,assemblage]), srt]
  z$LAD <- x[match(tabs$assemblage, x[,assemblage]), end]
  FAD <- x[,(srt)]
  LAD <- x[,(end)]
  tax_flag <- FAD
  tax_flag[] <- 0

  # find assemblages with errors
  to_do <- apply(tabs[,4:ncol(tabs)], 1, sum)
  to_do <- tabs$assemblage[which(to_do != 0)]
  # convert to data.table for rapid indexing
  x$rnum <- 1:nrow(x)
  x1 <- data.table::data.table(x)
  data.table::setkeyv(x1, assemblage)

  # for the incorrect assemblages, find the consensus if possible
  for(i in 1:length(to_do)) {

    # select assemblage occurrences and trim to unique values
    occs <- as.data.frame(x1[.(to_do[i]),])
    occs2 <- unique(occs[,c(taxon, srt, end)])
    occs2 <- occs2[stats::complete.cases(occs2),]
    zpos <- match(to_do[i], z$assemblage)

    # extract range chart taxa
    foo <- y[match(occs2[,taxon], y[,taxon]), c(taxon, srt, end)]
    foo <- foo[stats::complete.cases(foo),]

    # if there are range chart taxa, test for a full overlapping range
    if(nrow(foo) != 0) {
      # god-like solution from https://stackoverflow.com/questions/66754356/finding-the-overlapping-range-of-a-set-of-vectors-in-r/66758534#66758534
      dt <- data.table(lb = foo[[3]], ub = foo[[2]])
      mdt <- dt[, .(b = unique(unlist(.SD)))]
      sol <- dt[mdt, on = .(lb <= b, ub >= b), .N, by = .EACHI]
      tprop <- max(sol$N) / nrow(foo)

      # if the returned interval covers the proportion of taxa prop
      # (N = 1 means that there are no overlaps between any ranges)
      if(tprop >= prop) {

        # get the old age boundary
        ores <- as.vector(z[zpos, c("FAD", "LAD")])
        # get the assemblage span
        res <- rev(sol[N == max(N), range(lb)])
        sol2 <- rbind(ores, res)
        sol2 <- sol2[order(sol2[,1], decreasing = TRUE),]

        # if there is any overlap, including a pointwise overlap e.g. assemblage LAD = 21, overlap FAD = 21
        if(sol2[2,1] >= sol2[1,2] & sol2[2, 1] < sol2[1, 1]) {

          # if the overlap is pointwise..
          if(sol2[2,1] == sol2[1,2]) {

            # .. and if pointwise overlaps have been allowed
            if(allow.zero) {
              sol2 <- c(min(sol2[,1]), max(sol2[,2]))

            # otherwise take the revised age rather than the assemblage-overlap intersection
            } else {
              sol2 <- res
            }

          # if the overlap is not pointwise, take it by default
          } else {
            sol2 <- c(min(sol2[,1]), max(sol2[,2]))
          }

        # if there is no overlap, take the revised age
        } else {
          sol2 <- res
        }
        # update assemblage age and revision status
        z[zpos, c("FAD", "LAD")] <- sol2
        z[zpos, "revision"] <- "revised"
        z[zpos, "prop"] <- tprop
        FAD[occs$rnum] <- sol2[1]
        LAD[occs$rnum] <- sol2[2]

        # test for occurrences to flag
        if(tprop != 1) {
          foo2 <- foo[which((foo$max_ma > res[1] & foo$min_ma >= res[1]) |
                              (foo$max_ma <= res[2] & foo$min_ma < res[2])),taxon]
          tax_flag[occs$rnum[which(occs[,taxon] %in% foo2)]] <- 1
        }

      } else {
        z[zpos, "prop"] <- tprop
        z[zpos, "revision"] <- "unresolved"
        tax_flag[occs$rnum] <- 1
      }


    } else {
      z[zpos, "prop"] <- 0
      z[zpos, "revision"] <- "unresolved"
      tax_flag[occs$rnum] <- 1
    }

    # notify R
    if(verbose) {
      if(i != 1) {cat(paste0("\r"))}
      cat(paste0("Assemblage ", i, "/", length(to_do), " checked"))
    }
  }

  # return
  per_occ <- cbind.data.frame(FAD, LAD, tax_flag)
  colnames(per_occ) <- c("FAD", "LAD", "tax_flag")
  out <- list()
  out[[1]] <- z
  out[[2]] <- per_occ
  return(out)
}
