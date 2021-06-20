#' resolve_ranges
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

resolve_ranges <- function(x, y, assemblage = "collection_no", srt = "max_ma", end = "min_ma", taxon = "genus",
                           prop = 0.75, allow.zero = TRUE, verbose = TRUE) {

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
  #if(is.null(err)) {
  #  err <- "age_flag"
  #  x$age_flag <- rep("0R0", times = nrow(x))
  #} else {
  #  if(!err %in% colnames(x)) {
  #    stop("err must be a colum name in x")
  #  }
  #}
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
  #codes <- c("000", "R1R", "0R0", "00R", "R00", "0R1", "1R0")
  #tabs <- data.frame(unique(x[,assemblage]), NA, NA, NA, NA, NA, NA, NA)
  #for(i in 1:length(codes)) {
  #  vals <- tapply(x[,err], x[,assemblage], function(x) {length(which(x == codes[i]))})
  #  tabs[,i + 1] <- vals[match(tabs$unique.x...assemblage.., as.numeric(names(vals)))]
  #}
  #colnames(tabs) <- c("assemblage", codes)

  # objects to store output. This is initialised with the assemblage ages as given in the PBDB
  # and nodes for unrevised (revision = 0) and 100% congruency (prop = 1).
  # 1 revision = "revised", 2 revision = "unsolved"
  #z <- data.frame(tabs$assemblage, NA, NA, "valid", 1)
  z <- data.frame(unique(x[,assemblage]), NA, NA, NA, "000", NA)
  colnames(z) <- c("assemblage", "FAD", "LAD", "revision", "status", "prop")
  z$FAD <- x[match(z$assemblage, x[,assemblage]), srt]
  z$LAD <- x[match(z$assemblage, x[,assemblage]), end]
  FAD1 <- FAD <- x[,(srt)]
  LAD1 <- LAD <- x[,(end)]
  tax_flag <- FAD_diff <- LAD_diff <- FAD
  tax_flag[] <- "000"
  FAD_diff[] <- LAD_diff[] <- NA

  # find assemblages with errors
  #to_do <- apply(tabs[,4:ncol(tabs)], 1, sum)
  #to_do <- tabs$assemblage[which(to_do != 0)]
  # convert to data.table for rapid indexing
  x$rnum <- 1:nrow(x)
  x1 <- data.table::data.table(x)
  data.table::setkeyv(x1, assemblage)

  # for the incorrect assemblages, find the consensus if possible
  #for(i in 1:length(to_do)) {
  for(i in 1:nrow(z)) {

    # select assemblage occurrences and trim to unique values
    #occs <- as.data.frame(x1[.(to_do[i]),])
    occs <- as.data.frame(x1[.(z$assemblage[i]),])
    occs2 <- unique(occs[,c(taxon, srt, end)])
    occs2 <- occs2[stats::complete.cases(occs2),]

    # extract range chart taxa
    foo <- y[match(occs2[,taxon], y[,taxon]), c(taxon, srt, end)]
    foo <- foo[stats::complete.cases(foo),]
    tprop <- NA
    revise <- TRUE
    resolution <- "Revised"

    # if there are range chart taxa, test for an overlapping range
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
        ores <- as.vector(z[i, c("FAD", "LAD")])
        # get the assemblage span
        res <- rev(sol[N == max(N), range(lb)])
        sol2 <- rbind(ores, res)
        sol2 <- sol2[order(sol2[,1], decreasing = TRUE),]

        # if pointwise overlap in ages has been allowed, go ahead and assess overlap, taking the revised age if no overlap with assemblage age
        if(allow.zero) {

          if(sol2[2,1] >= sol2[1,2]) {
            sol2 <- c(min(sol2[,1]), max(sol2[,2]))
          } else {sol2 <- res}

        # if pointwise overlaps have not been allowed,
        } else {

          # if there are no pointwise overlaps
          if(res[1] != res[1] & sol2[2,1] != sol2[1,2]) {

            # if the assemblage and consensus ages overlap, take the overlap, otherwise take the revised age
            if(sol2[2,1] >= sol2[1,2]) {
              sol2 <- c(min(sol2[,1]), max(sol2[,2]))
            # otherwise take the consensus age
            } else {sol2 <- res}

          # otherwise pointwise overlaps have been disallowed so ignore
          } else {
            revise <- FALSE
            revision <- "Unresolved"
          }
        }

      # otherwise the proportion is not high enough so ignore
      } else {
        revise <- FALSE
        revision <- "Unresolved"
      }

    # otherwise there are no taxa to inform revision, so ignore
    } else {
      revise <- FALSE
      revision <- "Unrevised"
    }

    # if a solution has been found, revise
    if(revise) {

      # update assemblage age and revision status
      z[i, c("FAD", "LAD")] <- sol2
      z[i, "revision"] <- revision
      z[i, "prop"] <- tprop
      FAD[occs$rnum] <- sol2[1]
      LAD[occs$rnum] <- sol2[2]

      # flag occurrences if the consensus is not 100%
      if(tprop != 1) {

        # occurrences exceeding FAD or LAD
        foo2 <- foo[which((foo$max_ma > res[1] & foo$min_ma >= res[1])),taxon]
        foo3 <- foo[which((foo$max_ma <= res[2] & foo$min_ma < res[2])),taxon]
        rnum1 <- occs$rnum[which(occs[,taxon] %in% foo2)]
        rnum2 <- occs$rnum[which(occs[,taxon] %in% foo3)]

        # taxon-wise flagging
        tax_flag[rnum1] <- "00R"
        tax_flag[rnum2] <- "R00"
        FAD_diff[rnum1] <- FAD1[rnum1] - FAD[rnum1]
        LAD_diff[rnum2] <- LAD[rnum2] - LAD1[rnum2]

        # collection-wise flagging
        if(nrow(foo2) != 0 & nrow(foo3) != 0) {
          coll_flag <- "0R0"
          z[i, "status"] <- "0R0"
        } else {
          if(nrow(foo2) != 0) {
            coll_flag <- "00R"
            z[i, "status"] <- "00R"
          }
          if(nrow(foo3) != 0) {
            coll_flag <- "R00"
            z[i, "status"] <- "R00"
          }
        }

      # otherwise the revision is completely consistent
      } else {
        tax_flag[occs$rnum] <- "R1R"
        tax_flag[occs$rnum] <- "R1R"
        z[i, "status"] <- "R1R"
      }

    # otherwise assign unrevised
    } else {
      z[i, "prop"] <- tprop
      z[i, "revision"] <- revision
    }

    # notify R
    if(verbose) {
      if(i != 1) {cat(paste0("\r"))}
      cat(paste0("Assemblage ", i, "/", nrow(z), " checked"))
    }
  }

  # return
  per_occ <- cbind.data.frame(FAD, LAD, FAD_diff, LAD_diff, tax_flag)
  colnames(per_occ) <- c("FAD", "LAD", "FAD_diff", "LAD_diff", "tax_flag")
  out <- list()
  out[[1]] <- z
  out[[2]] <- per_occ
  return(out)
}
