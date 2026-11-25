#' Randomization test for singing coordination with random assignment of missing individual IDs
#'
#' Same as \code{\link{test_coordination}} but missing values (NAs) in the 'indiv' column are randomly assigned to one of the two individuals. The assignment is repeated several times (imputations) and each imputed data set is tested against its own set of randomized data sets.
#' @inheritParams test_coordination
#' @param  X Data frame containing columns for singing event (sing.event),
#' individual (indiv), and start and end time of signal (start and end). The 'indiv' column must contain exactly 2 individual IDs (besides NAs).
#' @param prob Numeric between 0 and 1. Probability of assigning a missing ('NA') cell in the 'indiv' column to the individual in \code{id}. The remaining cells are assigned to the other individual. Values are drawn from a binomial distribution (\code{\link[stats]{rbinom}}). Only applies to signals whose ID is not determined by overlapping signals (see details). Default is 0.5.
#' @param id Character string with the individual ID to which \code{prob} applies. Must be one of the two IDs in the 'indiv' column. If \code{NULL} (default) the first ID in alphabetical order is used.
#' @param imputations Number of imputed data sets (i.e. random assignments of missing IDs) per singing event. Default is 100. Events without missing IDs are imputed only once. Note that the total number of randomizations per event is \code{imputations * iterations}.
#' @param iterations number of iterations for shuffling and calculation of the expected number of overlaps for each imputed data set. Default is 1000.
#' @param seed Integer with the random seed. If \code{NULL} (default) no seed is set. Results are only fully reproducible when \code{parallel = 1}.
#' @return A data frame as in \code{\link{test_coordination}} in which values are averaged across imputations. The data frame includes an additional column 'imputations' with the number of imputed data sets that could be analyzed. Results for each imputation are returned as the attribute 'imputation.results'.
#' @export
#' @name test_coordination_impute
#' @details
#' The test is run in two nested levels:
#' \enumerate{
#'   \item Imputation: NAs in the 'indiv' column are assigned to one of the two individuals, ensuring that overlapping signals are assigned to different individuals. Within each singing event, signals are grouped into sets of (directly or indirectly) overlapping signals. In groups containing labeled signals the IDs of the NA signals are determined by alternating IDs from the labeled ones. In groups without labeled signals the first signal is assigned to the individual in \code{id} with probability \code{prob} (and to the other individual with probability \code{1 - prob}) and IDs alternate from there. Hence, NA signals that do not overlap any other signal are assigned independently following a binomial distribution. This is repeated \code{imputations} times, creating several "observed" data sets.
#'   \item Randomization: each imputed data set is randomized \code{iterations} times as in \code{\link{test_coordination}}, keeping the ID labels of that imputed data set. The observed statistic of the imputed data set is compared against its own randomized data sets.
#' }
#' Results are then averaged across imputations. The p-value is the mean of the p-values from all imputations, which is equivalent to the proportion of all randomized data sets (across all imputations) that are as extreme as their corresponding imputed observed data set. Hence, it incorporates the uncertainty in the identity of the unlabeled signals.
#'
#' Imputations in which an individual ends up with fewer signals than \code{cutoff} (or less than 2) are discarded and drawn again (up to 100 times). When \code{rm.solo = TRUE} solo singing is removed after each imputation.
#'
#' Note that the function assumes no overlap between signals belonging to the same individual. When this cannot be achieved by the imputation (e.g. 3 signals overlapping each other) a warning is returned.
#'
#' @examples{
#' data(sim_coor_sing)
#'
#' # add some NAs to individual IDs
#' X <- sim_coor_sing
#' X$indiv[c(3, 10, 25, 150, 300)] <- NA
#'
#' # assign NAs to "indiv.1" with probability 0.7
#' test_coordination_impute(X, prob = 0.7, id = "indiv.1", imputations = 10,
#' iterations = 100, pb = FALSE, seed = 1)
#' }
#'
#' @seealso \code{\link{test_coordination}}
#' @author Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})

test_coordination_impute <-
  function(
    X = NULL,
    prob = 0.5,
    id = NULL,
    imputations = 100,
    seed = NULL,
    iterations = 1000,
    ovlp.method = "count",
    randomization = "keep.gaps",
    less.than.chance = TRUE,
    parallel = 1,
    pb = TRUE,
    rm.incomp = FALSE,
    cutoff = 2,
    rm.solo = FALSE
  ) {
    #### set arguments from options
    # get function arguments
    argms <- methods::formalArgs(test_coordination_impute)

    # get warbleR options
    opt.argms <- if (!is.null(getOption("warbleR"))) {
      getOption("warbleR")
    } else {
      SILLYNAME <- 0
    }

    # remove options not as default in call and not in function arguments
    opt.argms <- opt.argms[
      !sapply(opt.argms, is.null) & names(opt.argms) %in% argms
    ]

    # get arguments set in the call
    call.argms <- as.list(base::match.call())[-1]

    # remove arguments in options that are in call
    opt.argms <- opt.argms[!names(opt.argms) %in% names(call.argms)]

    # set options left
    if (length(opt.argms) > 0) {
      for (q in seq_len(length(opt.argms))) {
        assign(names(opt.argms)[q], opt.argms[[q]])
      }
    }

    if (!is.data.frame(X)) {
      stop2("X is not a data frame")
    }

    if (!"indiv" %in% names(X)) {
      stop2("'indiv' column not found in X")
    }

    # stop if some cells are not labeled
    if (any(is.na(X$sing.event))) {
      stop2("NA's in singing event names ('sing.event' column) not allowed")
    }

    # if there are NAs in start or end stop
    if (any(is.na(c(X$end, X$start)))) {
      stop2("NAs found in start and/or end")
    }

    if (
      !is.numeric(prob) || length(prob) != 1 || is.na(prob) || prob < 0 ||
        prob > 1
    ) {
      stop2("'prob' must be a single number between 0 and 1")
    }

    if (
      !is.numeric(imputations) || length(imputations) != 1 ||
        is.na(imputations) || imputations < 1
    ) {
      stop2("'imputations' must be a single positive number")
    }

    imputations <- round(imputations)

    # if iterations is not vector or length==1 stop
    if (any(!is.vector(iterations), !is.numeric(iterations))) {
      stop2("'interations' must be a numeric vector of length 1")
    } else {
      if (!length(iterations) == 1) {
        stop2("'interations' must be a numeric vector of length 1")
      }
    }

    # round iterations
    iterations <- round(iterations)

    # interations should be positive
    if (iterations < 2) {
      stop2("'iterations' must be > 1")
    }

    # if parallel is not numeric
    if (!is.numeric(parallel)) {
      stop2("'parallel' must be a numeric vector of length 1")
    }
    if (any(!(parallel %% 1 == 0), parallel < 1)) {
      stop2("'parallel' should be a positive integer")
    }

    #### get IDs to impute
    X$indiv <- as.character(X$indiv)
    X$sing.event <- as.character(X$sing.event)

    ids <- sort(unique(X$indiv[!is.na(X$indiv)]))

    if (length(ids) != 2) {
      stop2(paste0(
        "'indiv' column must contain exactly 2 individual IDs (",
        length(ids),
        " found)"
      ))
    }

    if (is.null(id)) {
      id <- ids[1]
    }

    if (length(id) != 1 || !id %in% ids) {
      stop2(paste0(
        "'id' must be one of the IDs in the 'indiv' column: ",
        paste(ids, collapse = ", ")
      ))
    }

    other_id <- setdiff(ids, id)

    # flag cells to be imputed
    X$imputed <- is.na(X$indiv)

    #### check events (NAs could belong to any of the 2 individuals)
    # maximum possible number of signals per individual in each event
    qw <- sapply(ids, function(i) {
      tapply(
        X$indiv %in% i | X$imputed,
        factor(X$sing.event, levels = unique(X$sing.event)),
        sum
      )
    })
    qw <- matrix(
      qw,
      ncol = 2,
      dimnames = list(unique(X$sing.event), ids)
    )

    indiv.cnt <- apply(qw > 0, 1, sum)
    sng.cnt <- apply(qw, 1, function(x) any(x < cutoff))

    # complete singing events
    if (any(indiv.cnt < 2)) {
      if (rm.incomp) {
        X <- X[X$sing.event %in% names(indiv.cnt)[indiv.cnt == 2], ]
        warning2(
          "Some events didn't have 2 interacting individuals and were excluded"
        )
      } else {
        warning2(
          "Some singing events don't have 2 interacting individuals ('indiv' column)"
        )
      }
    }

    # deal with cutoff value
    if (any(sng.cnt)) {
      X <- X[X$sing.event %in% names(indiv.cnt)[!sng.cnt], ]
      warning2(
        "Some individuals didn't have more songs that the 'cutoff' and the events were excluded"
      )
    }

    # if nothing was left
    if (nrow(X) == 0) {
      stop2("No events left after removing incomplete events")
    }

    # remove solo singing from a single event
    soloFUN <- function(Y) {
      Y <- Y[order(Y$start), ]
      Y_list <- split(Y, Y$indiv)

      fst <- max(sapply(Y_list, function(x) which(Y$start == min(x$start)))) -
        1
      lst <- min(sapply(Y_list, function(x) which(Y$start == max(x$start)))) -
        1

      if (lst > nrow(Y)) {
        lst <- nrow(Y)
      }

      Y[fst:lst, ]
    }

    # find groups of overlapping signals in a single event (Y sorted by start)
    overlapFUN <- function(Y) {
      ov <- outer(Y$start, Y$end, "<") & outer(Y$end, Y$start, ">")
      diag(ov) <- FALSE
      nghbrs <- lapply(seq_len(nrow(Y)), function(i) which(ov[i, ]))

      # connected components of the overlap graph
      comp <- rep(NA_integer_, nrow(Y))
      for (i in seq_len(nrow(Y))) {
        if (is.na(comp[i])) {
          queue <- i
          comp[i] <- i
          while (length(queue) > 0) {
            u <- queue[1]
            queue <- queue[-1]
            new <- nghbrs[[u]][is.na(comp[nghbrs[[u]]])]
            comp[new] <- i
            queue <- c(queue, new)
          }
        }
      }

      list(nghbrs = nghbrs, comps = split(seq_len(nrow(Y)), comp))
    }

    # count imputations in which overlapping signals could not be assigned to different individuals
    n.conflicts <- 0

    # assign IDs to NAs in a single event so overlapping signals belong to different individuals
    drawFUN <- function(Y, ovlps) {
      indiv <- Y$indiv
      conflict <- FALSE

      for (cmp in ovlps$comps) {
        if (!any(Y$imputed[cmp])) {
          next
        }

        # labelled signals constrain the group, if none the first signal is drawn
        seeds <- cmp[!Y$imputed[cmp]]
        if (length(seeds) == 0) {
          seeds <- cmp[1]
          indiv[seeds] <- if (stats::rbinom(1, 1, prob) == 1) id else other_id
        }

        # propagate alternating IDs through overlapping signals
        queue <- seeds
        while (length(queue) > 0) {
          u <- queue[1]
          queue <- queue[-1]
          opposite <- setdiff(ids, indiv[u])

          for (v in ovlps$nghbrs[[u]]) {
            if (is.na(indiv[v])) {
              indiv[v] <- opposite
              queue <- c(queue, v)
            } else if (indiv[v] == indiv[u] && (Y$imputed[u] || Y$imputed[v])) {
              conflict <- TRUE
            }
          }
        }
      }

      if (conflict) {
        n.conflicts <<- n.conflicts + 1
      }

      Y$indiv <- indiv
      Y
    }

    # assign IDs to NAs in a single event
    imputeFUN <- function(Y, ovlps) {
      Y0 <- Y
      for (attempt in 1:100) {
        if (any(Y0$imputed)) {
          Y <- drawFUN(Y0, ovlps)
        }

        Y2 <- if (rm.solo && length(unique(Y$indiv)) > 1) soloFUN(Y) else Y

        # each individual needs enough signals to be randomized
        n.sgnls <- table(factor(Y2$indiv, levels = ids))
        if (all(n.sgnls >= max(cutoff, 2))) {
          return(Y2)
        }

        # no point in trying again if there is nothing to impute
        if (!any(Y0$imputed)) {
          break
        }
      }

      return(NULL)
    }

    #### create imputed data sets
    if (!is.null(seed)) {
      set.seed(seed)
    }

    jobs <- list()

    for (h in unique(X$sing.event)) {
      Y <- X[X$sing.event == h, ]
      Y <- Y[order(Y$start), ]

      # overlapping signals (only needed once per event)
      ovlps <- overlapFUN(Y)

      # events without NAs need a single imputation
      n.imp <- if (any(Y$imputed)) imputations else 1

      for (m in seq_len(n.imp)) {
        jobs[[length(jobs) + 1]] <- list(
          sing.event = h,
          imputation = m,
          Y = imputeFUN(Y, ovlps)
        )
      }
    }

    if (n.conflicts > 0) {
      warning2(paste0(
        "In ",
        n.conflicts,
        " imputation(s) some overlapping signals could not be assigned to different individuals (e.g. more than 2 signals overlapping)"
      ))
    }

    # randomization function
    rndmFUN <- function(Y) {
      Y <- Y[order(Y$start), ]
      Y_list <- split(Y, Y$indiv)

      # null model
      # duration of signals
      durs <- lapply(Y_list, function(x) x$end - x$start)

      # duration of gaps
      gaps <- lapply(Y_list, function(y) {
        sapply(1:(nrow(y) - 1), function(x) {
          y$start[x + 1] - y$end[x]
        })
      })

      # randomize
      rnd.dfs <- lapply(1:iterations, function(x) {
        # randomize gaps
        if (randomization %in% c("keep.gaps", "keep.song.order")) {
          gaps <- lapply(gaps, sample)
        }

        if (randomization == "sample.gaps") {
          # generate gaps from lognormal distribution
          gaps <- lapply(gaps, function(x) {
            stats::rlnorm(
              n = length(x),
              meanlog = mean(log(unlist(gaps))),
              sdlog = stats::sd(log(unlist(gaps)))
            )
          })
        }

        # randomize signals
        if (randomization %in% c("keep.gaps", "sample.gaps")) {
          durs <- lapply(durs, sample)
        }

        # put all back together as a sequence of signals and gaps
        ndfs_list <- lapply(names(Y_list), function(x) {
          nbt <- NULL
          for (i in 1:(length(durs[[x]]) - 1)) {
            nbt[i] <- durs[[x]][i] + gaps[[x]][i]
            if (i != 1) nbt[i] <- nbt[i] + nbt[i - 1]
          }
          nbt <- c(0, nbt)
          nbt <- nbt + min(Y_list[[x]]$start)
          net <- nbt + durs[[x]]

          ndf <- data.frame(
            indiv = x,
            start = nbt,
            end = net
          )

          return(ndf)
        })

        ndfs <- do.call(rbind, ndfs_list)
        ndfs <- ndfs[order(ndfs$start), ]

        rownames(ndfs) <- 1:nrow(ndfs)

        return(ndfs)
      })

      # add observed as the first element of list
      dfs <- c(list(Y), rnd.dfs)

      return(dfs)
    }

    # counting ovlp.method
    countFUN <- function(Z) {
      # order by time and add duration
      Z <- Z[order(Z$start), ]

      Z1 <- Z[Z$indiv == unique(Z$indiv)[1], ]

      Z2 <- Z[Z$indiv == unique(Z$indiv)[2], ]

      out <- sapply(1:nrow(Z1), function(i) {
        # target start and end
        trg.strt <- Z1$start[i]
        trg.end <- Z1$end[i]

        # get the ones that overlap
        return(sum(Z2$end > trg.strt & Z2$start < trg.end))
      })

      return(sum(out))
    }

    # time.overlap
    durFUN <- function(Z) {
      # order by time and add duration
      Z <- Z[order(Z$start), ]
      Z$duration <- Z$end - Z$start

      Z1 <- Z[Z$indiv == unique(Z$indiv)[1], ]

      Z2 <- Z[Z$indiv == unique(Z$indiv)[2], ]

      out <- sapply(1:nrow(Z1), function(i) {
        # target start and end
        trg.strt <- Z1$start[i]
        trg.end <- Z1$end[i]

        # get the ones that overlap
        Z2 <- Z2[Z2$end > trg.strt & Z2$start < trg.end, , drop = FALSE]

        if (nrow(Z2) > 0) {
          # set new start and end at edges of overlaping signals
          if (any(Z2$start < trg.strt)) {
            trg.strt <- max(Z2$end[Z2$start < trg.strt])
          }

          if (any(Z2$end > trg.end)) {
            trg.end <- min(Z2$start[Z2$end > trg.end])
          }

          # new duration
          no.ovlp.dur <- trg.end - trg.strt

          ovlp <- if (no.ovlp.dur > 0) {
            Z1$duration[i] - no.ovlp.dur
          } else {
            Z1$duration[i]
          }

          return(ovlp)
        } else {
          return(0)
        }
      })

      return(sum(out))
    }

    # time.closest ovlp.method
    closestFUN <- function(Z) {
      timediffs <- vapply(
        seq_len(nrow(Z)),
        function(i) {
          Z_others <- Z[Z$indiv != Z$indiv[i], ]

          timediff <- if (
            any(Z_others$end > Z$start[i] & Z_others$start < Z$end[i])
          ) {
            0
          } else {
            min(abs(c(Z_others$start - Z$end[i], Z_others$end - Z$start[i])))
          }

          return(timediff)
        },
        FUN.VALUE = numeric(1)
      )

      return(mean(timediffs))
    }

    # select function/ovlp.method
    coortestFUN <- if (ovlp.method == "count") {
      countFUN
    } else if (ovlp.method %in% c("time.overlap", "duration")) {
      # duration kept for compatibility with previous versions
      durFUN
    } else if (ovlp.method == "time.closest") {
      # time to closest call from other individuals
      closestFUN
    }

    # set clusters for windows OS
    if (Sys.info()[1] == "Windows" & parallel > 1) {
      cl <- parallel::makePSOCKcluster(getOption("cl.cores", parallel))
    } else {
      cl <- parallel
    }

    # run loop apply function over imputed data sets
    cote <- .pblapply(
      pbar = pb,
      X = jobs,
      cl = cl,
      message = "computing coordination statistics",
      total = 1,
      FUN = function(job) {
        ovlp <- if (is.null(job$Y)) {
          NULL
        } else {
          try(sapply(rndmFUN(job$Y), coortestFUN), silent = TRUE)
        }

        if (!is.null(ovlp) && !is(ovlp, "try-error")) {
          # get observed overlap (first element)
          obs.overlaps <- ovlp[1]

          # get random overlap (all except the first element)
          rov <- ovlp[-1]
          mean.random.ovlps <- mean(rov)

          # calculate p-value
          if (less.than.chance) {
            p <- length(rov[rov <= obs.overlaps]) / iterations
          } else {
            p <- length(rov[rov >= obs.overlaps]) / iterations
          }

          # coordination score
          if (obs.overlaps == 0 & mean.random.ovlps == 0) {
            coor.score <- 0
          } else {
            coor.score <- round(
              (obs.overlaps - mean.random.ovlps) / mean.random.ovlps,
              3
            )
          }

          # bounded coordination score
          if (obs.overlaps == 0 & mean.random.ovlps == 0) {
            coor.score.bounded <- 0
          } else {
            coor.score.bounded <- round(
              (obs.overlaps - mean.random.ovlps) /
                (obs.overlaps + mean.random.ovlps),
              3
            )
          }

          l <- data.frame(
            sing.event = job$sing.event,
            imputation = job$imputation,
            obs.ovlp = obs.overlaps,
            mean.random.ovlp = mean.random.ovlps,
            p.value = p,
            coor.score,
            coor.score.bounded
          )
        } else {
          l <- data.frame(
            sing.event = job$sing.event,
            imputation = job$imputation,
            obs.ovlp = NA,
            mean.random.ovlp = NA,
            p.value = NA,
            coor.score = NA,
            coor.score.bounded = NA
          )
        }

        return(l)
      }
    )

    imp.df <- do.call(rbind, cote)
    rownames(imp.df) <- NULL

    # average across imputations
    df_list <- lapply(unique(imp.df$sing.event), function(h) {
      W <- imp.df[imp.df$sing.event == h, ]
      W <- W[!is.na(W$obs.ovlp), ]

      if (nrow(W) == 0) {
        return(data.frame(
          sing.event = h,
          obs.ovlp = NA,
          mean.random.ovlp = NA,
          p.value = NA,
          coor.score = NA,
          coor.score.bounded = NA,
          imputations = 0
        ))
      }

      data.frame(
        sing.event = h,
        obs.ovlp = mean(W$obs.ovlp),
        mean.random.ovlp = mean(W$mean.random.ovlp),
        p.value = mean(W$p.value),
        coor.score = round(mean(W$coor.score), 3),
        coor.score.bounded = round(mean(W$coor.score.bounded), 3),
        imputations = nrow(W)
      )
    })

    df <- do.call(rbind, df_list)

    # add results per imputation
    attr(df, "imputation.results") <- imp.df

    return(df)
  }
