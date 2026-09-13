## ---------------------------------------------------------------------------
## maxlfq.R -- MaxLFQ protein/peptide quantification in pure base R.
##
## Reimplementation of the MaxLFQ algorithm (Cox et al. 2014,
## doi:10.1074/mcp.M113.031591) with no package dependencies at all.
##
## Public functions:
##   maxlfq()        -- long-format data (report) -> group x sample matrix
##   maxlfq_matrix() -- single ids x samples matrix -> one profile
##
## Usage:
##   source("maxlfq.R")
##   pg <- maxlfq(report)                                  # precursor -> protein
##   pep <- maxlfq(report, group.header = "Stripped.Sequence")  # -> peptide
##   report <- report[report$Q.Value <= 0.01 &
##                    report$PG.Q.Value <= 0.01, ]         # filter first
##
## Scale convention: MaxLFQ determines a profile only up to a constant factor
## per group. Here each connected component is anchored so that the mean of the
## fitted log-intensities equals the mean of the observed log-intensities of
## that component (the standard convention; also what the `iq` package uses).
## ---------------------------------------------------------------------------


## --- helpers ---------------------------------------------------------------

## Column-wise medians ignoring NA. One radix sort for the whole matrix instead
## of one median() call per column, which is what makes the pairwise-ratio step
## cheap.
col_medians <- function(m) {
  nr <- nrow(m); nc <- ncol(m)
  if (nr == 0L || nc == 0L) return(rep(NA_real_, nc))
  v <- as.numeric(m)
  na <- is.na(v)
  k <- as.integer(.colSums(!na, nr, nc))
  v[na] <- Inf                                         # NAs sort to the end
  col <- rep.int(seq_len(nc), rep.int(nr, nc))
  v <- v[order(col, v, method = "radix")]              # sorted within column
  start <- (seq_len(nc) - 1L) * nr
  lo <- start + (k + 1L) %/% 2L
  hi <- start + k %/% 2L + 1L
  empty <- k == 0L
  lo[empty] <- 1L; hi[empty] <- 1L
  res <- (v[lo] + v[hi]) / 2
  res[empty] <- NA_real_
  res
}

## Medians of m[, ii] - m[, jj] for all pairs, computed in blocks so that the
## temporary difference matrix never gets large.
pair_medians <- function(m, ii, jj, block = 2e6) {
  np <- length(ii)
  if (np == 0L) return(numeric(0))
  out <- numeric(np)
  nr <- nrow(m)
  step <- max(1L, as.integer(block %/% max(1L, nr)))
  for (s in seq.int(1L, np, by = step)) {
    e <- min(s + step - 1L, np)
    k <- s:e
    out[k] <- col_medians(m[, ii[k], drop = FALSE] - m[, jj[k], drop = FALSE])
  }
  out
}

## Connected components of a symmetric logical adjacency matrix, by vectorised
## BFS (2 iterations when the graph is fully connected, the usual case).
components <- function(adj) {
  n <- nrow(adj)
  comp <- integer(n)
  cur <- 0L
  repeat {
    todo <- which(comp == 0L)
    if (!length(todo)) break
    cur <- cur + 1L
    seen <- logical(n)
    front <- todo[1L]
    seen[front] <- TRUE
    repeat {
      nb <- which(.colSums(adj[front, , drop = FALSE], length(front), n) > 0 & !seen)
      if (!length(nb)) break
      seen[nb] <- TRUE
      front <- nb
    }
    comp[seen] <- cur
  }
  comp
}


## --- core solver -----------------------------------------------------------

## m: ids x samples matrix of *log* intensities, NA = missing.
## Returns a vector of length ncol(m) of log intensities.
maxlfq_solve <- function(m, margin = -10, min.pairs = 1L) {
  nr <- nrow(m); nc <- ncol(m)
  if (nr == 1L) return(as.numeric(m[1L, ]))
  present <- matrix(as.numeric(!is.na(m)), nr, nc)
  empty <- .colSums(present, nr, nc) == 0                # samples with no data
  if (any(empty)) {
    out <- rep(NA_real_, nc)
    if (!all(empty)) out[!empty] <- maxlfq_solve(m[, !empty, drop = FALSE],
                                                 margin, min.pairs)
    return(out)
  }
  if (nc == 1L) return(max(m, na.rm = TRUE))

  ## how many ids each pair of samples has in common -- one BLAS call
  shared <- crossprod(present)
  adj <- shared >= min.pairs
  diag(adj) <- FALSE

  e <- which(adj & upper.tri(adj), arr.ind = TRUE, useNames = FALSE)
  if (nrow(e) == 0L) return(colMeans(m, na.rm = TRUE))   # nothing connected
  ii <- e[, 1L]; jj <- e[, 2L]

  ## median log-ratio for every connected pair
  r <- pair_medians(m, ii, jj)
  ok <- !is.na(r)
  if (!all(ok)) {
    ii <- ii[ok]; jj <- jj[ok]; r <- r[ok]
    adj[] <- FALSE
    adj[cbind(ii, jj)] <- TRUE
    adj[cbind(jj, ii)] <- TRUE
  }

  ## least squares: minimise sum over connected pairs of (x_i - x_j - r_ij)^2
  ## normal equations are a graph Laplacian, solved per component with the
  ## scale fixed by a mean constraint.
  R <- matrix(NA_real_, nc, nc)
  R[cbind(ii, jj)] <- r
  R[cbind(jj, ii)] <- -r
  lap <- -(adj * 1)
  diag(lap) <- rowSums(adj)
  rhs <- rowSums(R, na.rm = TRUE)

  comp <- components(adj)
  out <- rep(NA_real_, nc)
  for (g in seq_len(max(comp))) {
    cc <- which(comp == g)
    obs <- m[, cc, drop = FALSE]
    anchor <- mean(obs[!is.na(obs)])
    k <- length(cc)
    if (k == 1L) { out[cc] <- anchor; next }            # isolated sample
    A <- matrix(0, k + 1L, k + 1L)                      # KKT system
    A[seq_len(k), seq_len(k)] <- lap[cc, cc]
    A[k + 1L, seq_len(k)] <- 1
    A[seq_len(k), k + 1L] <- 1
    sol <- tryCatch(solve(A, c(rhs[cc], k * anchor))[seq_len(k)],
                    error = function(e) col_medians(obs))
    out[cc] <- sol
  }
  out[out <= margin] <- NA_real_
  out
}


## --- user-facing -----------------------------------------------------------

#' MaxLFQ profile for a single group.
#'
#' @param m ids x samples matrix of intensities (linear scale, NA or 0 = missing)
#' @param margin intensities below exp(margin) are treated as missing
#' @param min.pairs minimum number of shared ids required to use a sample pair
#' @return named numeric vector, linear scale
maxlfq_matrix <- function(m, margin = -10, min.pairs = 1L) {
  m <- as.matrix(m)
  storage.mode(m) <- "double"
  if (any(m < 0, na.rm = TRUE)) stop("Only non-negative quantities accepted")
  m[m == 0] <- NA_real_
  m <- log(m)
  m[m <= margin] <- NA_real_
  keep <- rowSums(!is.na(m)) > 0
  if (!any(keep)) return(setNames(rep(NA_real_, ncol(m)), colnames(m)))
  setNames(exp(maxlfq_solve(m[keep, , drop = FALSE], margin, min.pairs)),
           colnames(m))
}

#' MaxLFQ quantification from a long-format report.
#' @param x data frame with sample / group / id / quantity columns
#' @param margin quantities below exp(margin) are treated as missing
#' @param min.pairs minimum number of shared ids required to use a sample pair
#' @return groups x samples matrix of intensities
#' @export
maxlfq <- function(x,
                   sample.header = "File.Name",
                   group.header = "Protein.Names",
                   id.header = "Precursor.Id",
                   quantity.header = "Precursor.Normalised",
                   margin = -10, min.pairs = 1L) {

  need <- c(sample.header, group.header, id.header, quantity.header)
  miss <- need[!need %in% names(x)]
  if (length(miss)) stop("Missing column(s): ", paste(miss, collapse = ", "))

  smp <- as.character(x[[sample.header]])
  grp <- as.character(x[[group.header]])
  idv <- as.character(x[[id.header]])
  q   <- as.numeric(x[[quantity.header]])

  if (any(q < 0, na.rm = TRUE)) stop("Only non-negative quantities accepted")
  if (margin > -2) { margin <- -2;  warning("margin reset to -2.0") }
  if (margin < -20) { margin <- -20; warning("margin reset to -20.0") }

  keep <- !is.na(q) & q > 0 & !is.na(grp) & grp != "" &
          !is.na(smp) & !is.na(idv) & idv != ""
  smp <- smp[keep]; grp <- grp[keep]; idv <- idv[keep]
  q <- log(q[keep])
  keep <- q > margin
  smp <- smp[keep]; grp <- grp[keep]; idv <- idv[keep]; q <- q[keep]
  if (!length(q)) stop("No quantities left after filtering")

  ## integer codes in order of first appearance (keeps the original row/column
  ## order of diann_maxlfq) -- all downstream work is on integers.
  groups  <- unique(grp); g <- match(grp, groups)
  samples <- unique(smp); s <- match(smp, samples)
  ids     <- unique(idv); p <- match(idv, ids)

  ## sort by group / id / sample, best quantity first, then drop duplicates:
  ## this is the "maximum of multiple quantities per id" rule, vectorised.
  o <- order(g, p, s, -q, method = "radix")
  g <- g[o]; p <- p[o]; s <- s[o]; q <- q[o]
  n <- length(q)
  if (n > 1L) {
    dup <- c(FALSE, g[-1L] == g[-n] & p[-1L] == p[-n] & s[-1L] == s[-n])
    if (any(dup)) {
      warning("Multiple quantities per id: the maximum of these will be used")
      g <- g[!dup]; p <- p[!dup]; s <- s[!dup]; q <- q[!dup]
    }
  }

  ## group blocks are now contiguous: no scanning of the whole table per group
  starts <- which(!duplicated(g))
  ends <- c(starts[-1L] - 1L, length(g))

  res <- matrix(NA_real_, length(groups), length(samples),
                dimnames = list(groups, samples))
  for (b in seq_along(starts)) {
    k <- starts[b]:ends[b]
    cols <- s[k]; rows <- p[k]
    ucols <- unique(cols); urows <- unique(rows)
    m <- matrix(NA_real_, length(urows), length(ucols))
    m[cbind(match(rows, urows), match(cols, ucols))] <- q[k]
    res[g[starts[b]], ucols] <- maxlfq_solve(m, margin, min.pairs)
  }
  exp(res)
}