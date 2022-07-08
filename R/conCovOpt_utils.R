
# selectMax
selectMax <- function(x, crit = quote(con*cov), cond = quote(TRUE),  warn = TRUE){
  stopifnot(inherits(x, "conCovOpt"))
  for (i in seq_along(x)){
    xi <- x[[i]]
    attribs <- attributes(xi)[c("reprodList", "exoGroups")]
    n <- nrow(xi)
    subsetCond <- sapply(seq_len(n), 
      function(i) eval(cond, xi[i, , drop = FALSE], parent.frame()))
    if (!(is.logical(subsetCond))){
      warning("[selectMax] 'cond' is expected to evaluate to logical, but doesn't.", 
                call. = FALSE)
      subsetCond <- as.logical(subsetCond)
    }
    if (warn && n>0 && !any(subsetCond)){
      warning("[selectMax] Outcome ", names(x)[[i]], ": no solution matches the condition '", 
              deparse(substitute(cond)), "'.", call. = FALSE)
    }
    xi <- subset(xi, subsetCond) 
    n <- nrow(xi)
    if (n>0){
      xi$critValue <- sapply(seq_len(nrow(xi)), 
        function(i) eval(crit, xi[i, , drop = FALSE], parent.frame()))
      if (!(is.numeric(xi$critValue))){
        warning("[selectMax] 'crit' is expected to evaluate to numeric, but doesn't.", 
                call. = FALSE)
        xi$critValue <- as.numeric(xi$critValue)
      }
      bestVal <- max(xi$critValue)
      xi <- xi[xi$critValue == bestVal, , drop = F]
      rownames(xi) <- NULL
    } else {
      xi$critValue <- numeric(0)
    }
    xi <- xi[c("con", "cov", "critValue", "id")]
    attributes(xi)[c("reprodList", "exoGroups")] <- attribs
    x[[i]] <- xi
  }
  attr(x, "parms") <- list(crit = crit, cond = cond)
  class(x) <- "selectMax"
  x
}
print.selectMax <- function(x, ...){
  pr <- do.call(rbind, x)
  pr <- data.frame(outcome = rep(names(x), vapply(x, nrow, integer(1))), 
                   pr, row.names = NULL, stringsAsFactors = FALSE)
  sortOrder <- with(pr, order(match(outcome, unique(outcome)), -critValue))
  pr <- pr[sortOrder, , drop = FALSE]
  rownames(pr) <- NULL
  print(pr, ...)
  cat("crit:", deparse(attr(x, "parms")$crit), "\n")
  invisible(x)
}

# multipleMax
multipleMax <- function(x, outcome){
  if (!inherits(x, "selectMax")) stop("multipleMax() expects an objects of class ", dQuote("selectMax"), call. = FALSE)
  if (!all(outcome %in% names(x))) stop("Invalid 'outcome' argument.")
  warning("multipleMax() is now obsolete. It essentially returns its input. See the remark in ?selectMax.")
  attribs <- attributes(x)[c("configTable", "class", "parms")]
  x <- x[outcome]
  attributes(x)[c("configTable", "class", "parms")] <- attribs
  x
}

# ------------------------------------------------------------------------------

# Auxiliary function getIndices
# Extract indices (applicable to reprodList) from id numbers stored in selctBase
# Note: Does something similiar to cna:::rowID()
getIndices <- function(i, ll) .getInd1(i-1L, ll) + 1L
.getInd1 <- function(i, ll){
  stopifnot(i-1 <= prod(ll))
  n <- length(ll)
  if (n == 1) return(matrix(i))
  pr1 <- prod(ll[-1])
  cbind((i) %/% pr1,
        .getInd1((i) %% pr1, ll[-1]))
}
if (F){
  getIndices(1:12, 12)
  getIndices(1:12, c(2, 6))
  getIndices(1:12, 3:4)
  getIndices(1:12, 4:3)
  getIndices(1:12, c(6, 2))
  getIndices(1:24, 2:4)
  getIndices(1:24, c(3, 2, 4))
}

# ------------------------------------------------------------

# reprodAssign
reprodAssign <- function(x, outcome = names(x), id = xi$id){
  stopifnot(inherits(x, c("selectMax", "conCovOpt")), outcome %in% names(x))
  if (length(outcome) > 1){
    warning("reprodAssign() expects a single outcome - taking the first of ", length(outcome), 
            " outcomes, ", outcome[[1]], call. = FALSE)
    outcome <- outcome[[1]]
  }
  xi <- x[[outcome]]
  nid <- length(id)
  if (nrow(xi) == 0L) 
    stop("No solution for outcome ", outcome, " found.")
  poss <- attr(xi, "reprodList")
  if (!all((id%%1==0) & id>=1 & id<= prod(lengths(poss)))) 
    stop("[reprodAssign] Invalid 'id' value specified")
  ii <- getIndices(id, lengths(poss))
  exoGroups <- attr(xi, "exoGroups")
  ll <- lengths(exoGroups)
  out <- matrix(NA_real_, nrow = sum(ll), ncol = nid)
  for (i in seq_len(nid)){
    lhsSc_g <- mapply("[", poss, ii[i, ])
    out[unlist(exoGroups), i] <- rep(lhsSc_g, ll)
  }
  out
}

# ------------------------------------------------------------------------------

# DNFbuild
DNFbuild <- function(x, outcome = names(x), 
                     reduce = c("ereduce", "rreduce", "none"), id = xi$id, maxCombs = 1e7){
  stopifnot(inherits(x, c("selectMax", "conCovOpt")), outcome %in% names(x))
  if (length(outcome) > 1){
    warning("DNFbuild() expects a single outcome - taking the first of ", length(outcome), 
            " outcomes, ", outcome[[1]], call. = FALSE)
    outcome <- outcome[[1]]
  }
  # resolve reduce arg
  if (isTRUE(reduce)) reduce <- "rreduce"
  if (isFALSE(reduce) | is.null(reduce)) reduce <- "none"
  reduce <- match.arg(reduce)
  # configTable
  ct <- attr(x, "configTable")
  type <- attr(ct, "type")
  if (type == "fs") stop("DNFbuild() has no implementation of the fs case.")
  xi <- x[[outcome]]
  if (nrow(xi) == 0L) 
    stop("There is no solution for outcome ", outcome, " stored in ", deparse(substitute(x)))
  poss <- attr(xi, "reprodList")
  lhsSc <- reprodAssign(x, outcome, id)
  d <- as.data.frame(ct, warn = FALSE)
  stopifnot(nrow(d) == nrow(lhsSc))
  outcomeVar <- if (type == "mv") sub("=.+", "", outcome) else outcome
  d[[outcomeVar]] <- NULL
  dups <- duplicated(d)
  d0 <- d[!dups, , drop = FALSE]
  lhsSc <- lhsSc[!dups, , drop = FALSE]
  out <- vector("list", length(id))
  for (i in seq_along(id)){
    d <- subset(d0, lhsSc[, i]==1)
    b <- matrix(colnames(d), nrow(d), ncol(d), byrow = TRUE)
    if (type == "cs"){
      b[d == 0] <- tolower(b[d == 0])
    } else if (type == "mv"){
      b <- array(paste0(b, "=", as.matrix(d)), dim(b))
    }
    out_i <- C_recCharList2char(list(split(b, row(b))), " + ")
    if (reduce == "rreduce") out_i <- rreduce(out_i, ct, full = FALSE)
    if (reduce == "ereduce") out_i <- ereduce(out_i, ct, full = FALSE, maxCombs = maxCombs)
    out[[i]] <- out_i
  }
  unique(setdiff(unlist(out), outcome))
}

