
ereduce <- function(cond, x = full.ct(cond), full = !missing(x),
                    simplify2constant = TRUE, maxCombs = 1e7){
  if (!inherits(x, "cti")) {
    if (!inherits(x, "configTable")) {
      x <- configTable(x, rm.dup.factors = FALSE, rm.const.factors = FALSE)
    }
  }
  if (attr(x, "type") == "fs") stop("Invalid use of data of type 'fs'." )
  if (full) x <- full.ct(x)
  cti <- ctInfo(x)
  stopifnot(length(cond) == 1L, nrow(cti) > 1L)
  cond <- noblanks(cond)
  sc <- cti$scores
  evalCond0 <- drop(qcond_bool(cond, sc))
  if (simplify2constant && all(evalCond0 == evalCond0[1L])) return(as.character(evalCond0[1L]))
  # convert cond to charList format 
  cond <- hstrsplit(cond, c("+", "*"), split.attr = FALSE)[[1]]
  # negative condition (as charList)
  cond_neg <- mat2charList(x, evalCond0 == 0)
  mhs <- MBproc(cond, cond_neg, sc, maxCombs = maxCombs)
  # Formulate result as char vector
  C_mconcat(mhs, sep = "+")
}

# Aux fns

# MB's minization procedure
MBproc <- function(cond_pos, cond_neg, sc, maxCombs = 1e7, verbose = TRUE){
  if (length(cond_neg)){
    d <- OUTER(cond_pos, cond_neg, setdiff)
    dd <- apply(d, 2, minimalHittingSets, maxCombs = maxCombs)
    mhs <- minimalHittingSets(lapply(dd, C_mconcat, sep = "*"), maxCombs = maxCombs, 
    													verbose = verbose)
  } else {
    mhs <- list(unique(unlist(cond_pos)))
  }
  # Check for redundant msc's and eliminate them, if any
  removeRedundantMsc(mhs, sc)
}  


# Find all minimal hitting sets for a collection of sets
minimalHittingSets <- function(x, maxCombs = 1e7, verbose = FALSE){
  if (length(x) == 0) return(character(0))
  elements <- sort(unique(unlist(x)))
  n <- length(elements)
  if (2^n > maxCombs){
  	stop(sprintf("[ereduce] Number of combinations to process (=%4.2e) exceeds maxCombs (=%4.2e)",
  	             2^n, maxCombs), 
  	     " - execution is interrupted.\n",
  	     "(you may increase ", sQuote("maxCombs"), ", but calculations may take long...)", 
  	     call. = FALSE)
  }
  el <- seq_len(n)
	xx <- lapply(x, match, elements)
  membership_x <- sapply(xx, function(x) el %in% x)
  dim(membership_x) <- c(n, length(xx))
  sol <- list()
  m <- 0L
	repeat {
		m <- m+1L
		membership_sol <- vapply(sol, function(x) ((1:n)-1) %in% x, FUN.VALUE = logical(n)) 
		dim(membership_sol) <- c(n, length(sol))
		newSol <- C_mhs_iteration(m, membership_x, membership_sol)
		sol <- c(sol, newSol)
		if (verbose) 
			cat("\rereduce() progress: ~", ceiling(sum(choose(n, 1:m)) / 2^n*100), "%", sep = "")
		if (m == n) break
	}
  if (verbose) cat("\r                         \r")
  lapply(sol, function(s) elements[s+1])
}


removeRedundantMsc <- function(mhs, sc){
  allMsc <- unique.default(unlist(mhs))
  cols <- lapply(mhs, match, allMsc)
  qc <- qcond_bool(allMsc, sc)
  sol <- list()
  repeat{
    red <- lapply(cols, function(i) C_redund(!qc[, i, drop = F]))
    anyRed <- lengths(red)>1 & vapply(red, any, logical(1))
    sol <- c(sol, cols[!anyRed])
    if (!any(anyRed)) break
    reductfn <- function(i, r) lapply(which(r), function(w) i[-w])
    lcols <- mapply(reductfn, cols, red, SIMPLIFY = FALSE)
    cols <- unique.default(unlist(lcols, recursive = FALSE))
  }
  C_mconcat(lapply(cols, function(v) allMsc[v]), 
            sep = "+")
}
