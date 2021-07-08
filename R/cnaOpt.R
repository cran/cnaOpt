
#   x           A data.frame or configTable (cs only!)
#   outcome     Name of the outcome (A column in x)
#   ...         Passed to configTable
#   crit, cond  Passed to selectMax
cnaOpt <- function(x, outcome, ..., reduce = c("ereduce", "rreduce", "none"), 
									 niter = 1, crit = quote(con * cov), cond = quote(TRUE), 
									 approx = FALSE, maxCombs = 1e7){
  time.init <- Sys.time()
  ct <- configTable(x, ...)
	if (attr(ct, "type") == "fs") 
		stop("cnaOpt() has no implementation of the fs case.")  
  outcomeVar <- if (attr(ct, "type") == "mv") sub("=.+", "", outcome) else (outcome)
  if (!(outcomeVar %in% names(ct)) || length(outcomeVar) != 1)
  	stop("Invalid specification of ", dQuote("outcome"))
  # Calculate con-cov optima for outcome
  cco <- conCovOpt(ct, outcome, approx = approx, maxCombs = maxCombs) 
  if (nrow(cco[[outcome]]) == 0 ||
  	  !any(is.finite(cco[[outcome]]$con) & is.finite(cco[[outcome]]$cov))) 
  	stop("conCovOpt() did not return valid consistency and coverage values.")
  # Calculate the con-cov maximum for outcome and scores
  best <- selectMax(cco, crit = crit, cond = cond)
  scores_best <- as.vector(reprodAssign(best, outcome = outcome))
  
  # get cond and cond_neg
	cond <- mat2charList(ct, scores_best == 1L, outcomeVar)
	cond_neg <- mat2charList(ct, scores_best == 0L, outcomeVar)

	# reduction procedure, depending on arg reduce  
  if (isTRUE(reduce)) 
      reduce <- "rreduce"
  if (isFALSE(reduce) | is.null(reduce)) 
      reduce <- "none"	
	reduce <- match.arg(reduce)
	if (!missing(niter) && reduce != "rreduce")
		warning(sQuote("niter"), " is ignored, as reduce!=\"rreduce\"")
	if (reduce == "ereduce"){
		mhs <- MBproc(cond, cond_neg, ctInfo(ct)$sc) # MB's minization procedure
		outStr <- C_mconcat(mhs, sep = "+")           # Formulate into strings
	} else {
		outStr <- paste0(C_mconcat(cond, "*"), collapse = "+")
		if (reduce == "rreduce"){
			if (niter == 1){
				outStr <- rreduce(outStr, 
													x = ct, full = FALSE,
													simplify2constant = FALSE)
			} else if (niter>1){
				# repeated attempts with rreduce()
				rred_attempts <- replicate(niter, 
						rreduce(outStr, x = ct, full = FALSE,
										simplify2constant = FALSE),
					simplify = TRUE)
				dups <- duplicated(stdCond(rred_attempts))
				outStr <- rred_attempts[!dups]
			}
		}
	}
  if (length(outStr)) outStr <- paste0(outStr, "<->", outcome)
  class(outStr) <- c("stdAtomic", "character")
  nsol <- length(outStr)

  # Output
  out <- data.frame(outcome = structure(rep(outcome, nsol), class = c("outcomeString", "character")),
                    condition = outStr,
                    consistency = best[[outcome]]$con,
                    coverage = best[[outcome]]$cov,
                    complexity = lengths(gregexpr("[\\+\\*]", outStr)) + 1L,
                    stringsAsFactors = FALSE)
  out <- out[order(out$complexity), , drop = FALSE]
  out <- structure(out,
                   ct = ct,
                   scores = scores_best,
                   timing = Sys.time() - time.init,
                   class = c("cnaOpt", "condTbl", "data.frame"))
  rownames(out) <- NULL
  out
}

# ==============================================================================

# AUXILIARY FUNCTIONS
# -------------------
# (used in the functions above)

# Extract cond in "charList" format from a matrix or data.frame
#   x must be a configTable!
mat2charList <- function(x, which, rmCol = NULL){
	if (!any(which)) return(list())
	x <- x[which, setdiff(colnames(x), rmCol)]
	unname(getCond(x, asf = NULL))
}

# Multiple application of a function 
# (non-vectorized/transposed version of base::outer)
# OUTER <- function(x, y, FUN, ...){
#   array(mapply(FUN, rep(x, each = length(y)), rep(y, length(x)), ...,
#                SIMPLIFY = FALSE),
#         dim = c(length(y), length(x)))
# }
OUTER <- function(x, y, FUN, ...){
	lx <- length(x)
	ly <- length(y)
	out <- vector("list", lx*ly)
	dim(out) <- c(ly, lx)
	FUN <- match.fun(FUN)
	for (i in seq_along(out)-1L){
		ix <- i %/% ly + 1
		iy <- i %% ly + 1
		out[[i+1]] <- FUN(x[[ix]], y[[iy]], ...)
	}
	out
}

# simple auxiliary function
contains <- function(x, y) all(y %in% x)

# Find all minimal hitting sets for a collection of sets
minimalHittingSets <- function(x){
  if (length(x) == 0) return(character(0))
  elements <- sort(unique(unlist(x)))
  l <- length(elements)
  sol <- list()
  i <- 0L
  repeat {
    i <- i+1L
    cand <- combn(elements, i, simplify = FALSE)
    toremove <- C_m_contains_one_of(cand, sol)
    if (all(toremove)) break
    cand <- cand[!toremove]
    ok <- C_m_intersects_with_all(cand, x)
    sol <- c(sol, cand[ok])
    if (i >= l) break
  }
  sol
} 
