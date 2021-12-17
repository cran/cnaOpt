
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
	if (outcome != tolower(outcome)) outcome <- toupper(outcome)
  outcomeVar <- if (attr(ct, "type") == "mv") sub("=.+", "", outcome) else toupper(outcome)
  if (!(outcomeVar %in% names(ct)) || length(outcomeVar) != 1)
  	stop("Invalid specification of ", dQuote("outcome"))
  # Calculate con-cov optima for outcome
  cco <- conCovOpt(ct, outcome, approx = approx, maxCombs = maxCombs) 
  if (nrow(cco[[outcome]]) == 0 ||
  	  !any(is.finite(cco[[outcome]]$con) & is.finite(cco[[outcome]]$cov))) 
  	stop("conCovOpt() did not return valid consistency and coverage values.")
  # Calculate the con-cov maximum for outcome and scores
  best <- selectMax(cco, crit = crit, cond = cond, warn = FALSE)
  if (nrow(best[[1]]) > 0){
	  scores_best <- reprodAssign(best, outcome = outcome)
	  outCond <- list()
	
	  for (i in seq_len(ncol(scores_best))){
	 	  # get cond and cond_neg
	  	sc_best_i <- scores_best[, i]
			cond <- mat2charList(ct, sc_best_i == 1L, outcomeVar)
			cond_neg <- mat2charList(ct, sc_best_i == 0L, outcomeVar)
		
			# reduction procedure, depending on arg reduce  
		  if (isTRUE(reduce)) 
		      reduce <- "rreduce"
		  if (isFALSE(reduce) | is.null(reduce)) 
		      reduce <- "none"	
			reduce <- match.arg(reduce)
			if (!missing(niter) && reduce != "rreduce")
				warning(sQuote("niter"), " is ignored, as reduce!=\"rreduce\"")
			if (reduce == "ereduce"){
				mhs <- MBproc(cond, cond_neg, ctInfo(ct)$sc, # MB's minization procedure
											maxCombs = maxCombs) 
				outStr <- C_mconcat(mhs, sep = "+")          # Formulate into strings
			} else {
				outStr <- paste0(C_mconcat(cond, "*"), collapse = "+")
				if (reduce == "rreduce"){
					outStr <- rreduce(outStr, x = ct, niter = niter, full = FALSE,
														simplify2constant = FALSE)
				}
			}
			outCond[[i]] <- outStr
	  }
	  ll <- lengths(outCond)
	  outCond <- do.call(c, outCond)
  }
  else {
  	scores_best <- matrix(numeric(0), nrow = nrow(ct), ncol = 0)
  	outCond <- character(0)
  	ll <- integer(0)
  }
	
  if (length(outCond)) outCond <- paste0(outCond, "<->", outcome)
  class(outCond) <- c("stdAtomic", "character")
  nsol <- length(outCond)

  # Output
  out <- data.frame(outcome = structure(rep(outcome, nsol), class = c("outcomeString", "character")),
                    condition = outCond,
                    consistency = rep(best[[outcome]]$con, ll), 
                    coverage = rep(best[[outcome]]$cov, ll),
                    complexity = lengths(gregexpr("[\\+\\*]", outCond)) + 1L,
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
