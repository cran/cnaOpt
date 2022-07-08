
conCovOpt <- function(x, outcome = NULL, ..., rm.dup.factors = FALSE, rm.const.factors = FALSE,
											maxCombs = 1e7, approx = FALSE, allConCov = FALSE){
  if (!missing(allConCov)) 
  	warning("Argument `allConCov` is deprecated in conCovOpt() and is ignored; see the remark in ?selectMax", 
  					call. = FALSE)
  eps <- 1e-12
  ct <- configTable(x, ..., rm.dup.factors = rm.dup.factors, rm.const.factors = rm.dup.factors)
  cti <- ctInfo(ct)
  if (is.null(outcome)){
    outcome <- cti$resp_nms
  } else {
  	notNeg <- (outcome != tolower(outcome))
		outcome[notNeg] <- toupper(outcome[notNeg])
    stopifnot(outcome %in% colnames(cti$scores))
  }
  f <- attr(ct, "n")

  out <- vector("list", length(outcome))
  names(out) <- outcome
  for (outc in outcome){
    out[[outc]] <- .conCovOpt1outcome(ct, cti$scores, f, outc, eps = eps, 
    																	maxCombs = maxCombs, approx = approx)
  }
  attr(out, "configTable") <- ct
  class(out) <- "conCovOpt"
  out
}

.conCovOpt1outcome <- function(ct, sc, f, outcome, eps, maxCombs, approx = FALSE, cc1_only = FALSE){
  type <- attr(ct, "type")
  x_df <- as.data.frame(ct, warn = FALSE)
  y <- sc[, outcome] 
  outcomeVar <- if (type == "mv") sub("=.+", "", outcome) else toupper(outcome)
  
  # step 1: grouping wrt lhs-factors --------------------------------------
  ct_without_outcome <- ct[-match(outcomeVar, names(ct)), 
  												 rm.dup.factors = FALSE, rm.const.factors = FALSE]
  sc <- ctInfo(ct_without_outcome)$scores
  noGroups <- nrow(ct_without_outcome) == nrow(ct)
  if (noGroups){
    y_g <- y
    g_freqs <- rep(1L, length(y))
    yUnique <- rep(TRUE, nrow(ct_without_outcome))
    ff <- C_relist_Int(f, g_freqs)
  } else {
  	cx <- do.call(paste, c(x_df[-match(outcomeVar, names(ct))], 
  												 list(sep = "\r")))
		cx <- as.integer(factor(cx, levels = unique(cx)))
		g_freqs <- tabulate(cx, max(cx))
		u_exoGroups <- seq_along(cx)[order(cx)]
		exoGroups <- C_relist_Int(u_exoGroups, g_freqs)
    y_g <- relist1(y[u_exoGroups], g_freqs)
    yUnique <- vapply(y_g, dplyr::n_distinct, integer(1)) == 1L
  	ff <- C_relist_Int(f[u_exoGroups], g_freqs)
  }
  if (cc1_only) return(all(yUnique))

  # step 2: possible lhs-values --------------------------------------
  y1 <- rep(NA, nrow(ct_without_outcome))
  y1[yUnique] <- if (noGroups){
    y 
  } else {
    vapply(y_g[yUnique], "[", 1, FUN.VALUE = y_g[[c(1, 1)]])
  }
  yMatch <- rowAnys(sc == y1)
  yMax <- y1 >= rowMaxs(sc)
  yMin <- y1 <= rowMins(sc)

  obvious <- yUnique & rowAnys(as.matrix(data.frame(yMatch, yMax, yMin)))
  n1 <- length(y_g)
    
  possible <- vector("list", n1)
  possible[yUnique & yMatch] <- as.list(y_g[yUnique & yMatch])
  possible[yUnique & !yMatch & yMax] <- as.list(rowMaxs(sc)[yUnique & !yMatch & yMax])
  possible[yUnique & !yMatch & yMin] <- as.list(rowMins(sc)[yUnique & !yMatch & yMin])
  possible[lengths(possible)>1] <- lapply(possible[lengths(possible)>1], "[", 1L)
  sel <- vapply(possible, is.null, logical(1))
  if (any(sel) && type %in% c("cs", "mv")){
    possible[sel] <- rep_len(list(0:1), sum(sel))
  } else if (any(sel) && type == "fs"){
    yranges <- vapply(y_g[sel], range, y_g[[1]][c(1, 1)])
    sc_sorted <- t(sc[sel, , drop = FALSE])
    sc_sorted[] <- as.vector(sc_sorted)[order(as.vector(col(sc_sorted)), as.vector(sc_sorted))]
    selected <- which(sel)
    for (i in seq_len(sum(sel))){
      possible[[selected[i]]] <- getPossibleValues(unique.default(sc_sorted[, i]), yranges[, i])
    }
  }

  # Approximate mode: reduce the number of possible values considered
  if (approx){
    ff <- C_relist_Int(if (noGroups) f else f[u_exoGroups], g_freqs)
    wm <- function(x, n) median(rep(x, n))
    findClosest <- function(x, m){
      dist <- abs(x - m)
      x[dist == min(dist)]
    }
    wms <- mapply(wm, y_g, ff, SIMPLIFY = TRUE)
    possible <- mapply(findClosest, possible, wms, SIMPLIFY = FALSE)
  }
  
  #possible
  n_reprodList <- round(product(lengths(possible)))
  if (n_reprodList > maxCombs){
    warning(sprintf("[conCovOpt] Outcome %s: n_reprodList (=%4.2e) exceeds maxCombs (=%4.2e)",
                    outcome, n_reprodList, maxCombs), " - no con and cov values are returned.\n",
    				"(you may increase ", sQuote("maxCombs"), ", but calculations may take long...)", 
    				call. = FALSE)
    out <- data.frame(con = numeric(0), cov = numeric(0), id = integer(0))
    attr(out, "reprodList") <- possible
    attr(out, "exoGroups") <- if (noGroups) seq_along(g_freqs) else exoGroups
    return(out)
  }

  # step 3: preparing structures for expanding --------------------------------
  s_f <- vapply(ff, sum, integer(1))
  Sx_base <- as.vector(vapply(possible, function(x) as.numeric(x[1]), numeric(1)) %*% s_f)
                         
  Sy <- sum(y*f)
  dx <- Map(function(x, f)(x-x[1])*f, possible, s_f)
  myfn <- function(v, yr, f){
    drop(f %*% (outer(yr, v, pmin) - v[1]))
    }
  dminxy <- Map(myfn, possible, y_g, ff)
	blksize <- 10000
  out <- data.frame(setNames(C_iterate2(dx, dminxy, Sx_base, Sy, blksize = blksize), 
  													 c("con", "cov", "id")))
  attr(out, "reprodList") <- possible
  attr(out, "exoGroups") <- if (noGroups) seq_along(ff) else exoGroups
  out
}

# Aux functions
getPossibleValues <- function(possibleValues, yrange){
  csm <- sign(possibleValues - yrange[[1]]) + sign(possibleValues - yrange[[2]])
  from <- if (any(csm == -1)) possibleValues[max(which(csm == -1))] else if (any(csm == -2)) possibleValues[max(which(csm == -2))] else 0
  to <- if (any(csm == 1)) possibleValues[min(which(csm == 1))] else if (any(csm == 2)) possibleValues[min(which(csm == 2))] else 1
  possibleValues[possibleValues >= from & possibleValues <= to]
}  

getOptim <- function(x, eps = 1e-12){
  x <- data.frame(x, id = seq_along(x[[1]]))
  n <- nrow(x)
  if (n <= 1L) return(x)
  ord <- order(-x[, 1], -x[, 2])
  x <- x[ord, , drop = FALSE]
  cum <- cummax(x[, 2])
  out <- c(TRUE, x[-n, 1] > x[-1, 1] & cum[-1] > cum[-n])
  x <- x[out, , drop = FALSE]
  if (any(eq <- (diff(x$cov) < eps))) x <- x[-(which(eq)+1), , drop = FALSE] # "Fuzzy Dedup" wrt cov
  if (any(eq <- (diff(x$con) > -eps))) x <- x[-which(eq), , drop = FALSE] # "Fuzzy Dedup" wrt con
  rownames(x) <- NULL
  x
}


# print method  
print.conCovOpt <- function(x, ...){
  cat("--- conCovOpt: optimal consistency-coverage pairs ---\n")
	for (outc in names(x)){
		cat("\nOutcome ", outc, ":\n", sep = "")
		print(x[[outc]], ...)
	}
	cat("\n")
	invisible(x)
}


# plot method  
plot.conCovOpt <- function(x, con = 1, cov = 1, ...){
  d <- do.call(rbind, 
    Map(function(cc, outc) data.frame(cc[c("con", "cov")], outcome = rep(outc, nrow(cc))), 
        x, names(x)))
  ggplot(d, aes_string("con", "cov", col = "outcome")) +
    geom_rect(data = data.frame(con, cov),
              aes_string(xmin = "con", ymin = "cov", xmax = 1, ymax = 1),
              col = "lightgray", alpha = 0.05) +
    geom_point() + geom_line()
}

