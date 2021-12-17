
# findOutcomes
findOutcomes <- function(x, con = 1, cov = 1, rm.dup.factors = FALSE, rm.const.factors = FALSE, ...){
  x <- configTable(x, rm.dup.factors = rm.dup.factors, rm.const.factors = rm.dup.factors)
  if (con == 1 && cov == 1){
  	return(con1cov1(x))
  }
  con_threshold <- con
  cov_threshold <- cov
  b <- conCovOpt(x, ...)
  sm <- selectMax(b, cond = bquote(con >= .(con_threshold) & cov >= .(cov_threshold)), 
                  warn = FALSE)
  out <- data.frame(Factor = names(b), 
                    outcome = vapply(sm, nrow, integer(1L)) > 0L)
  rownames(out) <- NULL
  out
}

con1cov1 <- function(ct){
  eps <- 1e-12
  cti <- ctInfo(ct)
  responses <- cti$resp_nms
  f <- attr(ct, "n")

  out <- logical(length(responses))
  names(out) <- responses
  for (outc in responses){
    out[[outc]] <- .conCovOpt1outcome(ct, cti$scores, f, outc, cc1_only = TRUE)
  }
  out <- data.frame(Factor = responses, outcome = out)
  rownames(out) <- NULL
  out
 }
