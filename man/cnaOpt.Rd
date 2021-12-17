\name{cnaOpt}
\alias{cnaOpt}

\title{Find atomic solution formulas with optimal consistency and coverage
}

\description{
\code{cnaOpt} attempts to find atomic solution formulas (asfs) for a given \code{outcome} (inferred from crisp-set, \code{"cs"}, or multi-value, \code{"mv"}, data) that are optimal with respect to the model fit parameters consistency and coverage (cf. Baumgartner and Ambuehl 2021). 
}

\usage{
cnaOpt(x, outcome, ..., reduce = c("ereduce", "rreduce", "none"), 
       niter = 1, crit = quote(con * cov), cond = quote(TRUE), 
			 approx = FALSE, maxCombs = 1e7)
}

\arguments{
  \item{x}{
A \code{\link{data.frame}} or \code{\link{configTable}} of type \code{"cs"} or \code{"mv"}.
}
  \item{outcome}{
A character string specifying one outcome, i.e. one factor value in \code{x}.
}
  \item{\dots}{
Additional arguments passed to \code{\link{configTable}}, for instance \code{rm.dup.factors}, \code{rm.dup.factors}, or \code{case.cutoff}.
}
  \item{reduce}{
A character string: if \code{"ereduce"} or \code{"rreduce"}, the canonical DNF realizing the con-cov optimum is freed of redundancies using \code{\link{ereduce}} or \code{\link{rreduce}} (possibly repeatedly, see \code{niter}), respectively; if \code{"none"}, the unreduced canonical DNF is returned. \code{reduce = TRUE} is interpreted as \code{"rreduce"}, \code{reduce = FALSE} and \code{reduce = NULL} as \code{"none"}.
}
  \item{niter}{
An integer value indicating the number of repetitive applications of \code{\link{rreduce}}. 
\code{niter} will be ignored (with a warning) if \code{reduce} is not equal to \code{"rreduce"}. 
Note that repeated applications may yield identical solutions and that duplicate solutions are eliminiated, so that the number of resulting solutions can be smaller than \code{niter}.
}
  \item{crit}{
Quoted expression specifying a numeric criterion to be maximized when selecting the best solutions among the ones that meet criterion \code{cond}, for example, \code{quote(min(con,cov))} or \code{quote(0.8*con + 0.2*cov)}, etc.}
  \item{cond}{
Quoted expression specifying a logical criterion to be imposed on the solutions inferred from \code{x} before selecting the best solutions on the basis of \code{crit}, for example, \code{quote(con > 0.85)} or \code{quote(con > cov)}, etc.
}
  \item{approx}{
As in \code{\link{conCovOpt}}.  
}
   \item{maxCombs}{
Maximal number of combinations that will be tested for optimality. If the number of necessary iterations exceeds \code{maxCombs}, \code{cnaOpt} will stop executing and return an error message stating the necessary number of iterations. Early termination can then be avoided by increasing \code{maxCombs} accordingly. This argument is passed to \code{\link{conCovOpt}} and \code{\link{ereduce}}.
}

}

\details{
\code{cnaOpt} implements a procedure introduced in Baumgartner and Ambuehl (2021). It infers causal models (atomic solution formulas, asf) for the \code{outcome} from data \code{x} that comply with the logical condition \code{cond} and maximize the numeric criterion \code{crit}. Data \code{x} may be crisp-set (\code{"cs"}) or multi-value (\code{"mv"}), but not fuzzy-set (\code{"fs"}). The function proceeds as follows:
\enumerate{
\item it calculates consistency and coverage optima (con-cov optima) for \code{x};
\item it selects the optima that meet \code{cond};
\item among those optima, it selects those that maximize \code{crit};
\item it builds the canonical disjunctive normal forms (DNF) of the selected optima
\item it generates all minimal forms of those canonical DNFs
}
%first calculates consistency and coverage optima (con-cov optima) for \code{x}, second, it selects the optima that meet \code{cond}, it, third, selects the optima that fare best according to \code{crit}, fourth, it builds the canonical disjunctive normal form (DNF) realizing the best optimum and, finally, generates all minimal forms of that canonical DNF. 
Roughly speaking, running \code{cnaOpt} amounts to sequentially executing \code{\link{configTable}}, \code{\link{conCovOpt}}, \code{\link{selectMax}}, \code{\link{DNFbuild}} and \code{\link{condTbl}}.

In the default setting, \code{cnaOpt} attempts to build all optimal solutions using \code{\link{ereduce}}. But that may be too computationally demanding because the space of optimal solutions can be very large. If the argument \code{reduce} is set to \code{"rreduce"}, \code{cnaOpt} builds one arbitrarily selected optimal solution, which typically terminates quickly. By giving the argument \code{niter} a non-default value, say, 20, the process of selecting one optimal solution under \code{reduce = "rreduce"} is repeated 20 times. As the same solutions will be generated on some iterations and duplicates are not returned, the output may contain less models than the value given to \code{niter}. If \code{reduce} is not set to \code{"rreduce"}, \code{niter} is ignored with a warning.

}

\value{
\code{cnaOpt} returns a \code{\link{data.frame}} with additional classes "cnaOpt" and "condTbl". See the "Value" section in \code{?\link{condTbl}} for details.
}

\seealso{
\code{\link{cna}}, \code{\link{conCovOpt}}
}

\references{
Baumgartner, Michael and Mathias Ambuehl. 2021. \dQuote{Optimizing Consistency and Coverage in Configurational Causal Modeling.} \emph{Sociological Methods & Research}.\cr doi:10.1177/0049124121995554. 
}

\examples{
# Example 1: Real-life crisp-set data, d.educate.
(res_opt1 <- cnaOpt(d.educate, "E"))

# Using the pipe operator (\%>\%), the steps processed by cnaOpt in the 
# call above can be reproduced as follows:
library(dplyr)
conCovOpt(d.educate, "E") \%>\% selectMax \%>\% DNFbuild(reduce = "ereduce") \%>\% 
  paste("<-> E") \%>\% condTbl(d.educate)

# Example 2: Simulated crisp-set data.
dat1 <- data.frame(
  A = c(1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0), 
  B = c(0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0), 
  C = c(0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0), 
  D = c(1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1), 
  E = c(1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1), 
  F = c(0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1)
)

(res_opt2 <- cnaOpt(dat1, "E"))

# Change the maximality criterion.
cnaOpt(dat1, "E", crit = quote(min(con, cov)))
# Change the selection condition.
cnaOpt(dat1, "E", cond = quote(con >= 0.9))
# Build all con-cov optima with coverage above 0.9 that maximize min(con, cov).
cnaOpt(dat1, "E", crit = quote(min(con, cov)), cond = quote(cov > 0.9))
# Different values of the reduce argument.
cnaOpt(dat1, "E", reduce = "none") # canonical DNF
cnaOpt(dat1, "E", reduce = "rreduce") # one randomly drawn optimal solution
# Iterate random solution generation 10 times.
cnaOpt(dat1, "E", reduce = "rreduce", niter = 10) 

# Example 3: All logically possible configurations.
(res_opt3 <- cnaOpt(full.ct(4), "D"))  # All combinations are equally bad.

# Example 4: Real-life multi-value data, d.pban.
cnaOpt(d.pban, outcome = "PB=1")
cnaOpt(d.pban, outcome = "PB=1", crit = quote(0.8*con + 0.2*cov))
cnaOpt(d.pban, outcome = "PB=1", cond = quote(con > 0.9))
cnaOpt(d.pban, outcome = "PB=0")
cnaOpt(d.pban, outcome = "PB=0", cond = quote(con > 0.9))
cnaOpt(d.pban, outcome = "F=2")
cnaOpt(d.pban, outcome = "F=2", crit = quote(0.8*con + 0.2*cov))

# Example 5: High computational demand.
dat2 <- configTable(d.performance[,1:8], frequency = d.performance$frequency)
try(cnaOpt(dat2, outcome = "SP"))   # error because too computationally demanding
# The following call does not terminate because of reduce = "ereduce".
try(cnaOpt(dat2, outcome = "SP", approx = TRUE))
# We could increase maxCombs, as in the line below
\dontrun{cnaOpt(dat2, outcome = "SP", approx = TRUE, maxCombs = 1.08e+09) }
# but this takes very long to terminate.
# Alternative approach: Produce one (randomly selected) optimal solution using reduce = "rreduce".
cnaOpt(dat2, outcome = "SP",  approx = TRUE, reduce = "rreduce")
# Iterate the previous call 10 times.
\donttest{cnaOpt(dat2, outcome = "SP", approx = TRUE, reduce = "rreduce", niter = 10)}
# Another alternative: Use ereduce for minimization but introduce a case.cutoff.
cnaOpt(dat2, outcome = "SP", case.cutoff = 10)
}
