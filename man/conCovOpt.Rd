\name{conCovOpt}
\alias{conCovOpt}
\alias{print.conCovOpt}
\alias{plot.conCovOpt}

\title{
Find consistency and coverage optima for configurational data
}

\description{
\code{conCovOpt} issues pairs of optimal consistency and coverage scores that atomic solution formulas (asf) of an outcome inferred from configurational data can possibly reach (cf. Baumgartner and Ambuehl 2021).
}

\usage{
conCovOpt(x, outcome = NULL, ..., rm.dup.factors = FALSE, rm.const.factors = FALSE,
          maxCombs = 1e+07, approx = FALSE, allConCov)
\method{print}{conCovOpt}(x, ...)
\method{plot}{conCovOpt}(x, con = 1, cov = 1, ...)
}

\arguments{
  \item{x}{
In \code{conCovOpt}: a \code{\link{data.frame}} or \code{\link{configTable}}. In the \code{print}- and \code{plot}-method: an output of \code{conCovOpt}.
}
  \item{outcome}{
A character vector of one or several factor values in \code{x}.
}
 \item{\dots}{
In \code{conCovOpt}: arguments passed to \code{\link{configTable}}, e.g. \code{case.cutoff}. 
 The \sQuote{\code{\dots}} are currently not used in \code{plot.conCovOpt}.
}
  \item{rm.dup.factors}{Logical; defaults to \code{FALSE} (which is different from \code{\link{configTable}}). If \code{TRUE}, all but the first of a set of factors with identical values in \code{x} are removed.}
  \item{rm.const.factors}{Logical;  defaults to \code{FALSE} (which is different from \code{\link{configTable}}). If \code{TRUE}, factors with constant values in \code{x} are removed.}
 
  \item{maxCombs}{
Maximal number of combinations that will be tested for optimality. If the number of necessary iterations exceeds \code{maxCombs}, \code{conCovOpt} will stop executing and return an error message stating the necessary number of iterations. Early termination can then be avoided by increasing \code{maxCombs} accordingly.
}
  \item{approx}{
Logical; if TRUE, an exhaustive search is only approximated; if FALSE, an exhaustive search is conducted. 
}
  \item{allConCov}{
Defunct argument (as of package version 0.5.0). See the remark in \code{?\link{multipleMax}}.
}
  \item{con, cov}{Numeric scalars between 0 and 1 indicating consistency and coverage thresholds marking the area of "good" models in a square drawn in the plot. Points within the square correspond to models reaching these thresholds.
}
}

\details{
\code{conCovOpt} implements a procedure introduced in Baumgartner and Ambuehl (2021). It  calculates consistency and coverage optima for models (i.e. atomic solution formulas, asf) of an \code{outcome} inferred from data \code{x} prior to actual CNA or QCA analyses.

An ordered pair (con, cov) of consistency and coverage scores is a \strong{con-cov optimum} for outcome Y=k in data \code{x} iff it is not excluded (based e.g. on the data structure) for an asf of Y=k inferred from \code{x} to reach (con, cov) but excluded to score better on one element of the pair and at least as well on the other.

\code{conCovOpt} calculates con-cov optima by executing the following steps:\enumerate{\item if \code{x} is a data frame, aggregate \code{x} in a \code{configTable}, \item build exo-groups with constant values in all factors other than the \code{outcome}, \item assign output values to each exo-group that reproduce the behavior of \code{outcome} as closely as possible, \item  calculate con-cov scores for each assignment resulting in step 3, \item eliminate all non-optimal scores.}

The implementation of step 4 calculates con-cov scores of about 10 million output value assignments in reasonable time, but step 3 may result in considerably more assignments. In such cases, the argument \code{approx} may be set to its non-default value \code{"TRUE"}, which determines that step 4 is only executed for those assignments closest to the \code{outcome}'s median value. This is an efficient approach for finding many, but possibly not all, con-cov optima.

In case of crisp-set and multi-value data, at least one actual model (asf) inferrable from \code{x} and reaching an optimum's consistency and coverage scores is guaranteed to exist for every con-cov optimum. The function \code{\link{DNFbuild}} can be used to build these optimal models. The same does not hold for fuzzy-set data. In fuzzy-set data it merely holds that the existence of a model reaching an optimum's consistency and coverage scores cannot be excluded prior to an actual application of \code{\link{cna}}. 

}

\value{
An object of class 'conCovOpt'. The exo-groups resulting from step 2 are stored as attribute \code{"exoGroups"}, the lists of output values resulting from step 3 are stored as attribute \code{"reprodList"} (reproduction list).
}

\references{
Baumgartner, Michael and Mathias Ambuehl. 2021. \dQuote{Optimizing Consistency and Coverage in Configurational Causal Modeling.} \emph{Sociological Methods & Research}.\cr doi:10.1177/0049124121995554. 
}



\seealso{
\code{\link{configTable}}, \code{\link{selectMax}}, \code{\link{DNFbuild}}
}

\examples{
(cco.irrigate <- conCovOpt(d.irrigate))
conCovOpt(d.irrigate, outcome = c("R","W"))
# Plot method.
plot(cco.irrigate)
plot(cco.irrigate, con = .8, cov = .8)

dat1 <- d.autonomy[15:30, c("EM","SP","CO","AU")]
(cco1 <- conCovOpt(dat1, outcome = "AU"))

print(cco1, digits = 3, row.names = TRUE)
plot(cco1)

# Exo-groups (configurations with constant values in all factors other than the outcome).
attr(cco1$A, "exoGroups")

# Rep-list (list of values optimally reproducing the outcome).
attr(cco1$A, "reprodList")

dat2 <- d.pacts
# Maximal number of combinations exceeds maxCombs.
(cco2 <- conCovOpt(dat2, outcome = "PACT")) # Generates a warning
# Increase maxCombs.
\donttest{(cco2_full <- try(conCovOpt(dat2, outcome = "PACT", 
  maxCombs=1e+08))) # Takes a long time to terminate}
# Approximate an exhaustive search.
(cco2_approx1 <- conCovOpt(dat2, outcome = "PACT", approx = TRUE))
selectMax(cco2_approx1)
# The search space can also be reduced by means of a case cutoff.
(cco2_approx2 <- conCovOpt(dat2, outcome = "PACT", case.cutoff=2))
selectMax(cco2_approx2)

}
