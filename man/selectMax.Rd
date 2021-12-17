\name{selectMax}
\alias{selectMax}
\alias{multipleMax}

\title{
Select the con-cov optima from a '\code{conCovOpt}' object that maximize a specified optimality criterion
}

\description{
\code{selectMax} selects the optima from a '\code{conCovOpt}' object that maximize a specified optimality criterion (cf. Baumgartner and Ambuehl 2021). 
}

\usage{
selectMax(x, crit = quote(con * cov), cond = quote(TRUE), warn = TRUE)
multipleMax(x, outcome)
}

\arguments{
  \item{x}{
An object output by \code{\link{conCovOpt}}.
}
  \item{crit}{
Quoted expression specifying a numeric criterion to be maximized when selecting from the con-cov optima that meet criterion \code{cond}, for example, \code{min(con,cov)} or \code{0.8*con + 0.2*cov}, etc.}
  \item{cond}{ 
Quoted expression specifying a logical criterion to be imposed on the con-cov optima in \code{x} before selecting the optima maximizing \code{crit}, for example, \code{con > 0.85} or \code{con > cov}, etc.
}
  \item{warn}{
Logical; if \code{TRUE}, \code{selectMax()} returns a warning if no solution is found.
}
  \item{outcome}{
A character string specifying a single outcome value in the original data.
}
}

\details{
While \code{\link{conCovOpt}} identifies \emph{all} con-cov optima in an analyzed data set, \code{selectMax} selects those optima from a '\code{conCovOpt}'  object \code{x} that comply with a logical condition \code{cond} and fare best according to the numeric optimality criterion \code{crit}. The default is to select so-called \emph{con-cov maxima}, meaning con-cov optima with highest product of consistency and coverage. % : An ordered pair (con, cov) of consistency and coverage scores is a \strong{con-cov maximum} for outcome Y=k in data \eqn{\delta} iff (con, cov) is a con-cov optimum for Y=k in \eqn{\delta} with \strong{highest product of consistency and coverage} (con-cov product).
But the argument \code{crit} allows for specifying any other numeric optimality criterion, e.g. \code{min(con, cov)}, \code{max(con, cov)}, or \code{0.8*con + 0.2*cov}, etc. (see Baumgartner and Ambuehl 2021). If \code{x} contains multiple outcomes, the selection of the best con-cov optima is done separately for each outcome.

As of package version 0.5.0, the function \code{multipleMax} is obsolete. It is kept for backwards compatibility only. 

Via the column \code{id} in the output of \code{selectMax} it is possible to select one among many equally good maxima, for instance, by means of \code{\link{reprodAssign}} (see the examples below).
}

%\section{Remark}{
%From version 0.5.0 of the package, \code{selectMax} returns multiple best lines if the highest \code{crit} value is attained in several rows of the '\code{conCovOpt}' object. As a consequence, the function \code{multipleMax} and the argument \code{allConCov} in \code{\link{conCovOpt}} are now obsolete. They are kept in the package to ensure backward compatibility of existing code. 
%}

\value{
\code{selectMax} returns an object of class 'selectMax'. 
}

\seealso{
\code{\link{conCovOpt}}, \code{\link{reprodAssign}}

See also examples in \code{\link{conCovOpt}}.
}

\references{
Baumgartner, Michael and Mathias Ambuehl. 2021. \dQuote{Optimizing Consistency and Coverage in Configurational Causal Modeling.} \emph{Sociological Methods & Research}.\cr doi:10.1177/0049124121995554. 
}

\examples{
dat1 <- d.autonomy[15:30, c("EM","SP","CO","AU")]
(cco1 <- conCovOpt(dat1, outcome = "AU"))
selectMax(cco1)
selectMax(cco1, cond = quote(con > 0.95))
selectMax(cco1, cond = quote(cov > 0.98))
selectMax(cco1, crit = quote(min(con, cov)))
selectMax(cco1, crit = quote(max(con, cov)), cond = quote(cov > 0.9))

# Multiple equally good maxima.
(cco2 <- conCovOpt(dat1, outcome = "AU")) 
(sm2 <- selectMax(cco2, cond = quote(con > 0.93)))
# Each maximum corresponds to a different rep-assignment, which can be selected
# using the id argument.
reprodAssign(sm2, "AU", id = 10)
reprodAssign(sm2, "AU", id = 11)
reprodAssign(sm2, "AU", id = 13)
}