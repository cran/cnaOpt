\name{ereduce}
\alias{ereduce}

\title{Find all minimal disjunctive normal forms (DNF) of an input DNF}

\description{\code{ereduce} builds all minimal disjunctive normal forms corresponding to an input DNF. It is similar to \code{\link[cna]{rreduce}}, which, however, only builds one minimal DNF at random.
}

\usage{
ereduce(cond, x = full.ct(cond), full = !missing(x), 
        simplify2constant = TRUE, maxCombs = 1e7)
}

\arguments{
  \item{cond}{
A character string specifying a disjunctive normal form (DNF); can be either crisp-set or multi-value.
}
  \item{x}{A 
\code{configTable} or \code{data.frame}; can be either crisp-set or multi-value.
}
  \item{full}{
Logical; if \code{TRUE}, redundancies are eliminated relative to \code{full.ct(x)}, otherwise relative to \code{x}.
}
  \item{simplify2constant}{
Logical; if \code{TRUE} (the default), a tautologous or contradictory \code{cond} is reduced to a constant \code{"1"} or \code{"0"}, respectively. If \code{FALSE}, a minimal tautology or contradiction, i.e. \code{"A+a"} or \code{"A*a"}, will result.
}
  \item{maxCombs}{
Maximal number of iterations that will be ran in the most time-consuming step. If the number of necessary iterations exceeds \code{maxCombs}, \code{ereduce} will stop executing and return an error message stating the necessary number of iterations. Early termination can then be avoided by increasing \code{maxCombs} accordingly.

}
}

\details{
\code{ereduce} eliminates conjuncts and disjuncts from a DNF \code{cond} as long as the result of \code{condition(cond, x)} remains the same. The only required argument is \code{cond}. If \code{x} is not provided, redundancies are eliminated relative to \code{full.ct(cond)}. 

\code{ereduce} generates all redundancy-free forms of \code{cond}, while \code{\link[cna]{rreduce}} only returns one randomly chosen one. \code{rreduce} is faster than \code{ereduce}, but often incomplete. In a nutshell,  \code{ereduce} searches for minimal hitting sets in \code{cond} preventing \code{cond} from being false in data \code{x}.
}

\value{
A vector of redundancy-free disjunctive normal forms (DNF).
}

\seealso{
\code{\link[cna]{rreduce}}, \code{\link{full.ct}}, \code{\link{conCovOpt}}, \code{\link{DNFbuild}}.
}


\examples{
# Logical redundancies.
cond1 <- "A*b + a*B + A*C + B*C"
ereduce(cond1)
rreduce(cond1) # repeated calls generate different outputs
cond2 <- "A*b + a*B + A*B + a*b"
ereduce(cond2)
ereduce(cond2, simplify2constant = FALSE)

# Redundancy elimination relative to simulated cs data.
dat1 <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 0, 1), 
  B = c(0, 1, 0, 1, 1, 0, 0, 0), 
  C = c(1, 1, 0, 1, 1, 0, 1, 1), 
  D = c(0, 0, 0, 0, 0, 1, 1, 1))
cco1 <- conCovOpt(dat1, "D")
best1 <- selectMax(cco1)
(formula1 <- DNFbuild(best1, outcome = "D", reduce = FALSE))
# ereduce
ereduce(formula1, dat1, full = FALSE)
# rreduce
rreduce(formula1, dat1, full = FALSE)

# Redundancy elimination relative to simulated mv data.
dat2 <- data.frame(
  A = c(3,2,1,1,2,3,2,2,2,1,1,2,3,2,2,2,1,2,3,3,3,1,1,1,3,1,2,1,2,3,3,2,2,2,1,2,2,3,2,1,2,1,3,3),
  B = c(1,2,3,2,1,1,2,1,2,2,3,1,1,1,2,3,1,3,3,3,1,1,3,2,2,1,1,3,3,2,3,1,2,1,2,2,1,1,2,2,3,3,3,3),
  C = c(1,3,3,3,1,1,1,2,2,3,3,1,1,2,2,2,3,1,1,2,1,2,2,3,3,1,2,2,2,3,2,1,1,2,2,2,1,1,1,2,2,1,1,2),
  D = c(3,1,2,2,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,1,1,1,1,1,2,2,2,2,2,3,1,1,1,1,1,2,2,2,2,2,3,3,3),
  E = c(3,2,2,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3)
)
cco2 <- conCovOpt(dat2, "D=3")
best2 <- selectMax(cco2)
(formula2 <- DNFbuild(best2, outcome = "D=3", reduce = FALSE))
# ereduce
ereduce(formula2, dat2, full = FALSE)
# rreduce
rreduce(formula2, dat2, full = FALSE)

# Any Boolean expressions.
cond <- "!(A*B*C)*!(a*b*c)" # or "A + B*!(D + e) <-> C" 
x <- selectCases(cond) 
(cond <- cna:::getCond(x)) # returns a DNF equivalent to cond, but with many redundancies
ereduce(cond)
rreduce(cond)
}
