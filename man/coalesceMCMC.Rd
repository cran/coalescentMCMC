\name{coalescentMCMC}
\alias{coalescentMCMC}
\alias{getMCMCtrees}
\title{Run MCMC for Coalescent Trees}
\description{
  This function runs a Markov chain Monte Carlo (MCMC) algorithm to
  generate a set of trees which is returned with their likelihoods.

  \code{getMCMCtrees} extracts the trees from the last MCMC run.
}
\usage{
coalescentMCMC(x, ntrees = 3000, burnin = 1000, frequency = 1,
               tree0 = NULL, quiet = FALSE)
getMCMCtrees()
}
\arguments{
  \item{x}{a set of DNA sequences, typically an object of class
    \code{"DNAbin"} or \code{"phyDat"}.}
  \item{ntrees}{the number of trees to output.}
  \item{burnin}{the number of trees to discard as ``burn-in''.}
  \item{frequency}{the frequency at which trees are sampled.}
  \item{tree0}{the initial tree of the chain; by default, a UPGMA
    tree with a JC69 distance is generated.}
  \item{quiet}{a logical specifying whether to print the numbers of trees
    proposed and accepted.}
}
\details{
  A relatively simple MCMC algorithm is programmed using at each step
  the ``neighborhood rearrangement'' operation (Kuhner et al., 1995) and
  Hastings's ratio for acceptance/rejection of the proposed tree. This
  can be modified easily.

  The number of generations of the chain is determined by: `frequency'
  multiplied by `ntrees' + `burnin'. Only the `ntrees' trees are output
  whereas all the log-likelihood values are output.
}
\value{
  an object of class \code{"coda"} with the log-likelihood of each tree.

  The list of trees is returned in the hidden object \code{.TREES} which
  can be extracted with \code{getMCMCtrees}.
}
\references{
  Hastings, W. K. (1970) Monte Carlo sampling methods using Markov
  chains and their applications. \emph{Biometrika}, \bold{57}, 97--109.

  Kuhner, M. K., Yamato, J. and Felsenstein, J. (1995) Estimating
  effective population size and mutation rate from sequence data using
  Metropolis-Hastings sampling. \emph{Genetics}, \bold{140}, 1421--1430.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{dcoal}}, \code{\link{treeOperators}}
}
\examples{
\dontrun{
data(woodmouse)
system.time(out <- coalescentMCMC(woodmouse)) # circa 12 sec.
plot(out)
getMCMCtrees() # returns 3000 trees
}}
\keyword{models}
