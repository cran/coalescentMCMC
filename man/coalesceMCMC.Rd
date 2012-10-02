\name{coalescentMCMC}
\alias{coalescentMCMC}
\title{Run MCMC for Coalescent Trees}
\description{
  This function runs a simple Markov chain Monte Carlo (MCMC) algorithm
  to generate a set of trees which is returned with their likelihoods.
}
\usage{
coalescentMCMC(x, ntrees = 1e4, burnin = ntrees,
               tree0 = NULL, quiet = FALSE)

}
\arguments{
  \item{x}{a set of DNA sequences, typically an object of class
    \code{"DNAbin"} or \code{"phyDat"}.}
  \item{ntrees}{the number of trees to output.}
  \item{burnin}{the number of trees to discard as ``burn-in''.}
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
}
\value{
  a list with the following elements:

  \item{trees}{an object of class \code{"multiPhylo"} containing the
  `ntrees' trees.}
  \item{logLik}{a numerical vector with the log-likelihood of each tree.}
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
system.time(out <- coalescentMCMC(woodmouse)) # circa 2.5 mins
plot(out$logLik, type = "l")
}}
\keyword{models}
