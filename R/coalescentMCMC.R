## coalescentMCMC.R (2013-07-19)

##   Run MCMC for Coalescent Trees

## Copyright 2012-2013 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

coalescentMCMC <-
    function(x, ntrees = 3000, burnin = 1000, frequency = 1,
             tree0 = NULL, quiet = FALSE)
{
    if (is.null(tree0)) {
        d <- dist.dna(x, "JC69")
        X <- phangorn::phyDat(x)
        tree0 <- as.phylo(hclust(d, "average"))
    }

    n <- Ntip(tree0)
    nodeMax <- 2*n - 1
    nOut <- ntrees# + burnin

    getlogLik <- function(phy, X) phangorn::pml(phy, X)$logLik #phangorn:::pml6(phy, X)

    TREES <- vector("list", nOut)
    LL <- numeric(nOut)
    TREES[[1L]] <- tree0
    lnL0 <- getlogLik(tree0, X)
    LL[[1L]] <- lnL0

    i <- 2L
    j <- 0L # number of accepted trees
    k <- 0L # number of sampled trees

    if (!quiet) {
        cat("Running the Markov chain:\n")
        cat("  Number of trees to output:", ntrees, "\n")
        cat("  Burn-in period:", burnin, "\n")
        cat("  Sampling frequency:", frequency, "\n")
        cat("  Number of generations to run:", ntrees * frequency + burnin, "\n")
        cat("Generation    Nb of accepted trees\n")
    }

    while (k < nOut) {
        if (!quiet)
            cat("\r  ", i, "                ", j, "           ")

        tr.b <- NeighborhoodRearrangement(tree0, n, nodeMax)
        if (!(i %% frequency) && i > burnin) {
            k <- k + 1L
            TREES[[k]] <- tr.b
        }
        ## do TipInterchange() every 10 steps:
        ## tr.b <-
        ##     if (!i %% 10) TipInterchange(TREES[[i]], n)
        ##     else NeighborhoodRearrangement(TREES[[i]], n, nodeMax)
        lnL.b <- getlogLik(tr.b, X)
        LL[[i]] <- lnL.b
        i <- i + 1L
        ACCEPT <- if (is.na(lnL.b)) FALSE else {
            if (lnL.b >= lnL0) TRUE
            else rbinom(1, 1, exp(lnL.b - lnL0))
        }
        if (ACCEPT) {
            j <- j + 1L
            lnL0 <- lnL.b
            tree0 <- tr.b
        }
    }

    dim(LL) <- c(i - 1, 1)
    colnames(LL) <- "logLik"
    LL <- mcmc(LL, start = 1, end = i - 1)

    class(TREES) <- "multiPhylo"
    assign(".TREES", TREES, envir = .coalescentMCMCenv)
    ## TREES <- .compressTipLabel(TREES)
    if (!quiet) cat("\nDone.\n")
    #list(mcmc = LL, trees = TREES)
    LL
}

getMCMCtrees <- function() get(".TREES", envir = .coalescentMCMCenv)
