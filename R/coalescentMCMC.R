## coalescentMCMC.R (2012-09-26)

##   Run MCMC for Coalescent Trees

## Copyright 2012 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

coalescentMCMC <-
    function(x, ntrees = 1e4, burnin = ntrees,
             tree0 = NULL, quiet = FALSE)
{
    if (is.null(tree0)) {
        d <- dist.dna(x, "JC69")
        X <- phangorn::phyDat(x)
        tree0 <- as.phylo(hclust(d, "average"))
    }

    n <- Ntip(tree0)
    nodeMax <- 2*n - 1

    TREES <- vector("list", ntrees)
    LL <- numeric(ntrees)

    lnL0 <- phangorn::pml(tree0, X)$logLik

    i <- j <- 0L # j: number of tree proposals

    if (!quiet) {
        cat("Running the Markov chain:\n")
        cat("Nb of proposed trees    Nb of accepted trees\n")
    }

    while (i <= ntrees) {
        if (!quiet)
            cat("\r     ", j, "                   ", i, "           ")
        j <- j + 1L
        tr.b <- NeighborhoodRearrangement(tree0, n, nodeMax)
        ## do TipInterchange() every 10 steps:
        ## tr.b <-
        ##     if (!i %% 10) TipInterchange(TREES[[i]], n)
        ##     else NeighborhoodRearrangement(TREES[[i]], n, nodeMax)
### see the note about CladeInterchange() in 'dcoal.R'
        lnL.b <- phangorn::pml(tr.b, X)$logLik
        ACCEPT <- if (is.na(lnL.b)) FALSE else {
            if (lnL.b >= lnL0) TRUE
            else rbinom(1, 1, exp(lnL.b - lnL0))
        }
        if (ACCEPT) {
            lnL0 <- lnL.b
            tree0 <- tr.b
            if (j > burnin) {
                i <- i + 1L
                LL[i] <- lnL0
                TREES[[i]] <- tree0
            }
        }
    }

    class(TREES) <- "multiPhylo"
    ## TREES <- .compressTipLabel(TREES)
    if (!quiet) cat("\nDone.\n")
    list(trees = TREES, logLik = LL)
}
