## coalescentMCMC.R (2013-10-18)

##   Run MCMC for Coalescent Trees

## Copyright 2012-2013 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

coalescentMCMC <-
    function(x, ntrees = 3000, burnin = 1000, frequency = 1,
             tree0 = NULL, model = NULL, quiet = FALSE)
{
    if (is.null(tree0)) {
        d <- dist.dna(x, "JC69")
        tree0 <- as.phylo(hclust(d, "average"))
    }

    X <- phyDat(x)
    n <- length(tree0$tip.label)
    nodeMax <- 2*n - 1
    nOut <- ntrees
    nOut2 <- ntrees * frequency + burnin

    getlogLik <- function(phy, X) pml(phy, X)$logLik

    TREES <- vector("list", nOut)
    LL <- numeric(nOut2)
    TREES[[1L]] <- tree0
    lnL0 <- getlogLik(tree0, X)
    LL[1L] <- lnL0

    if (is.null(model)) {
        np <- 1L
        para.nms <- "theta"
        ## quantities to calculate THETA:
        two2n <- 2:n
        K4theta <- length(two2n)
        tmp <- choose(two2n, 2)
        getparams <- function(phy, bt) {
            x4theta <- rev(diff(c(0, sort(bt))))
            sum(x4theta * tmp)/K4theta
        }
        f.theta <- function(t, p) p
    } else { # only "time" (ie, exponential model) for the moment
        np <- 2L
        para.nms <- c("theta0", "rho")
        getparams <- function(phy, bt) { # 'bt' is not used but is needed to have the same arguments than above
            out <- nlminb(c(0.02, 0),
                          function(p) -dcoal.time(phy, p[1], p[2], log = TRUE))
            out$par
        }
        f.theta <- function(t, p) p[1] * exp(p[2] * t)
    }
    params <- matrix(0, nOut2, np)

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

    ##
    bt0 <- branching.times(tree0)
    params[1L, ] <- para0 <- getparams(tree0, bt0)

    nodesToSample <- (n + 2):nodeMax

    while (k < nOut) {
        if (!quiet)
            cat("\r  ", i, "                ", j, "           ")

        ## select one internal node excluding the root:
        target <- sample(nodesToSample, 1L) # target node for rearrangement
        THETA <- f.theta(bt0[target - n], para0) # the value of THETA at this node

        tr.b <- NeighborhoodRearrangement(tree0, n, nodeMax, target, THETA, bt0)
        ## do TipInterchange() every 10 steps:
        ## tr.b <-
        ##     if (!i %% 10) TipInterchange(tree0, n)
        ##     else NeighborhoodRearrangement(tree0, n, nodeMax, target, THETA, bt0)
        if (!(i %% frequency) && i > burnin) {
            k <- k + 1L
            TREES[[k]] <- tr.b
        }
        lnL.b <- getlogLik(tr.b, X)
        LL[i] <- lnL.b
        ## calculate theta for the proposed tree:
        bt <- branching.times(tr.b)
        params[i, ] <- para <- getparams(tr.b, bt)
        i <- i + 1L
        ACCEPT <- if (is.na(lnL.b)) FALSE else {
            if (lnL.b >= lnL0) TRUE
            else rbinom(1, 1, exp(lnL.b - lnL0))
        }
        if (ACCEPT) {
            j <- j + 1L
            lnL0 <- lnL.b
            tree0 <- tr.b
            para0 <- para
            bt0 <- bt
        }
    }

    #dim(LL) <- c(i - 1, 1)
    LL <- cbind(LL, params)
    colnames(LL) <- c("logLik", para.nms)
    LL <- mcmc(LL, start = 1, end = i - 1)

    ## compress the list of trees:
    attr(TREES, "TipLabel") <- TREES[[1L]]$tip.label
    for (i in seq_len(nOut)) TREES[[i]]$tip.label <- NULL
    class(TREES) <- "multiPhylo"

    j <- 1
    list.trees <- ls(envir = .coalescentMCMCenv)
    if (l <- length(list.trees))
        j <- 1 + as.numeric(sub("TREES_", "", list.trees[l]))
    assign(paste("TREES", j, sep = "_"), TREES, envir = .coalescentMCMCenv)
    if (!quiet) cat("\nDone.\n")
    LL
}

getMCMCtrees <- function() {
    list.trees <- ls(envir = .coalescentMCMCenv)
    l <- length(list.trees)
    if (!l) return(NULL)
    if (l == 1)
        return(get(list.trees, envir = .coalescentMCMCenv))
    ## l > 1:
    cat("Several lists of MCMC trees are stored:\n\n")
    for (i in 1:l) cat(i, ":", list.trees[i], "\n")
    cat("\nReturn which number? ")
    i <- as.numeric(readLines(n = 1))
    get(paste("TREES", i, sep = "_"), envir = .coalescentMCMCenv)
}
