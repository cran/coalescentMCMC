## treeOperators.R (2012-09-19)

##   Trees Operators for Running MCMC

## Copyright 2012 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

NeighborhoodRearrangement <- function(phy, n, nodeMax)
{
    ##phy <- tr <- rcoal(n <- 10); nodeMax <- 2*n - 1
    ##plot(tr); nodelabels(); axisPhylo()

    bt <- c(rep(0, n), branching.times(phy))
    target <- sample((n + 2):nodeMax, size = 1)
    i1 <- which(phy$edge[, 2] == target)
    parent <- phy$edge[i1, 1]
    i2 <- which(phy$edge[, 1] == target) # the 2 descendants from 'target'
    i3 <- which(phy$edge[, 1] == parent) # this includes i1, so:
    i3 <- i3[!i3 %in% i1]
    sister <- phy$edge[i3, 2]
    select <- sample(c(TRUE, FALSE))
    i2.move <- i2[select]
    i2.stay <- i2[!select]
    phy$edge[i3, 2] <- child2move <- phy$edge[i2.move, 2]
    child2stay <- phy$edge[i2.stay, 2]
    phy$edge[i2.move, 2] <- sister
    ## now adjust branch lengths:
    ## adjust the branch length that was subtending 'sister':
    phy$edge.length[i3] <- bt[parent] - bt[child2move]
    ## random age for 'target' between the ones of 'sister' and 'parent':
    ## newage <- runif(1, bt[sister], bt[parent])
    newage <- (bt[parent] + max(bt[sister], bt[child2stay]))/2
    phy$edge.length[i1] <- bt[parent] - newage
    phy$edge.length[i2.move] <- newage - bt[sister]
    ## adjust the branch length below the child that has not been moved:
    phy$edge.length[i2.stay] <- newage - bt[child2stay]
    attr(phy, "order") <- NULL
    phy <- reorder(phy)
    newNb <- integer(nodeMax)
    newNb[n + 1L] <- n + 1L
    sndcol <- phy$edge[, 2] > n
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <- (n + 2):nodeMax
    phy$edge[, 1] <- newNb[phy$edge[, 1]]
    phy
    ##plot(phy); nodelabels(); axisPhylo()
}

TipInterchange <- function(phy, n)
{
    getIndexEdge <- function(tip, edge)
        ## 'integer(1)' mustn't be substituted by '0L' except if 'DUP = TRUE':
        .C("get_single_index_integer", as.integer(edge[, 2L]),
           as.integer(tip), integer(1L),
           PACKAGE = "coalescentMCMC", NAOK = TRUE, DUP = FALSE)[[3L]]

    repeat {
        chosen <- sample(n, size = 2L)
        i1 <- getIndexEdge(chosen[1], phy$edge)
        i2 <- getIndexEdge(chosen[2], phy$edge)
        ## check that the two tips in 'chosen' are not sisters
        if (phy$edge[i1, 1L] != phy$edge[i2, 1L]) break
    }

    ##phy$edge[c(i2, i1), 2] <- chosen
    phy$tip.label[chosen] <- phy$tip.label[rev(chosen)]
    phy
}

EdgeLengthJittering <- function(phy)
### all edge lengths are added to a random value on U[-MIN, MAX]
### (the ultrametric nature of the tree is kept)
{
   z <- range(phy$edge.length)
   MIN <- z[1]
   MAX <- z[2]
   x <- runif(1, -MIN, MAX) # should be OK even if MIN=0
   phy$edge.length <- phy$edge.length + x
   phy
}
