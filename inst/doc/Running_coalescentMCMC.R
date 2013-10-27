### R code from vignette source 'Running_coalescentMCMC.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Running_coalescentMCMC.Rnw:19-20
###################################################
options(width=60)


###################################################
### code chunk number 2: Running_coalescentMCMC.Rnw:68-70
###################################################
library(coalescentMCMC)
body(coalescentMCMC)


###################################################
### code chunk number 3: Running_coalescentMCMC.Rnw:87-89
###################################################
data(woodmouse)
out <- coalescentMCMC(woodmouse, ntrees = 300, burnin = 100, quiet = TRUE)


###################################################
### code chunk number 4: Running_coalescentMCMC.Rnw:112-113
###################################################
plot(out)


###################################################
### code chunk number 5: Running_coalescentMCMC.Rnw:118-120
###################################################
TR <- getMCMCtrees()
TR


###################################################
### code chunk number 6: Running_coalescentMCMC.Rnw:124-126
###################################################
dim(out)
colnames(out)


###################################################
### code chunk number 7: Running_coalescentMCMC.Rnw:131-132
###################################################
out2 <- coalescentMCMC(woodmouse, ntrees = 300, burnin = 100, model = "time", quiet = TRUE)


###################################################
### code chunk number 8: Running_coalescentMCMC.Rnw:151-152
###################################################
plot(out2)


###################################################
### code chunk number 9: Running_coalescentMCMC.Rnw:155-157
###################################################
dim(out2)
colnames(out2)


###################################################
### code chunk number 10: Running_coalescentMCMC.Rnw:172-173
###################################################
TR2 <- get("TREES_2", envir = .coalescentMCMCenv)


###################################################
### code chunk number 11: Running_coalescentMCMC.Rnw:185-189
###################################################
tr <- TR[-(1:200)]
THETA <- out[-(1:300), 2]
logLik0 <- mapply(dcoal, phy = tr, theta = THETA, log = TRUE)
summary(logLik0)


###################################################
### code chunk number 12: Running_coalescentMCMC.Rnw:193-198
###################################################
tr2 <- TR2[-(1:200)]
THETA0 <- out2[-(1:300), 2]
RHO <- out2[-(1:300), 3]
logLik1 <- mapply(dcoal.time, phy = tr2, theta = THETA0, rho = RHO, log = TRUE)
summary(logLik1)


###################################################
### code chunk number 13: Running_coalescentMCMC.Rnw:202-203
###################################################
print(histogram(~c(logLik0, logLik1) | gl(2, 100, labels = c("H0", "H1"))))


###################################################
### code chunk number 14: Running_coalescentMCMC.Rnw:209-212
###################################################
treeloglik <- out2[-(1:300), 1]
(theta0ML <- weighted.mean(THETA0, treeloglik))
(rhoML <- weighted.mean(RHO, treeloglik))


###################################################
### code chunk number 15: Running_coalescentMCMC.Rnw:216-219
###################################################
x <- seq(0, 0.01, 0.0001)
y <- theta0ML * exp(rhoML * x)
plot(-x, y, "l", xlab = "Time", ylab = expression(Theta))


