### R code from vignette source 'Running_coalescentMCMC.Rnw'

###################################################
### code chunk number 1: Running_coalescentMCMC.Rnw:19-20
###################################################
options(width=60)


###################################################
### code chunk number 2: Running_coalescentMCMC.Rnw:56-58
###################################################
library(coalescentMCMC)
args(coalescentMCMC)


###################################################
### code chunk number 3: Running_coalescentMCMC.Rnw:79-81
###################################################
data(woodmouse)
out <- coalescentMCMC(woodmouse, ntrees = 300, printevery = 0)


###################################################
### code chunk number 4: Running_coalescentMCMC.Rnw:101-102
###################################################
plot(out)


###################################################
### code chunk number 5: Running_coalescentMCMC.Rnw:107-109
###################################################
TR <- getMCMCtrees()
TR


###################################################
### code chunk number 6: Running_coalescentMCMC.Rnw:113-115
###################################################
dim(out)
colnames(out)


###################################################
### code chunk number 7: Running_coalescentMCMC.Rnw:120-121
###################################################
out2 <- coalescentMCMC(woodmouse, ntrees = 300, model = "time", printevery = 0)


###################################################
### code chunk number 8: Running_coalescentMCMC.Rnw:137-138
###################################################
plot(out2)


###################################################
### code chunk number 9: Running_coalescentMCMC.Rnw:141-143
###################################################
dim(out2)
colnames(out2)


###################################################
### code chunk number 10: Running_coalescentMCMC.Rnw:158-159
###################################################
TR2 <- getMCMCtrees(2)


###################################################
### code chunk number 11: Running_coalescentMCMC.Rnw:164-165
###################################################
getMCMCstats()


