## -----------------------------------------------------------------------------
library("equateMultiple")
data("data2pl", package = "equateIRT")

## ----message=FALSE, results='hide'--------------------------------------------
library("mirt")
m1 <- mirt(data2pl[[1]], SE = TRUE)
m2 <- mirt(data2pl[[2]], SE = TRUE)
m3 <- mirt(data2pl[[3]], SE = TRUE)
m4 <- mirt(data2pl[[4]], SE = TRUE)
m5 <- mirt(data2pl[[5]], SE = TRUE)

## -----------------------------------------------------------------------------
mlist<- list(m1, m2, m3, m4, m5)
test <- paste("test", 1:5, sep = "")
mods <- modIRT(est.mods = mlist, names = test, display = FALSE)

## -----------------------------------------------------------------------------
lplan<-linkp(mods = mods)
lplan

## -----------------------------------------------------------------------------
eqMM <- multiec(mods = mods, method = "mean-mean")
summary(eqMM)

## -----------------------------------------------------------------------------
eqMGM <- multiec(mods = mods, method = "mean-gmean")
summary(eqMGM)

## -----------------------------------------------------------------------------
eqIRF<-multiec(mods = mods, method = "irf")
summary(eqIRF)

## -----------------------------------------------------------------------------
eqMGM <- multiec(mods = mods, method = "mean-gmean", se = FALSE)
eqIRF<-multiec(mods = mods, method = "irf", start = eqMGM)
summary(eqIRF)

## -----------------------------------------------------------------------------
eqTRF<-multiec(mods = mods, method = "trf")
summary(eqTRF)

## -----------------------------------------------------------------------------
eqLIK <- multiec(mods = mods, method = "lik")
summary(eqLIK)

## -----------------------------------------------------------------------------
eqLIK <- multiec(mods = mods, method = "lik", base = 5)
summary(eqLIK)

## -----------------------------------------------------------------------------
item.common(eqIRF)

## -----------------------------------------------------------------------------
scTSE<-score(eqIRF)
round(scTSE,3)

## -----------------------------------------------------------------------------
scOSE<-score(eqIRF, method = "OSE", se = FALSE)
round(scOSE,3)

