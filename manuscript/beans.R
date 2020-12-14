library("data.table")
library("beanplot")

setA <- fread("5times6subsinbs3single.csv")
setB <- fread("10times10subsinbs5single.csv")
setC <- fread("67889_3single.csv")
setD <- fread("5.5.8.8.10.10.12.12.15.15bs5single.csv")

sapply(setA, range)
tail(sapply(setA, function (x) table(round(x))))
sum(setB$randBin > min(setB$SBA))
sum(setB$randBin >= min(setB$SBA))
sapply(setC, range)
sapply(setD, range)

setAm <- fread("5times6subsinbs3maratime.csv")
setBm <- fread("10times10subsinbs5maratime.csv")
setCm <- fread("67889_3maratime.csv")
setDm <- fread("5.5.8.8.10.10.12.12.15.15bs5maratime.csv")

sapply(setAm, range)
sapply(setBm, range)
sum(setBm$randBin > min(setBm$SBA))
sapply(setCm, range)
sapply(setCm, function (x) table(round(x)))
sapply(setDm, range)
sum(setDm$RBA >= min(setDm$SBA))
sum(setDm$RBA >= median(setDm$SBA))

## colourscheme from
## https://colorbrewer2.org/#type=diverging&scheme=PuOr&n=4
PuOr1 <- "#e66101"
PuOr2 <- "#fdb863"
PuOr3 <- "#b2abd2"
PuOr4 <- "#5e3c99"
png("beanplots.png", width=800, height=800)
par(mfrow=c(2,2))
beanplot(setA[, c("RBA", "SBA")],
         what=c(0,1,0,1),
         col = list(c(PuOr4, PuOr3), c(PuOr1, PuOr2)),
         ll=0.0005,
         method = "stack", innerborder=NA,
         main = "Balanced incomplete-block",
         ylab="D-criterium", names=c("RBA", "SBA"))
mtext("A", side=3, line=1.75, at=0.15)
mtext("1000 Allocations", side=1, line=3, at=1.45)
beanplot(cbind(setB[,c("RBA", "SBA")], setBm[,c("RBA", "SBA")]),
         what=c(0,1,0,1),
         col = list(c(PuOr4, PuOr3, PuOr3),
                    c(PuOr1, PuOr2, PuOr2)),
         ll=0.005,
         method = "stack", innerborder=NA,
         main = "Equi-replicate incomplete-block",
         ylab="D-criterium", names=c("RBA", "SBA", "RBA", "SBA"))
mtext("B", side=3, line=1.75, at=-0.3)
mtext("1000 Allocations", side=1, line=3, at=1.475)
mtext("1000 Runs", side=1, line=3, at=3.5)
beanplot(cbind(setC[,c("RBA", "SBA")], setCm[,c("RBA", "SBA")]),
         what=c(0,1,0,1),
         col = list(c(PuOr4, PuOr3, PuOr3),
                    c(PuOr1, PuOr2, PuOr2)),
         ll=0.005,
         method = "stack", innerborder=NA,
         main = "Incomplete-block not equi-replicate",
         ylab="D-criterium", names=c("RBA", "SBA", "RBA", "SBA"))
mtext("C", side=3, line=1.75, at=-0.2)
mtext("1000 Allocations", side=1, line=3, at=1.475)
mtext("1000 Runs", side=1, line=3, at=3.5)
beanplot(cbind(setD[,c("RBA", "SBA")], setDm[,c("RBA", "SBA")]),
         what=c(0,1,0,1),
         col = list(c(PuOr4, PuOr3, PuOr3),
                    c(PuOr1, PuOr2, PuOr2)),
         ll=0.005,
         method = "stack", innerborder=NA,
         main = "Incomplete-block not equi-replicate",
         ylab="D-criterium", names=c("RBA", "SBA", "RBA", "SBA"))
mtext("D", side=3, line=1.75, at=-0.3)
mtext("1000 Allocations", side=1, line=3, at=1.475)
mtext("1000 Runs", side=1, line=3, at=3.5)
dev.off()
