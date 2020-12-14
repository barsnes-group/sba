
library("data.table")

## include("contrastvariances.jl")
## contrastvariances()

setC <- fread("67889_3contrastvar.csv")
colnames(setC) <- c("fun", "contr", "variance")
unique(setC)
unique(setC[order(contr, fun)])
for (xcontr in unique(setC$contr)) {
    for (xfun in unique(setC$fun)) {
        print(paste(xfun, xcontr))
        print(table(setC[contr == xcontr & fun == xfun, variance]))
    }
}

merge(setC[fun=="rba", .(minim=min(variance)), by=.(contr)],
      setC[fun=="sba", .(maxim=max(variance)), by=.(contr)])[, maxim/minim]
merge(setC[fun=="rba", .(maxim=max(variance)), by=.(contr)],
      setC[fun=="sba", .(minim=min(variance)), by=.(contr)])[, minim/maxim]
