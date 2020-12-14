library("data.table")

## include("spaceresults.jl")
## optimalrun(6, "optimal4andon.csv")
## runspace(6, "mbs4andon.csv")
sizetime <- fread("mbs4andon.csv")
besttime <- fread("optimal4andon.csv")

besttime[order(time)]
besttime[batchsize==4 & nsubs==30,]

timings <- merge(sizetime, besttime)
timings$timecolour <- "#009392" # up to 1 second
timings[time > 1, timecolour := "#72aaa1"] # 1 - 10 seconds
timings[time > 10, timecolour := "#b1c7b3"] # 10 - 60 seconds
timings[time > 60, timecolour := "#e5b9ad"]  # 1 - 5 minutes
timings[time > 300, timecolour := "#d98994"] # 5 - 30 minutes
timings[time > 1800, timecolour := "#d0587e"] # > 30 minutes

## Compare actual vs estimated search space
plot(log10(timings$naive), log10(timings$space),
     main="Timing vs # calculations",
     ylab="# Calculations",
     xlab="Estimated # calculations",
     xaxt="n", yaxt="n",
     type="n", col="white")
abline(0,1,col="lightgrey")
abline(v=8.5,col="lightgrey")
abline(v=7.5,col="lightgrey")
points(log10(timings$naive), log10(timings$space),
       pch=20, col=timings$timecolour)
legend("topleft",
       legend=c("Time to exhaustively calculate\n optimal allocation",
                "< 1 sec", "1 - 10 sec", "10 - 60 sec",
                "1 - 5 min", "5 - 30 min", "> 30 min"),
       col = c("white", unique(timings[order(time), timecolour])),
       pch=20)
axis(1, at=seq(0,10),
     labels=c(expression(10^0), expression(10^1), expression(10^2),
              expression(10^3), expression(10^4), expression(10^5),
              expression(10^6), expression(10^7), expression(10^8),
              expression(10^9), expression(10^10)))
axis(2, at=seq(0,10),
     labels=c(expression(10^0), expression(10^1), expression(10^2),
              expression(10^3), expression(10^4), expression(10^5),
              expression(10^6), expression(10^7), expression(10^8),
              expression(10^9), expression(10^10)))

## compare actual vs estimated size and timing
scatter.smooth(timings$naive, timings$time/60, log="x",
               family = c("symmetric", "gaussian"),
               main="Time vs # calculations",
               xlab="# Calculations",
               ylab="Time (min) for exhaustive search",
               pch=20, col="white",
               xaxt="n", yaxt="n",
               lpars=list(col=rgb(256, 0, 0, 100, maxColorValue = 256)),
               evaluation=100)
lline <- loess.smooth(timings$space, timings$time/60,
                      evaluation = 100)
abline(h=5, col="lightgrey", lty="dashed")
abline(v=10^c(7.5), col="lightgrey", lty="dashed")
abline(v=10^8.5, col="lightgrey", lty="solid")
lines(lline, col=rgb(160, 32, 240, 100, maxColorValue = 256))
points(timings$naive, timings$time/60,
       pch=20, col="red")
points(timings$space, timings$time/60,
       pch=20, col="purple")
axis(1, at=10^seq(0,10),
     labels=c(expression(10^0), expression(10^1), expression(10^2),
              expression(10^3), expression(10^4), expression(10^5),
              expression(10^6), expression(10^7), expression(10^8),
              expression(10^9), expression(10^10)))
axis(2, at=seq(0,60,15), labels=seq(0,60,15))
legend("topleft", legend=c("# Calculations", "Estimated", "Actual"),
       col=c("white", "red", "purple"), pch=20)

## Figure CR 3
singlecolour <- rgb(0, 147, 146, max=256)
singlecolourline <- rgb(0, 147, 146, 64, max=256)
png("timevscalcexpo.png", width=600, height=600)
scatter.smooth(timings$naive, timings$time/60, log="x",
               family = c("symmetric", "gaussian"),
               main="Time vs # calculations",
               xlab="Estimated # calculations",
               ylab="Time (min) for exhaustive search",
               pch=20, col="white",
               xaxt="n", yaxt="n",
               ylim=c(0,75),
               lpars=list(col=singlecolourline),
               evaluation=100)
abline(v=10^9, col="grey", lty="dashed")
points(timings$naive, timings$time/60,
       pch=20, col=singlecolour)
axis(1, at=10^seq(0,10),
     labels=c(expression(10^0), expression(10^1), expression(10^2),
              expression(10^3), expression(10^4), expression(10^5),
              expression(10^6), expression(10^7), expression(10^8),
              expression(10^9), expression(10^10)))
axis(2, at=seq(0,75,15), labels=seq(0,75,15))
dev.off()

