###################################################
### chunk number 1: ReutersTickData
###################################################
gbp <- read.table(file="data/gbp_usd_tick_29_10_2007.csv", 
 sep=",", col.names=c("time","last"))
gbp$time <- as.POSIXct(strptime(gbp$time, "%d/%m/%Y%H:%M:%S"))
plot(gbp$time[2:length(gbp$time)], diff(gbp$time), type='h', 
 main="GBP/USD Inter-tick Arrival Times", 
 ylab="Intertick Time (sec)", xlab="Time")


###################################################
### chunk number 2: EBSHighFreqData
###################################################
require(zoo)
# Create a date/time conversion function
# that uses fractional seconds
convertDate <- function(x) {
	# 2007-11-23 14:48:43.140
	as.POSIXct(strptime(x, "%Y-%m-%d %H:%M:%OS"))
}

# Convenience function (from zoo quickref)
is.weekend <- function(x) {
 x <- as.POSIXlt(x)
 x$wday > 5 | x$wday < 1
}

# Load the data file 
# (Approx 1 month EBS EUR/USD tick data)
eurusd <- read.zoo(file="data/eurusd_EBS_11_07.dat", 
	sep="|", FUN=as.POSIXct,
	col.names=c("time","bid","ask"))

# EBS ticks can be one-sided (bid or ask = 0)
# Replace zeros with NA and use na.approx() to 
# interpolate ticks using previous values
coredata(eurusd)[,"bid"][coredata(eurusd)[,"bid"]==0] = NA
coredata(eurusd)[,"ask"][coredata(eurusd)[,"ask"]==0] = NA
eurusd <- na.locf(eurusd)

# Now strip out the weekend data
eurusd <- eurusd[!is.weekend(time(eurusd))]

# This is a bad way to do this (TODO find out how to
# aggregate and extract series based on time)
ts <- eurusd[as.Date(time(eurusd))=="2007-11-12"]
plot(time(ts)[2:length(time(ts))], diff(time(ts)), type='h', 
 main="EBS EUR/USD Tick Interarrival Times\n12-11-2007", 
 xlab="Time", ylab="Interrarival Time (sec)")



###################################################
### chunk number 3: TimeSeriesModels
###################################################
par(mfrow=c(2,2))
# Create MA(1) Series
# `$ y_t = \epsilon_t + \beta \epsilon_{t-1}$`
ma1 <- arima.sim(list(order=c(0,0,1), ma=.8), n=100)
plot(ma1, type="o", main=expression(paste("MA(1) ",beta,"=0.8")),
 xlab="", ylab="")
maAcf <- acf(ma1, plot=FALSE)
plot(maAcf, type='h', main="MA(1) ACF", lwd=2)
points(maAcf$lag, maAcf$acf, pch=20, cex=2)
# Create AR(1) Series
# `$y_t = \alpha y_{t-1} + \epsilon_t$`
ar1 <- arima.sim(list(order=c(1,0,0), ar=.8), n=100)
plot(ar1, type="o", main=expression(paste("AR(1) ",alpha,"=0.8")),
 xlab="", ylab="")
arAcf <- acf(ar1, plot=FALSE)
plot(arAcf, type='h', main="AR(1) ACF", lwd=2)
points(arAcf$lag, arAcf$acf, pch=20, cex=2)


###################################################
### chunk number 4: QuasiRandomNumbers
###################################################
library(gsl)
q <- qrng_alloc(type="sobol", 2)
rs <- qrng_get(q,1000)
par(mfrow=c(3,1))
plot(rnorm(1000), rnorm(1000), pch=20, main="~N(0,1)",
 ylab="", xlab="")
plot(rs, pch=20, main="Sobol",
 ylab="", xlab="")
plot(rcauchy(1000), rcauchy(1000),pch=20,
main="~C(0,1)", ylab="",xlab="")


###################################################
### chunk number 5: 3dVariates
###################################################
library(gsl)
library(lattice)
q <- qrng_alloc(type="sobol", 3)
npoints <- 200
rs <- qrng_get(q,npoints)
ltheme <- canonical.theme(color = FALSE)
ltheme$strip.background$col <- "transparent"
lattice.options(default.theme = ltheme)
trellis.par.set(layout.heights =
list(top.padding = -20,
main.key.padding = 1,
key.axis.padding = 0,
axis.xlab.padding = 0,
xlab.key.padding = 0,
key.sub.padding = 0,
bottom.padding = -20))

# Plot the normal variates in a 3-dim cube
p1 <- cloud(rnorm(npoints) ~ rnorm(npoints) + rnorm(npoints), xlab="x", ylab="y",
	zlab="z", pch=20, main="~N(0,1)")
p2 <- cloud(rs[,1] ~ rs[,2] + rs[,3], xlab="x", ylab="y",
	zlab="z", pch=20, main="Sobol")
print(p1, split=c(1,1,2,1), more=TRUE)
print(p2, split=c(2,1,2,1))
detach("package:gsl")   # gsl overloads stats::sd function
unloadNamespace("gsl")  # so unload it for now


###################################################
### chunk number 6: BinomialLattice
###################################################
# Generates drawing instructions for dot
# from an array of lattice values
dotlattice <- function(S, labels=FALSE) {
shape <- ifelse(labels == TRUE, "plaintext", "point")

cat("node[shape=",shape,", samehead, sametail];","\n", sep="")
cat("rankdir=LR;","\n")

cat("edge[arrowhead=none];","\n")

# Create a dot node for each element in the lattice
for (i in 1:length(S)) {
    cat("node", i, "[label=\"", S[i], "\"];", "\n", sep="")
}

# The number of levels in a binomial lattice of length N
# is `$\frac{\sqrt{8N+1}-1}{2}$`
L <- ((sqrt(8*length(S)+1)-1)/2 - 1)

k<-1
for (i in 1:L) {
tabs <- rep("\t",i-1)
j <- i
while(j>0) {
 cat("node",k,"->","node",(k+i),";\n",sep="")
 cat("node",k,"->","node",(k+i+1),";\n",sep="")
 k <- k + 1
 j <- j - 1
}
}
}

# Generate a binomial lattice
# for a given up, down, start value and number of steps
genlattice <- function(X0=100, u=1.1, d=.75, N=5) {
X <- c()
X[1] <- X0
count <- 2

for (i in 1:N) {
for (j in 0:i) {
 X[count] <- X0 * u^j * d^(i-j)
 count <- count + 1
}
}
return(X)
}


###################################################
### chunk number 7: LatticeGraph
###################################################
dotlattice(genlattice(N=8), labels=FALSE)


###################################################
### chunk number 8: BernoulliTrials
###################################################
n <- 10
y <- c(0:10)
op <- par(mfrow=c(4,3))
# Compute `${n \choose y}w^y(1-w)^{n-y}, 0 \leq w \leq 1 $`
for (w in seq(0,1,.1)) {
       barplot(choose(n,y)*(w)^y*(1-w)^(n-y),
        ylab="", xlab="y",main=paste("w =",w))
}
par(op)


###################################################
### chunk number 9: MaximumLikelihood
###################################################
# Calculate maximum likelihood value of w
w <- seq(0,1,.001)
y <- 7
# f is the function `${n \choose y}w^y(1-y)^{n-y}$`
f <- function(w) {choose(n,y)*w^y*(1-w)^(n-y)}
maxlk <- optimize(f, interval=c(0,1), maximum=TRUE)$maximum
plot(w,f(w), main="Likelihood Function",ylab="f(w|n,y)",xlab="w")


###################################################
### chunk number 10: ReutersFXDensity
###################################################
# Function that reads Reuters CSV tick data and converts Reuters dates
# Assumes format is Date,Tick
readRTD <- function(filename) {
tickData <- read.csv(file=filename, header=TRUE, col.names=c("Date","Tick"))
tickData$Date <- as.POSIXct(strptime(tickData$Date, format="%d/%m/%Y %H:%M:%S"))
tickData
}

# Boilerplate function for Reuters FX tick data transformation and density plot
plot.reutersFXDensity <- function() {
filenames <- c("data/eur_usd_tick_26_10_2007.csv",
	"data/eur_usd_1min_26_10_2007.csv",
	"data/eur_usd_5min_26_10_2007.csv",
	"data/eur_usd_hourly_26_10_2007.csv",
	"data/eur_usd_daily_26_10_2007.csv")
labels <- c("Tick", "1 Minute", "5 Minutes", "Hourly", "Daily")

par(mfrow=c(length(filenames), 2),mar=c(0,0,2,0), cex.main=2)
tickData <- c()
i <- 1
for (filename in filenames) {
 tickData[[i]] <- readRTD(filename)
 # Transform: `$Y = \nabla\log(X_i)$`
 logtick <- diff(log(tickData[[i]]$Tick))
 # Normalize: `$\frac{(Y-\mu_Y)}{\sigma_Y}$`
 logtick <- (logtick-mean(logtick))/sd(logtick)
 # Theoretical density range: `$\left[\lfloor\mathrm{min}(Y)\rfloor,\lceil\mathrm{max}(Y)\rceil\right]$`
 x <- seq(floor(min(logtick)), ceiling(max(logtick)), .01)
 plot(density(logtick), xlab="", ylab="", axes=FALSE, main=labels[i])
 lines(x,dnorm(x), lty=2)
 #legend("topleft", legend=c("Empirical","Theoretical"), lty=c(1,2))
 plot(density(logtick), log="y", xlab="", ylab="", axes=FALSE, main="Log Scale")
 lines(x,dnorm(x), lty=2)
 i <- i + 1
}
par(op)
}

plot.reutersFXDensity()


###################################################
### chunk number 11: SimpleLinearRegression
###################################################
trees.lm <- lm(Girth ~ Height, data=trees)
plot(Girth ~ Height, data=trees, 
	main=expression(paste("Linear Regression Example",
	beta,"=",
	trees.lm$coefficients[2])))
abline(coef=trees.lm$coef)


###################################################
### chunk number 12: gbpUsdVariogram
###################################################
# We have ~ 1 working day of tick values versus 25 years of daily values
timescales <- c("tick","1min","5min","10min","30min","60min",
 "daily","weekly","monthly")

filenames <- c()
# Create the vector of filenames for this exercise
for (i in seq(1,along.with=timescales)) {
 filenames[i] <- paste("data/gbp_usd_",timescales[i],"_29_10_2007.csv",sep="")
}

# Variogram Function
plot.tickDatavariogram <- function(filenames, labels) {
 n <- length(filenames)
 vars <- numeric(n)
 # Read and log-transform the data, and calculate sample variance
 for (i in 1:n) {
 data <- readRTD(filenames[i])
 vars[i] <- var(diff(log(data$Tick)))
 }

 # Now create a custom plot
 par(las=2)
 plot(vars, axes=FALSE, ylab="", xlab="", main="GBP/USD Variogram",
 cex.lab=0.8, type='o', pch=20, log="y", sub="Sample Variance (log scale)")
 axis(1,at=c(1:length(labels)), labels=labels)
 axis(2)
}

plot.tickDatavariogram(filenames, timescales)


###################################################
### chunk number 13: activationFunctions
###################################################
x <- seq(-5,5,.01)
a <- pnorm(x)           # Cumulative Gaussian
b <- tanh(x)            # tansig
c <- 1/(1+exp(-x))      # Logsimoid
plot(a, ylim=range(a,b,c), ylab="Output", xlab="Input",
 main="Activation Functions", type="l")
lines(b, lty=2)
lines(c, lty=3)
legend("bottomright", lty=c(1,2,3), legend=c("Gaussian CDF", "Tansig", "Logsimoid"))



