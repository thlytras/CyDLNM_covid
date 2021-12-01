# Load required packages
library(readxl)
library(Rivets)
library(FluMoDL)

source("attrdl2.R")

limwk <- 202145  #
limmo <- 202110  # Limit data to complete months

# Read in the data
dat.raw <- as.data.frame(do.call(rbind, lapply(list.files("input", "xlsx$", full=TRUE), read_excel)))[,c("DoD","Age","Sex","DoR")]
names(dat.raw) <- c("date","age","sex","DoR")
# Process the data
dat.raw$date <- as.Date(dat.raw$date)
dat.raw$obs <- 1
dat.raw$wk <- isoweek(dat.raw$date)
dat.raw <- subset(dat.raw, wk>=201001 & wk<=202145) # Limit data to avoid periods with reporting delays

# Aggregate series of daily deaths
dat <- aggregate(dat.raw[,"obs",drop=FALSE], dat.raw[,"date",drop=FALSE], sum)

# Merge in temperature data, and add further terms for DLNM
load("input/cyTemps.RData")
dat <- merge(dat, aggregate(cyTemps[,"temp",drop=FALSE], cyTemps[,"date",drop=FALSE], mean), all.x=TRUE)
dat$temp <- linterp(dat$temp) # Linearly interpolate any missing temperatures
dat$dow <- as.factor(format(dat$date,"%u"))
dat$doy <- as.integer(format(dat$date, "%j"))
dat$t <- 1:nrow(dat)

# Split into training and prediction datasets
dat.train <- subset(dat, date<"2020-1-1")
dat.predict <- subset(dat, date>=as.Date("2020-1-1")-32)

# Now fit the DLNM

# Maximum lag is 30 days
lg <- 30
# Set up the crossbasis for temperature
basis.temp <- crossbasis(dat.train$temp, lag=lg, 
    argvar=list(fun="bs", degree=2, knots=quantile(dat.train$temp, c(0.1,0.75,0.9))),
    arglag=list(fun="ns", knots=logknots(lg,3)))
# Fit the model
model <- glm(obs ~ basis.temp + dow + t + pbs(doy, knots=c(61,183,305)), dat=dat.train, family="quasipoisson")
# Get the exposure-lag-response associations for temperature, and calculate minimum mortality point
predTemp <- crosspred(basis.temp, model, at=ceiling(min(dat.train$temp)):floor(max(dat.train$temp)), bylag=0.2, cen=20, cumul=TRUE)
MMP <- as.integer(names(which(predTemp$allfit==min(predTemp$allfit)))) # Minimum mortality point
predTemp <- crosspred(basis.temp, model, at=ceiling(min(dat.train$temp)):floor(max(dat.train$temp)), bylag=0.2, cumul=TRUE,
    cen=MMP)
od <- max(1, sum(model$weights * model$residuals^2)/model$df.r)  # Overdispersion parameter phi


# Now calculate expected daily mortality for prediction dataset
datP <- (function(){
  basis.temp0 <- basis.temp  # Stash basis.temp from the global environment
  basis.temp <<- crossbasis(dat.predict$temp, lag=lg, 
      argvar=list(fun="bs", degree=2, knots=quantile(dat.train$temp, c(0.1,0.75,0.9))),
      arglag=list(fun="ns", knots=logknots(lg,3)))
  prL <- predict(model, dat.predict, se.fit=TRUE, type="link")
  prR <- predict(model, dat.predict, se.fit=TRUE, type="response")
  datP <- dat.predict[,c("date","obs","temp")]
  names(datP)[2] <- "obs"
  datP$pred <- prR$fit
  datP$pred.se <- prR$se.fit
  datP$predL <- prL$fit
  datP$predL.se <- prL$se.fit
  datP$attrTemp <- attrdl(datP$temp, basis.temp, dat.predict$obs, model, cen=MMP, type="an", tot=FALSE)
  datP$wk <- isoweek(datP$date)
  datP$mo <- as.integer(format(datP$date, "%Y%m"))
  datP <- subset(datP, wk>=202001)
  basis.temp <<- basis.temp0  # Restore `basis.temp` in the global environment
  datP
})()


datPW <- as.data.frame(do.call(rbind, by(datP, datP[,"wk"], function(x) {
  c(wk = unique(x$wk), obs=sum(x$obs), pred=round(sum(x$pred)), se=sqrt(sum(x$pred.se^2)), attrTemp=sum(x$attrTemp))
})))
# # Alternative: calculate weekly SE with bootstrapping
# t <- system.time({
# datPW <- as.data.frame(do.call(rbind, by(datP, datP[,"wk"], function(x) {
#   c(wk = unique(x$wk), obs=sum(x$obs), pred=round(sum(x$pred)), 
#     se=sd(colSums(exp(matrix(rnorm(nrow(x)*10^5, mean=x$predL, sd=x$predL.se), nrow=nrow(x))))), 
#     attrTemp=sum(x$attrTemp))
# })))
# })


datPM <- as.data.frame(do.call(rbind, by(subset(datP,date>="2020-1-1"), subset(datP,date>="2020-1-1")[,"mo"], function(x) {
  c(mo = unique(x$mo), obs=sum(x$obs), pred=round(sum(x$pred)), se=sqrt(sum(x$pred.se^2)), attrTemp=sum(x$attrTemp))
})))
# # Alternative: calculate weekly SE with bootstrapping
# t <- system.time({
# datPMb <- as.data.frame(do.call(rbind, by(subset(datP,date>="2020-1-1"), subset(datP,date>="2020-1-1")[,"mo"], function(x) {
#   c(mo = unique(x$mo), obs=sum(x$obs), pred=round(sum(x$pred)), 
#     se=sd(colSums(exp(matrix(rnorm(nrow(x)*10^6, mean=x$predL, sd=x$predL.se), nrow=nrow(x))))), 
#     attrTemp=sum(x$attrTemp))
# })))
# })

# Helper function to assist in plotting excess mortality...
plotExcessPoly <- function(o, p=NULL, col=c("pink","skyblue")) {
  if (is.null(p)) p <- rep(0, length(o))
  if (length(o)!=length(p)) stop("Observed and predicted vectors should have the same length.")
  x <- 1:length(o)
  y.lower <- pmin(o, p)
  y.upper <- pmax(o, p)
  is <- which(diff(o>p)!=0)+1
  ip <- sapply(is, function(i) {
    xp <- -(p[i-1] - o[i-1])/((p[i]-p[i-1]) - (o[i]-o[i-1]))
    yp <- p[i-1] + (p[i]-p[i-1])*xp
    xp <- xp + x[i-1]
    c(xp, yp)
  })
  x <- c(x, ip[1,])
  y.lower <- c(y.lower, ip[2,])[order(x)]
  y.upper <- c(y.upper, ip[2,])[order(x)]
  x <- sort(x)
  polygon(c(x,length(o):1), c(y.lower,rev(p)), col=col[2], border=NA)
  polygon(c(x,length(o):1), c(y.upper,rev(p)), col=col[1], border=NA)
}

# Functions to calculate Farrington prediction limits and z-scores
farrPI <- function(pred, se, Z=2, od=1) {
  pred * (1 + (2/3) * Z * ( (od + se^2/pred) /pred ) ^(1/2))^(3/2)
}
farrZ <- function(obs, pred, se, od=1) {
  ((obs/pred)^(2/3) - 1) / ((2/3) * ( (od + se^2/pred) /pred ) ^(1/2))
}

datPW$UPI2 <- farrPI(datPW$pred, datPW$se, Z=2, od=od)
datPW$UPI4 <- farrPI(datPW$pred, datPW$se, Z=4, od=od)

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}


load("input/cyCOVID.RData")

cyCOVID$wk <- isoweek(cyCOVID$date)
cyCOVIDw <- aggregate(cyCOVID[,c("new_cases","new_deaths")], cyCOVID[,"wk",drop=F], sum, na.rm=TRUE)
datPW <- merge(datPW, cyCOVIDw, all.x=TRUE)
datPW$new_deaths[is.na(datPW$new_deaths)] <- 0
datPW$new_cases[is.na(datPW$new_cases)] <- 0


cairo_pdf("output/cy_excess.pdf", width=10, height=6, pointsize=8)

(function() {
lwd <- 1.5
par(family="Fira Sans", mar=c(4,4,2,2))
layout(matrix(1:2), heights=c(1.5,1))
ofs <- floor(nrow(datPW)*0.04)
plot(datPW$pred, type="n", xaxs="i", xlim=c(1,nrow(datPW))+c(-ofs,ofs),
  ylim=range(datPW[,c("obs","pred","UPI4")]) + c(0,diff(range(datPW[,c("obs","pred","UPI4")]))*0.1), 
  xlab=NA, ylab="Weekly number of deaths", xaxt="n",
  main="Observed and expected all-cause mortality, Cyprus, 2020-2021")
mtext("Week number (ISO)", side=1, line=2.5)
abline(v=which((datPW$wk %% 100) %in% c(1,20,40)), col="grey", lty="dotted")
axis(1, at=which((datPW$wk %% 100) %in% c(1,20,40)), labels=datPW$wk[which((datPW$wk %% 100) %in% c(1,20,40))])
plotExcessPoly(datPW$obs, datPW$pred, col=addalpha(c("red","skyblue"), 0.3))
#polygon(c(1:nrow(datPW),nrow(datPW):1), with(datPW, c(obs-new_deaths, rev(obs))), density=30)
polygon(c(1:nrow(datPW),nrow(datPW):1), with(datPW, c(obs-new_deaths, rev(obs))), border=NA, col=addalpha("green", 0.4))
points(datPW$UPI2, type="l", col="orange", lwd=lwd, lty="dashed")
points(datPW$UPI4, type="l", col="red", lwd=lwd, lty="dashed")
points(datPW$pred, type="l", col="grey", lwd=lwd, lty="dashed")
points(datPW$obs, type="l", lwd=lwd)
legend("topright", c("Observed", "Expected", "Expected + 2 SD", "Expected + 4 SD", 
    "Excess mortality", "Mortality deficit", "COVID-19 deaths"),
  col=c("black", "grey", "orange","red",addalpha(c("red","skyblue"), 0.3), addalpha("green", 0.4)), 
  lwd=lwd, lty=c("solid",rep("dashed",3), NA, NA, NA),
  pch=c(rep(NA,4),15,15,15), pt.cex=2, ncol=2, bg="white", seg.len=3.5, adj=0.05)
with(datP, {
  plot(date, temp, type="l", xaxs="i", xlim=range(date)+c(-1,1)*(ofs*7-3), xaxt="n", lwd=lwd*2/3,
    main="Mean daily temperature, Cyprus, 2020-2021", xlab=NA, ylab="Temperature (celsius)")
  daxs <- sapply(datPW$wk[which((datPW$wk %% 100) %in% c(1,20,40))], function(w) date[match(w, wk)+3])
  abline(v=daxs, col="grey", lty="dotted")
  axis(1, at=daxs, label=as.Date(daxs, origin="1970-1-1"))
})
})()

dev.off()


png("output/Cy_temp_DLNM.png", width=1000, height=700, res=140)
plot(predTemp, "overall", lwd=2, ylab="Mortality Rate Ratio (RR)", xlab="Mean daily temperature (celsius)",
  main="Temperature-mortality association, across all lags")
mtext("(maximum lag = 30 days)", side=3, line=0)
dev.off()
