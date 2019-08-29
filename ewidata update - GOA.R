library(tidyverse)
library(ncdf4)
library(chron)
library(fields)
library(maps)
library(mapdata)
library(zoo)

dat <- read.csv("ewidata.csv", row.names = 1)

levels(dat$system)

# restrict to GOA 

goa.dat <- dat %>%
  filter(system %in% c("EGOA", "WGOA"))

ggplot(goa.dat, aes(year, value)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  facet_wrap(~code, scales="free_y") 

# extend SST data
# load ERSSTv5 data for the N. Pacific
# 20º-70ºN, 120º-250ºE, 1854-present

# identify latest year and month needed
year <- 2019
month <- "07"

URL <- paste("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1854-01-01):1:(", year, "-", month, "-01T00:00:00Z)][(0.0):1:(0.0)][(20):1:(70)][(120):1:(250)]", sep="")

download.file(URL, "North.Pacific.ersst")
nc <- nc_open("North.Pacific.ersst")

# extract dates
ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
sst.d <- dates(h, origin = c(1,1,1970))

sst.x <- ncvar_get(nc, "longitude")
sst.y <- ncvar_get(nc, "latitude")

# save months and years for use later on
m <- months(sst.d)
yrs <- years(sst.d)

# and set a couple functions for standardizing below
f1 <- function(x) tapply(x, m, mean)
f2 <- function(x) tapply(x, m, sd)

SST <- ncvar_get(nc,  "sst")
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
sst.lat <- rep(sst.y, length(sst.x))   
sst.lon <- rep(sst.x, each = length(sst.y))   
dimnames(SST) <- list(as.character(sst.d), paste("N", sst.lat, "E", sst.lon, sep=""))

##############################################################################
# Western GOA
# 54-62N and 198-216E
keep <- sst.lat %in% 54:62 & sst.lon %in% 198:216
SST <- SST[,keep]
sst.lat <- sst.lat[keep]
sst.lon <- sst.lon[keep]
sst.y <- sst.y[sst.y %in% 54:62]
sst.x <- sst.x[sst.x %in% 198:216]

#blank out Bristol Bay and far offshore areas
blank <- c("N56E198", "N58E198", "N60E198","N62E198",  "N56E200", "N58E200", "N58E202", "N54E216",
           "N56E216", "N54E214", "N56E214", "N54E212", "N54E210", "N54E208")
SST[,blank] <- NA

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(40,66))
contour(sst.x, sst.y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# looks good!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(sst.d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)
add <- length(sst.d)-12*floor(length(sst.d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(sst.d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# check with some plots!
par(las=1)
wgoa.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1854,1), frequency=12)
plot(wgoa.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1854,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("WGOA SST & 13-mo rolling mean")

# now calculate seasonal means
# define seasons for Alaska
win <- c("Nov", "Dec", "Jan", "Feb","Mar")
spr <- c("Apr", "May", "Jun")

# and define winter years
win.yrs <- as.numeric(as.character(yrs)) # change winter years to a numeric object
# now we need to assign Nov and Dec data to the year corresponding to Jan
win.yrs[m %in% c("Nov", "Dec")] <- win.yrs[m %in% c("Nov", "Dec")]+1

# calculate winter mean SST
wgoa.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
wgoa.mean <- wgoa.mean[use] # select winter means only
win.yrs <- win.yrs[use] # restrict win.yrs to winter months only
wgoa.win <- tapply(wgoa.mean, win.yrs, mean)

# plot to check
plot(names(wgoa.win), wgoa.win, type="b")

# calculate spring mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% spr # logical vector that's true for spring months
SST.mean <- SST.mean[use] # select spring means only
spr.yrs <- yrs[use] # restrict years vector to spring months
wgoa.spr <- tapply(SST.mean, spr.yrs, mean)

plot(names(wgoa.spr), wgoa.spr, type="b")

##########################################################################################
# now the Eastern GOA
SST <- ncvar_get(nc,  "sst")
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  
sst.x <- ncvar_get(nc, "longitude")
sst.y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc,  "sst")
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
sst.lat <- rep(sst.y, length(sst.x))   
sst.lon <- rep(sst.x, each = length(sst.y))   
dimnames(SST) <- list(as.character(sst.d), paste("N", sst.lat, "E", sst.lon, sep=""))

# 52-60N and 218-232E
keep <- sst.lat %in% 52:60 & sst.lon %in% 218:232
SST <- SST[,keep]
sst.lat <- sst.lat[keep]
sst.lon <- sst.lon[keep]
sst.y <- sst.y[sst.y %in% 52:60]
sst.x <- sst.x[sst.x %in% 218:232]

#blank out the far offshore areas
blank <- c("N52E218", "N54E218", "N56E218", "N52E220", "N54E220", "N52E222")
SST[,blank] <- NA

# plot mean temperature pattern to check
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(sst.y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(sst.x,sst.y,z, col=tim.colors(64), xlim=c(160,240), ylim=c(40,66))
contour(sst.x, sst.y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)
# looks good!

# now remove seasonal signal and scale 
mu <- apply(SST, 2, f1)	# Compute monthly means for each time series (location)
mu <- mu[rep(1:12, floor(length(sst.d)/12)),]  # Replicate means matrix for each year at each location 

#now need to account for trailing months (i.e., fractions of a year that are available with updated data)
add <- length(sst.d)-12*floor(length(sst.d)/12)

# add in additional required months
mu <- rbind(mu, mu[1:add,])
# check
identical(nrow(mu), nrow(SST)) #true!

std <- apply(SST, 2, f2)	# Compute monthly sd for each time series (location)
# and stack as for mu
std <- std[rep(1:12, floor(length(sst.d)/12)),] 
# add in additional required months
std <- rbind(std, std[1:add,])

# now calculate anomalies...
SST.anom <- (SST - mu)/std

# check with some plots!
par(las=1)
egoa.sst <- ts(rowMeans(SST.anom, na.rm=T), start=c(1854,1), frequency=12)
plot(egoa.sst, type="l", col="grey", lwd=0.8, xlab="", ylab="SST anomaly")
sm <- ts(rollmean(rowMeans(SST.anom, na.rm=T), 13, fill=NA), start=c(1854,1), frequency=12)
lines(sm, col="red", lwd=1.5)
abline(h=0)
mtext("EGOA SST & 13-mo rolling mean")

# calculate winter mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% win # logical vector that's true for winter months
SST.mean <- SST.mean[use] # select winter means only

egoa.win <- tapply(SST.mean, win.yrs, mean)

# calculate spring mean SST
SST.mean <- rowMeans(SST.anom, na.rm=T) # means in every month
use <- m %in% spr # logical vector that's true for spring months
SST.mean <- SST.mean[use] # select spring means only
egoa.spr <- tapply(SST.mean, spr.yrs, mean)

# collect new sst values
head(dat)
new.sst <- data.frame(year=rep(1950:2019, 4), 
                      code=rep(c("AKCLIM_egoa.spr.sst", "AKCLIM_egoa.win.sst", "AKCLIM_wgoa.spr.sst", "AKCLIM_wgoa.win.sst"), each=length(1950:2019)),
                      value=c(egoa.spr[names(egoa.spr) %in% 1950:2019], egoa.win[names(egoa.win) %in% 1950:2019], wgoa.spr[names(wgoa.spr) %in% 1950:2019], wgoa.win[names(wgoa.win) %in% 1950:2019]),
                      system=rep(c("EGOA", "WGOA"), each=length(1950:2019)), subtype=NA)

drop <- c(grep("egoa", dat$code), grep("wgoa", dat$code))
dat <- dat[-drop,]

dat <- rbind(dat, new.sst)

goa.dat <- dat %>%
  filter(system %in% c("EGOA", "WGOA"))

ggplot(goa.dat, aes(year, value)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  facet_wrap(~code, scales="free_y") 

# load covariates from AK salmon analysis
covar.dat <- read.csv("salmon.covariates.csv")
head(covar.dat)
