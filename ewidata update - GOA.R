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

# get Papa updates
papa <- read.csv("xtra.papa.csv")
papa

covar.dat$Papa[54:56] <- papa$papa[3:5]

# now load updated upwelling station data
uw <- read.csv("upwelling.csv")

levels(uw$POSITION)

uw <- uw %>%
  filter(POSITION %in% c("54N.134W", "57N.137W", "60N.146W", "60N.149W")) %>%
  select(1,2,7:9) %>%
  gather(key="MONTH", value, -POSITION, -YEAR) %>%
  arrange(POSITION, YEAR, MONTH)

sum.uw <- uw %>%
  group_by(POSITION, YEAR) %>%
  summarise(mean=mean(value)) %>%
  filter(YEAR >= 1950)

unique(goa.dat$code)

# I'm gonna make a new df for GOA climate data!

names(sum.uw) <- c("code", "year", "value")

sum.uw <- tapply(sum.uw$value, list(sum.uw$year, sum.uw$code), identity)
drop <- is.na(colMeans(sum.uw))
sum.uw <- sum.uw[,!drop]

# load stress / ssh / slp gradient
goa.clim <- read.csv("long-term goa ssh stress slp gradient.csv")

goa.clim <- cbind(goa.clim, sum.uw)

keep <- grep("sst", goa.dat$code)
temp <- goa.dat[keep,]
temp <- tapply(temp$value, list(temp$year, temp$code), identity)
drop <- is.na(colMeans(temp))
temp <- temp[,!drop]

goa.clim <- cbind(goa.clim, temp)
names(goa.clim)[10:13] <- c("egoa.spr.sst", "egoa.win.sst", "wgoa.spr.sst", "wgoa.win.sst")

temp <- goa.dat %>%
  filter(code=="AKCLIM_GOA_PAPA")

goa.clim$papa.index <- c(temp$value, papa$papa[3:5])

plot.clim <- goa.clim %>%
  gather(key="code", value="value", -year)

ggplot(plot.clim, aes(year, value)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  facet_wrap(~code, scale="free_y")

write.csv(goa.clim, "updated goa climate data.csv", row.names=F)

#######################################
# now run the dang climate DFA model! #
#######################################

library(Rcpp)
library(ewidata)
library(knitr)
library(reshape2)
library(rstan)
library(tibble)
library(dplyr)
library(gtools)
# devtools::install_github("fate-ewi/bayesdfa", force=T)

library(bayesdfa)

mcmc_iter = 4000
max_trends = 3
mcmc_chains = 3
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# doing this in a clunky way -
# find dfa trends with 2 or 3 candidate trends
# and original data (as in previous iteration of paper)
# or expanded data (adding SSH and wind stress)
goa.clim <- read.csv("updated goa climate data.csv")

sub_data <- goa.clim %>%
  gather(key="code", value="value", -year)

# set the name for the model output!
name <- "GOA_clim_3_trends_expanded_data"

# reshape data
melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

# do the trend search
set.seed(99)
dfa_summary = find_dfa_trends(
  y = Y,
  kmax = min(max_trends, nrow(Y)),
  iter = mcmc_iter,
  compare_normal = FALSE,
  variance = c("unequal", "equal"),
  chains = mcmc_chains
)
saveRDS(dfa_summary, file = paste0(name, ".rds"))

# Make default plots (currently work in progress)
pdf(paste0(name, "_plots.pdf"))
rotated = rotate_trends(dfa_summary$best_model)
# trends
print(plot_trends(rotated, years = as.numeric(colnames(Y))))
# loadings
print(plot_loadings(rotated, names = names))

if(ncol(rotated$Z_rot_mean)==2) {
  plot(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], col="white",
       xlab="Loading 1", ylab = "Loading 2")
  text(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], names, cex=0.3)
  lines(c(-10,10),c(0,0))
  lines(c(0,0), c(-10,10))
}
# predicted values with data
print(plot_fitted(dfa_summary$best_model,names=names) +
        theme(strip.text.x = element_text(size = 6)))

# table of AIC and std errors
summary_table<-dfa_summary$summary
capture.output(summary_table, file = paste0(name, "_summary.txt"))

dev.off()

# now two trends with expanded data
max_trends = 2
name <- "GOA_clim_2_trends_expanded_data"

# reshape data
melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

# do the trend search
set.seed(99)
dfa_summary = find_dfa_trends(
  y = Y,
  kmax = min(max_trends, nrow(Y)),
  iter = mcmc_iter,
  compare_normal = FALSE,
  variance = c("unequal", "equal"),
  chains = mcmc_chains
)
saveRDS(dfa_summary, file = paste0(name, ".rds"))

# Make default plots (currently work in progress)
pdf(paste0(name, "_plots.pdf"))
rotated = rotate_trends(dfa_summary$best_model)
# trends
print(plot_trends(rotated, years = as.numeric(colnames(Y))))
# loadings
print(plot_loadings(rotated, names = names))

if(ncol(rotated$Z_rot_mean)==2) {
  plot(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], col="white",
       xlab="Loading 1", ylab = "Loading 2")
  text(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], names, cex=0.3)
  lines(c(-10,10),c(0,0))
  lines(c(0,0), c(-10,10))
}
# predicted values with data
print(plot_fitted(dfa_summary$best_model,names=names) +
        theme(strip.text.x = element_text(size = 6)))

# table of AIC and std errors
summary_table<-dfa_summary$summary
capture.output(summary_table, file = paste0(name, "_summary.txt"))

dev.off()

## ############################
# now the original data set
sub_data <- goa.clim %>%
  select(-ssh, -wind.stress) %>%
  gather(key="code", value="value", -year)

# ONE! trend model
# with original data
max_trends = 1
name <- "GOA_clim_1_trend_original_data"

# reshape data
melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

# do the trend search
set.seed(99)
dfa_summary = find_dfa_trends(
  y = Y,
  kmax = min(max_trends, nrow(Y)),
  iter = mcmc_iter,
  compare_normal = FALSE,
  variance = c("unequal", "equal"),
  chains = mcmc_chains
)
saveRDS(dfa_summary, file = paste0(name, ".rds"))

# Make default plots (currently work in progress)
pdf(paste0(name, "_plots.pdf"))
rotated = rotate_trends(dfa_summary$best_model)
# trends
print(plot_trends(rotated, years = as.numeric(colnames(Y))))
# loadings
print(plot_loadings(rotated, names = names))

if(ncol(rotated$Z_rot_mean)==2) {
  plot(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], col="white",
       xlab="Loading 1", ylab = "Loading 2")
  text(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], names, cex=0.3)
  lines(c(-10,10),c(0,0))
  lines(c(0,0), c(-10,10))
}
# predicted values with data
print(plot_fitted(dfa_summary$best_model,names=names) +
        theme(strip.text.x = element_text(size = 6)))

# table of AIC and std errors
summary_table<-dfa_summary$summary
capture.output(summary_table, file = paste0(name, "_summary.txt"))

dev.off()


######################
# two trends with original data
max_trends = 2
name <- "GOA_clim_2_trends_original_data"

# reshape data
melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

# do the trend search
set.seed(99)
dfa_summary = find_dfa_trends(
  y = Y,
  kmax = min(max_trends, nrow(Y)),
  iter = mcmc_iter,
  compare_normal = FALSE,
  variance = c("unequal", "equal"),
  chains = mcmc_chains
)
saveRDS(dfa_summary, file = paste0(name, ".rds"))

# Make default plots (currently work in progress)
pdf(paste0(name, "_plots.pdf"))
rotated = rotate_trends(dfa_summary$best_model)
# trends
print(plot_trends(rotated, years = as.numeric(colnames(Y))))
# loadings
print(plot_loadings(rotated, names = names))

if(ncol(rotated$Z_rot_mean)==2) {
  plot(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], col="white",
       xlab="Loading 1", ylab = "Loading 2")
  text(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], names, cex=0.3)
  lines(c(-10,10),c(0,0))
  lines(c(0,0), c(-10,10))
}
# predicted values with data
print(plot_fitted(dfa_summary$best_model,names=names) +
        theme(strip.text.x = element_text(size = 6)))

# table of AIC and std errors
summary_table<-dfa_summary$summary
capture.output(summary_table, file = paste0(name, "_summary.txt"))

dev.off()

######################
# one-trend models fit separately to first 15 and last 15 years of the TS
sub_data <- goa.clim %>%
  select(-ssh, -wind.stress) %>%
  filter(year %in% 1972:1986) %>%
  gather(key="code", value="value", -year)

# ONE! trend model
# with original data
max_trends = 1
# name <- "GOA_clim_1_trend_original_data_1950_1964"
name <- "GOA_clim_1_trend_original_data_1972_1986"
# reshape data
melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

# do the trend search
set.seed(99)
dfa_summary = find_dfa_trends(
  y = Y,
  kmax = min(max_trends, nrow(Y)),
  iter = mcmc_iter,
  compare_normal = FALSE,
  variance = c("unequal", "equal"),
  chains = mcmc_chains
)
saveRDS(dfa_summary, file = paste0(name, ".rds"))

# Make default plots (currently work in progress)
pdf(paste0(name, "_plots.pdf"))
rotated = rotate_trends(dfa_summary$best_model)
# trends
print(plot_trends(rotated, years = as.numeric(colnames(Y))))
# loadings
print(plot_loadings(rotated, names = names))

if(ncol(rotated$Z_rot_mean)==2) {
  plot(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], col="white",
       xlab="Loading 1", ylab = "Loading 2")
  text(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], names, cex=0.3)
  lines(c(-10,10),c(0,0))
  lines(c(0,0), c(-10,10))
}
# predicted values with data
print(plot_fitted(dfa_summary$best_model,names=names) +
        theme(strip.text.x = element_text(size = 6)))

# table of AIC and std errors
summary_table<-dfa_summary$summary
capture.output(summary_table, file = paste0(name, "_summary.txt"))

dev.off()


# last 15
sub_data <- goa.clim %>%
  select(-ssh, -wind.stress) %>%
  filter(year >= 2005) %>%
  gather(key="code", value="value", -year)


# ONE! trend model
# with original data
max_trends = 1
name <- "GOA_clim_1_trend_original_data_2005_2019"

# reshape data
melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

# do the trend search
set.seed(99)
dfa_summary = find_dfa_trends(
  y = Y,
  kmax = min(max_trends, nrow(Y)),
  iter = mcmc_iter,
  compare_normal = FALSE,
  variance = c("unequal", "equal"),
  chains = mcmc_chains
)
saveRDS(dfa_summary, file = paste0(name, ".rds"))

# Make default plots (currently work in progress)
pdf(paste0(name, "_plots.pdf"))
rotated = rotate_trends(dfa_summary$best_model)
# trends
print(plot_trends(rotated, years = as.numeric(colnames(Y))))
# loadings
print(plot_loadings(rotated, names = names))

if(ncol(rotated$Z_rot_mean)==2) {
  plot(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], col="white",
       xlab="Loading 1", ylab = "Loading 2")
  text(rotated$Z_rot_mean[,1], rotated$Z_rot_mean[,2], names, cex=0.3)
  lines(c(-10,10),c(0,0))
  lines(c(0,0), c(-10,10))
}
# predicted values with data
print(plot_fitted(dfa_summary$best_model,names=names) +
        theme(strip.text.x = element_text(size = 6)))

# table of AIC and std errors
summary_table<-dfa_summary$summary
capture.output(summary_table, file = paste0(name, "_summary.txt"))

dev.off()

########################################################
# make an era-specific plot of loadings!

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# join posteriors for loadings in the two eras
# load the first 15 year results
# load the 1-trend object
GOA.clim.era1 <- readRDS("GOA_clim_1_trend_original_data_1972_1986.rds") # read in GOA climate model

# recover the years and names for plots!
goa.clim <- read.csv("updated goa climate data.csv")

sub_data <- goa.clim %>%
  select(-ssh, -wind.stress) %>%
  gather(key="code", value="value", -year)

melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

rotated.1 = rotate_trends(GOA.clim.era1$best_model)

# set new names
new.names <- c("East spring SST", "East winter SST", "SLP gradient", "Papa advection",
               "GAK1 salinity", "West spring SST", "West winter SST", "Downwell 54 134", "Downwell 57 137", "Downwell 60 146", "Downwell 60 149")

# back to the bespoke code!
loadings.1 <- as.data.frame(rotated.1$Z_rot[,,1])
names(loadings.1)  <- new.names

# now....era2
GOA.clim.era2 <- readRDS("GOA_clim_1_trend_original_data_2005_2019.rds") # read in GOA climate model
rotated.2 = rotate_trends(GOA.clim.era2$best_model)
loadings.2 <- as.data.frame(rotated.2$Z_rot[,,1])
names(loadings.2)  <- new.names

loadings.1$era <- "1972-1986"
loadings.2$era <- "2005-2019"


# plot
plot.load <- rbind(loadings.1, loadings.2) %>%
  gather(key, value, -era)


rank <- tapply(plot.load$value, plot.load$key, mean, na.rm=T)
plot.load$rank <- rank[match(plot.load$key, names(rank))]
plot.load$key <- reorder(plot.load$key, -plot.load$rank)

# set pallette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot.a <- ggplot(plot.load, aes(key, value, fill=era)) +
  geom_violin(color=NA) +
  scale_fill_manual(values=cb[c(2,4)]) +
  theme_bw() +
  ylim(-2, 2) + geom_hline(yintercept = 0, lty=2) +
  theme(legend.title = element_blank(), axis.title.y = element_blank(), legend.position = 'top',
        title = element_text(size=9)) +
  coord_flip() + ylab("Loading") + ggtitle("a) Era-specific loadings")

plot.trend <- data.frame(year=1972:1986, trend=as.vector(rotated.1$trends_mean),
                         lo=as.vector(rotated.1$trends_lower),
                         hi=as.vector(rotated.1$trends_upper))

plot.b <- ggplot(plot.trend, aes(year, trend)) +
  theme_bw() +
  geom_line(color=cb[2]) +
  geom_ribbon(aes(ymin=lo, ymax=hi), fill=cb[2], alpha=0.4) +
  geom_hline(yintercept = 0) + ylab("Trend value") + xlab("Year") +
  ggtitle("b) Trend 1972-1986") +
  theme(title = element_text(size=8))

plot.trend <- data.frame(year=2005:2019, trend=as.vector(rotated.2$trends_mean),
                         lo=as.vector(rotated.2$trends_lower),
                         hi=as.vector(rotated.2$trends_upper))

plot.c <- ggplot(plot.trend, aes(year, trend)) +
  theme_bw() +
  geom_line(color=cb[4]) +
  geom_ribbon(aes(ymin=lo, ymax=hi), fill=cb[4], alpha=0.4) +
  geom_hline(yintercept = 0) + ylab("Trend value") + xlab("Year") +
  ggtitle("c) Trend 2005-2019") +
  theme(title = element_text(size=8))

plot.null <- ggplot() + theme_void()

png("era-specific climate dfa.png", 6,6, units="in", res=300)

ggpubr::ggarrange(plot.a,
                  ggpubr::ggarrange(plot.null, plot.b, plot.c, ncol=1, nrow=3, heights = c(0.2,1,1)),
                  ncol=2, widths=c(1,0.8))

dev.off()

