library(maps)
library(mapdata)


# get SST grids
# 54-62N and 198-216E
west.x <- seq(198, 216, 2)
west.y <- seq(54, 62, 2)
lat <- rep(west.y, length(west.x))   # Vector of latitudes
lon <- rep(west.x, each = length(west.y))   # Vector of longitudes

west <- data.frame(lat=lat, lon=lon, name=paste("N", lat, "E", lon, sep=""))

#blank out Bristol Bay and far offshore areas
blank <- c("N56E198", "N58E198", "N60E198","N62E198",  "N56E200", "N58E200", "N58E202", "N54E216",
           "N56E216", "N54E214", "N56E214", "N54E212", "N54E210", "N54E208")


drop <- west$name %in% blank
keep <- grep(F, drop)
west <- west[keep,]

xx.1 <- c(197, 197, 201, 205, 208, 217, 217, 213, 213, 207, 207, 197)
yy.1 <- c(53, 55, 56, 58, 61, 61, 57, 57, 55, 55, 53, 53)

east.x <- seq(218, 232, 2)
east.y <- seq(52, 60, 2)
lat <- rep(east.y, length(east.x))   # Vector of latitudes
lon <- rep(east.x, each = length(east.y))   # Vector of longitudes

east <- data.frame(lat=lat, lon=lon, name=paste("N", lat, "E", lon, sep=""))

#blank out the far offshore areas
blank <- c("N52E218", "N54E218", "N56E218", "N52E220", "N54E220", "N52E222")

drop <- east$name %in% blank
keep <- grep(F, drop)
east <- east[keep,]

xx.2 <- c(217, 219, 219, 221, 221, 223, 223, 233, 232, 230, 228, 221, 217, 217)
yy.2 <- c(57, 57, 55, 55, 53, 53, 51, 51, 54, 56, 58, 61, 61, 57)

png("study site figure.png", 4, 6, units="in", res=300)
par(mfrow=c(3,1), tcl=0.2, cex.lab=0.8, cex.axis=1, mar=c(2,2.1,2,0.4),mgp=c(1.5,0.3,0), las=1)

land.col <- "lightyellow3"
xlim <- c(190, 233)
ylim <- c(49,62)
dummy.x <- 190:200
dummy.y <- 50:60

plot(dummy.x, dummy.y, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", yaxt="n", xaxt="n")
axis(1, at=seq(190, 230, 10), labels=c("170ºW", "160º", "150º", "140º", "130º"))
axis(2, at=seq(50, 62, 2), labels=c("50º", "52º", "54º", "56º", "58º", "60º", "62ºN"))

map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)

lines(xx, yy)
lines(xx.2, yy.2)

box()
mtext("a) Climate", adj=0, cex=0.8)

xlim <- c(195, 229)
ylim <- c(53,61)

plot(dummy.x, dummy.y, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", yaxt="n", xaxt="n")
axis(1, at=seq(195, 225, 10), labels=c("165ºW", "155º", "145º", "140º"))
axis(2, at=seq(54, 60, 2), labels=c("54º", "56º", "58º", "60ºN"))

map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)

box()
mtext("b) Biology", adj=0, cex=0.8)

xlim <- c(195, 229)
ylim <- c(53,61)

plot(dummy.x, dummy.y, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", yaxt="n", xaxt="n")
axis(1, at=seq(195, 225, 10), labels=c("165ºW", "155º", "145º", "140º"))
axis(2, at=seq(54, 60, 2), labels=c("54º", "56º", "58º", "60ºN"))

map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col=land.col)

box()
mtext("c) Lower trophic level", adj=0, cex=0.8)

dev.off()
#
