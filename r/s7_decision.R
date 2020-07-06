#
# IKI Bangladesh (MIOASI): S7 Decision Theory
#
# Define loss function and use decsision theory to find
# optimal action to take given the posterior prediction from S5
#
# Author: HS with TE
# Created: Jan 2020

library(RNetCDF)
library(reshape2)
library(fields)
library(evd)
library(mgcv)
library(plyr)
library(rogme)
library(abind)
library(doParallel)

STORMNAMES <- c("AILA", "AKASH", "BOB01", "BOB07", "BULBUL", "FANI", "MORA", "RASHMI", "ROANU", "SIDR", "TC01B", "VIYARU")
VAR <- "fg.T1Hmax"
RES <- "4p4"
INDIR <- ""

# Load preds
load("grand_fp_preds.Rdata")


# Load base netcdf
filename <- paste("fp.", VAR, ".", STORMNAMES[1], ".", RES, "km.nc", sep = "")
stormsnc <- open.nc(paste(INDIR, filename, sep = "/"))

# Extract netcdf variables
lat <- var.get.nc(stormsnc, "latitude")
lon <- var.get.nc(stormsnc, "longitude")
nLat <- length(lat)
nLon <- length(lon)

# Make lon-lat-ens grid
lonlat <- expand.grid(lon = lon, lat = lat)

# Issue warnings
loss <- cbind(c(0,10,20,50),    # green column
              c(30,10,20,40),   # orange column
              c(80,50,20,40),   # red column
              c(100,100,80,20)) # black column

loss = t(loss)

getProb <- function(x){
  x1 <- mean(x<13.9)
  x2 <- mean(x>=13.9 & x<16.9)
  x3 <- mean(x>=16.9 & x <24.8)
  x4 <- mean(x>=24.8)
c(x1,x2,x3,x4)
}

warning <- function(x){
  p <- getProb(x)
  LP <- apply(loss*p,2,sum)
  which.min(LP)
}

warnings <- apply(preds,2,warning)
warnCols <- c("green","orange","red", "black")
breaks <- c(0.5,1.5,2.5,3.5,4.5)
MatrixWarn <- matrix(warnings,nrow=nLon,ncol=nLat)
x11(width=6,height=5)
par(mfrow=c(1,1),mar = c(4, 4, 1, 3),cex=1.3,lwd=1)
image(lon,lat,MatrixWarn,breaks=breaks,col=warnCols,xlab="",ylab="")
map(add=T, col="white")
title("Severe weather warnings")
dev.print(pdf,"warnings.pdf")