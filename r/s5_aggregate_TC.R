#
# IKI Bangladesh (MIOASI): S5 Aggregate TC event ensemble data
#
# This script aggregate all individual TC GAM models into a grandmodel for
# looking at combined TC hazard within our 99 member ensmeble.
#
# This loads each individual RData file and pools the posterior samples, then
# calculates statistics as in s3 and saves a single netCDFfile.
#
#
# Author: HS
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

USEWEIGHTS <- FALSE

QUANTILES <- c(0.01, 0.05, 0.10, 0.5, 0.9, 0.95, 0.99)

registerDoParallel(cores = 4)

n.sims <- 1000  # No. of simulations from posterior predictive distribution

# Set lat-lon limits c(min, max)
LAT <- c(20.503502, 27.483002)
LON <- c(87.5555, 92.942)
num.knots <- 600


# Load base netcdf
filename <- paste("fp.", VAR, ".", STORMNAMES[1], ".", RES, "km.nc", sep = "")
stormsnc <- open.nc(paste(INDIR, filename, sep = "/"))

# Extract netcdf variables
lat <- var.get.nc(stormsnc, "latitude")
lon <- var.get.nc(stormsnc, "longitude")

# Make lon-lat-ens grid
lonlat <- expand.grid(lon = lon, lat = lat)

### Get prediction intervals around the GAM predictions, using a Bayesian interpretation of the model
rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig);m <- ncol(L);
  t(mu + L %*% matrix(rnorm(m*n),m,n))
}

# Set-up out file
outfile <- "grand_fp_preds.Rdata"

if (!(file.exists(outfile))) {
  # Loop over storms
  preds <- matrix(ncol = nrow(lonlat),nrow = n.sims*length(STORMNAMES))
  for (s in seq_along(STORMNAMES)) {
    # Get RData file
    if (isTRUE(USEWEIGHTS)) {
      inmodel <-
        paste(
          "/fp.wmodel.",
          VAR,
          ".",
          STORMNAMES[s],
          ".",
          RES,
          "km.Rdata",
          sep = ""
        )
    } else {
      inmodel <-
        paste(
          "fp.model.",
          VAR,
          ".",
          STORMNAMES[s],
          ".",
          RES,
          "km.Rdata",
          sep = ""
        )
    }
      # Load RData file, this loads a 'model' variable
      load(inmodel)

      betas <- rmvn(n.sims, coef(model), model$Vc)
      X <- predict(model, newdata = lonlat,type = "lpmatrix")
      Mean <- X[,1:num.knots] %*% t(betas[,1:num.knots])  ## 1000 values from the posterior of the mean (of the Gaussian)
      LPsd <- X[,(num.knots+1):(num.knots*2)] %*% t(betas[,(num.knots+1):(num.knots*2)])
      Sd <- exp(LPsd + 0.01) # as per ?gaulss

      # Make sequence of indicies for the large pred array for each loop along STORMNAMES
      # ie. first storm = i{1,1000}
      # second storm = i{1001, 2000} etc.
      sim.num <- seq(n.sims*(s-1)+1, n.sims*s)
      # simulate from predictive distribution (noting that the conditional is Log-Normal)
      for (i in seq_along(sim.num)) {
        for (j in 1:nrow(lonlat)) {
          preds[sim.num[i],j] <- rlnorm(1,meanlog = Mean[j,i],sdlog = Sd[j,i])
        }}
      # Do garbage collection
      rm(model)
      gc()
  }
  # Save predictions
  save(preds, file = outfile)
} else {
  # Load the model file
  load(outfile)
}


# Plot maps of mean estimate and upper, lower 95% pred. intervals
nLon <- length(lon)
nLat <- length(lat)
nCells <- nLon*nLat

# Calculate posterior stats
Mean <- apply(preds,2,mean)
MatrixMean <- matrix(Mean,nrow = nLon, ncol = nLat)

# Use the simulated values at each gridpoint to calculate 95% prediction interval
# Do in parallel with %dopar% and write out to list qlist
# Take transpose of matrix to ensure that lat-lon order is correct for writing to netcdf files later
qlist <-
  foreach(i = seq_along(QUANTILES)) %dopar% {
    matrix( apply(preds,2,hd, q = QUANTILES[i]), nrow = nLon, ncol = nLat)
  }

# Convert list to 3D matrix
gustquantile <- abind(qlist, along = 0)


# Write out netCDF file
if (isTRUE(USEWEIGHTS)) {
  outfile <-
    paste(
      "fpgrandw.",
      VAR,
      ".",
      RES,
      "km.nc",
      sep = ""
    )
} else {
  outfile <-
    paste(
      "fpgrand.",
      VAR,
      ".",
      RES,
      "km.nc",
      sep = ""
    )

}

# Make netCDF file
outnc <- create.nc(outfile, large = T)
# Dimensions
dim.def.nc(outnc, "latitude", nLat)
dim.def.nc(outnc, "longitude", nLon)
dim.def.nc(outnc, "quantiles", length(QUANTILES))
# Variables
var.def.nc(outnc, "latitude", "NC_FLOAT", "latitude")
var.def.nc(outnc, "longitude", "NC_FLOAT", "longitude")
var.def.nc(outnc, "quantiles", "NC_FLOAT", "quantiles")
var.def.nc(outnc, "meangust", "NC_FLOAT", c("longitude", "latitude"))
var.def.nc(outnc, "latitude_longitude", "NC_CHAR", NA)

var.def.nc(outnc, "gust_credible_interval", "NC_FLOAT", c("quantiles", "longitude", "latitude"))

# Define some attributes
att.put.nc(outnc, "latitude", "axis", "NC_CHAR", "Y")
att.put.nc(outnc, "latitude", "units", "NC_CHAR", "degrees_north")
att.put.nc(outnc, "latitude", "standard_name", "NC_CHAR", "latitude")

att.put.nc(outnc, "longitude", "axis", "NC_CHAR", "X")
att.put.nc(outnc, "longitude", "units", "NC_CHAR", "degrees_east")
att.put.nc(outnc, "longitude", "standard_name", "NC_CHAR", "longitude")

att.put.nc(outnc, "quantiles", "long_name", "NC_CHAR", "prediction quantiles")

att.put.nc(outnc, "meangust", "long_name", "NC_CHAR", "posterior mean wind_speed_of_gust")
att.put.nc(outnc, "meangust", "units", "NC_CHAR", "m s-1")
att.put.nc(outnc, "meangust", "um_stash_source", "NC_CHAR", "m01s03i463")
att.put.nc(outnc, "meangust", "standard_name", "NC_CHAR", "wind_speed_of_gust")
att.put.nc(outnc, "meangust", "grid_mapping", "NC_CHAR", "latitude_longitude")

att.put.nc(outnc, "gust_credible_interval", "long_name", "NC_CHAR", "posterior credible intervals of wind_speed_of_gust")
att.put.nc(outnc, "gust_credible_interval", "units", "NC_CHAR", "m s-1")
att.put.nc(outnc, "gust_credible_interval", "um_stash_source", "NC_CHAR", "m01s03i463")
att.put.nc(outnc, "gust_credible_interval", "grid_mapping", "NC_CHAR", "latitude_longitude")

att.put.nc(outnc, "latitude_longitude", "grid_mapping_name", "NC_CHAR", "latitude_longitude")
att.put.nc(outnc, "latitude_longitude", "longitude_of_prime_meridian", "NC_INT", 0)
att.put.nc(outnc, "latitude_longitude", "earth_radius", "NC_INT", 6371229.)
att.put.nc(outnc, "latitude_longitude", "proj4", "NC_CHAR", "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

att.put.nc(outnc, "NC_GLOBAL", "comment", "NC_CHAR", "Supported by the International Climate Initiative (IKI) and the Federal Ministry for the Environment, Nature Conservation and Nuclear Safety, based on a decision of the Germany Bundestag")
att.put.nc(outnc, "NC_GLOBAL", "contact", "NC_CHAR", "enquiries@metoffice.gov.uk")
att.put.nc(outnc, "NC_GLOBAL", "data_type", "NC_CHAR", "grid")
att.put.nc(outnc, "NC_GLOBAL", "date_created", "NC_CHAR", format(Sys.time(), "%Y%m%dT%H:%M:%S"))
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lat_max", "NC_FLOAT", max(lat))
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lat_min", "NC_FLOAT", min(lat))
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lat_resolution", "NC_FLOAT", 0.04)
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lat_units", "NC_CHAR", "degrees_north")
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lon_max", "NC_FLOAT", max(lon))
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lon_min", "NC_FLOAT", min(lon))
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lon_resolution", "NC_FLOAT", 0.04)
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lon_units", "NC_CHAR", "degrees_east")
att.put.nc(outnc, "NC_GLOBAL", "history", "NC_CHAR", "(1.0) Initial release")
att.put.nc(outnc, "NC_GLOBAL", "id", "NC_CHAR", paste("fpgrandw.", VAR, ".", RES, "km.nc", sep = ""))
att.put.nc(outnc, "NC_GLOBAL", "institution", "NC_CHAR", "Met Office, UK")
att.put.nc(outnc, "NC_GLOBAL", "licence", "NC_CHAR", "Creative Commons Attribution 4.0 International (CC BY 4.0)")
att.put.nc(outnc, "NC_GLOBAL", "product_version", "NC_CHAR", "v1.0")
att.put.nc(outnc, "NC_GLOBAL", "project", "NC_CHAR", "Oasis Platform for Climate and Catastrophe Risk Assessment – Asia")
att.put.nc(outnc, "NC_GLOBAL", "institution", "NC_CHAR", "Met Office, UK")
att.put.nc(outnc, "NC_GLOBAL", "keywords", "NC_CHAR", "Bangladesh, footprint, quantiles, Met Office")
att.put.nc(outnc, "NC_GLOBAL", "source", "NC_CHAR", "Met Office UM RA2T CON")
att.put.nc(outnc, "NC_GLOBAL", "spatial_resolution", "NC_CHAR", "4.4km")
att.put.nc(outnc, "NC_GLOBAL", "summary", "NC_CHAR", "Tropical cyclone footprint mean and prediction quantiles over Bangladesh")
att.put.nc(outnc, "NC_GLOBAL", "keywords", "NC_CHAR", "Bangladesh, footprint, posterior, credible intervals, Met Office")

# Write in data
var.put.nc(outnc, "latitude", lat)
var.put.nc(outnc, "longitude", lon)
var.put.nc(outnc, "quantiles", QUANTILES)
var.put.nc(outnc, "meangust", MatrixMean)
var.put.nc(outnc, "gust_credible_interval", gustquantile)
close.nc(outnc)