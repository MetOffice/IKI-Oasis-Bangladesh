#
# IKI Bangladesh (MIOASI): S3 Sample TC event ensemble data
#
# Make netcdf of estimates of mean and credible intervals of gust speed. This
# process fits a GAM to the ensemble data, estimates a posterior distribution,
# and samples them to find the mean and credible interval based on a set of
# quantiles.
#
# Author: HS with TE
# Created: 12/8/19
# QA: KB (2/3/20) & TP (13/2/20)

library(RNetCDF)
library(reshape2)
library(fields)
library(evd)
library(mgcv)
library(plyr)
library(rogme)
library(abind)
library(doParallel)

# Set global variables
ROOTDIR = ""
STORMNAME <- 'VIYARU'
VAR <- 'fg.T1Hmax'
RES <- '4p4'
INDIR <- paste(ROOTDIR, "netcdf", sep = "")

USEWEIGHTS <- FALSE

QUANTILES <- c(0.5, 0.8, 0.9, 0.95, 0.96, 0.98, 0.99, 0.995)

registerDoParallel(cores = 4)

# Set lat-lon limits c(min, max)
LAT <- c(20.503502, 27.483002)
LON <- c(87.5555, 92.942)
num.knots <- 600

# Get weights
if (isTRUE(USEWEIGHTS)) {
  ensr <- read.csv(paste(ROOTDIR, "/fpens_r.csv", sep = ""))
  w <- subset(ensr, Name == STORMNAME)$r
} else {
  # Use equal weights
  w <- rep(1, 9)
}


# Load netcdf
filename <- paste('fpens.', VAR, '.', STORMNAME, '.', RES, 'km.nc', sep = '')
stormsnc <- open.nc(paste(INDIR, filename, sep = '/'))

# Extract netcdf variables
gust <- var.get.nc(stormsnc, 'wind_speed_of_gust') # an array of (810, 790, 9)
lat <- var.get.nc(stormsnc, 'latitude')
lon <- var.get.nc(stormsnc, 'longitude')

# Make lon-lat-ens grid
lonlatens <- expand.grid(lon, lat, seq(1, 9, 1))

# Make gust speed df
gustdf <- melt(gust, varnames = c('lon', 'lat', 'ens'), value.name = 'gust')
gustdf[1:3] <- lonlatens

# Extract area of interest to reduce data size
gustdf <- subset(gustdf, lon >= LON[1] & lon <= LON[2] & lat >= LAT[1] & lat <= LAT[2])
# Add weights information
gustdf$weights <- mapvalues(gustdf$ens,
                            from = seq_len(9),
                            to = w / mean(w))
# head(gustdf)
# Optional: plot distribution gust across space and ens. members
# plot(density(gustdf$gust))
# Quite skewed so go for a log-normal GAM
gustdf$lgust <- log(gustdf$gust)

# Set-up file names
if (isTRUE(USEWEIGHTS)) {
  outmodel <- paste(ROOTDIR, "fp.wmodel.", VAR, ".", STORMNAME, ".", RES, "km.Rdata", sep = "")
} else {
  outmodel <- paste(ROOTDIR, "fp.model.", VAR, ".", STORMNAME, ".", RES, "km.Rdata", sep = "")
}

# Check to see if GAM model file already exists, if not then run the model, otherwise load
if (!(file.exists(outmodel))) {
  # Run GAM: use Gaussian location scale model fittig across spatial lon-lat dimension
  # We use an isotropic smoother s() as we assume that smoothing will be similar in both
  # lon and lat dimensions (ie. equally as smooth in both dimensions)
  model <- gam(list(lgust ~ s(lon, lat, k=num.knots),
                    ~ s(lon, lat, k=num.knots)),
               family = gaulss(), method = "REML", weights = weights,
               data = gustdf)
  # Save the model
  save(model, file = outmodel)
  # Model checking
  print('Model Knots Check')
  printCoefmat(mgcv:::k.check(model, subsample = 5000, n.rep = 200), digits = 3)
  print('Model Summary')
  summary(model)
} else {
  # Load the model file
  load(outmodel)
}

### Get prediction intervals around the GAM predictions, using a Bayesian interpretation of the model
rmvn <- function(n,mu,sig) { ## MVN random deviates
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m*n),m,n))
}

### Simulate 1000 values from the posterior predictive distr. of gusts at each gridpoint
n.sims <- 1000
# 1000 replicate param. vectors
# betas <- rmvn(n.sims,coef(model),model$Vp)
betas <- rmvn(n.sims,coef(model),model$Vc) ## Vc version = Vp corrected covariance matrix for smoothing parameter uncertainty

# Define new dataframe with location data for predicitions, based on gustdf dataframe
newd <- gustdf[gustdf$ens == 1, 1:2]
# Make predictions
X <- predict(model,newdata = newd,type = "lpmatrix") ## model matrix

Mean <- X[, 1:num.knots] %*% t(betas[, 1:num.knots])  ## 1000 values from the posterior of the mean (of the Gaussian)
LPsd <- X[,(num.knots+1):(num.knots*2)] %*% t(betas[,(num.knots+1):(num.knots*2)])
Sd <- exp(LPsd) + 0.01 # as per ?gaulss

# simulate from predictive distribution (noting that the conditional is Log-Normal)
preds <- matrix(ncol = nrow(newd), nrow = n.sims)
for (simnum in 1:n.sims) {
    for (rownum in 1:nrow(newd)) {
      preds[simnum,rownum] <- rlnorm(1, meanlog = Mean[rownum,simnum], sdlog = Sd[rownum,simnum])
    }
  }


# Plot maps of mean estimate and upper, lower 95% pred. intervals
lons <- unique(newd$lon)
nLon <- length(lons)
lats <- unique(newd$lat)
nLat <- length(lats)
nCells <- nLon*nLat

# Calculate posterior stats
predMean <- apply(preds, 2, mean)
MatrixMean <- matrix(predMean, nrow = nLon, ncol = nLat)

# Use the simulated values at each gridpoint to calculate quantile prediction interval
# Do in parallel with %dopar% and write out to list qlist
qlist <-
  foreach(quant = seq_along(QUANTILES)) %dopar% {
    matrix( apply(preds, 2, hd, q = QUANTILES[quant]), nrow = nLon, ncol = nLat)
    }
# Convert list to 3D matrix
gustquantile <- abind(qlist, along = 0)

# Write out netCDF file
if (isTRUE(USEWEIGHTS)) {
  outfile <- paste(ROOTDIR, "fpw.", VAR, ".", STORMNAME, ".", RES, "km.nc", sep = "")
} else {
  outfile <- paste(ROOTDIR, "fp.", VAR, ".", STORMNAME, ".", RES, "km.nc", sep = "")
}

outnc <- create.nc(outfile, large = T)
# Dimensions
dim.def.nc(outnc, 'latitude', nLat)
dim.def.nc(outnc, 'longitude', nLon)
dim.def.nc(outnc, 'quantiles', length(QUANTILES))
# Variables
var.def.nc(outnc, 'latitude', 'NC_FLOAT', 'latitude')
var.def.nc(outnc, 'longitude', 'NC_FLOAT', 'longitude')
var.def.nc(outnc, 'quantiles', 'NC_FLOAT', 'quantiles')
var.def.nc(outnc, 'meangust', 'NC_FLOAT', c('longitude', 'latitude'))
var.def.nc(outnc, 'latitude_longitude', 'NC_CHAR', NA)
var.def.nc(outnc, 'gust_credible_interval', 'NC_FLOAT', c('quantiles', 'longitude', 'latitude'))

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
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lat_max", "NC_FLOAT", max(lats))
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lat_min", "NC_FLOAT", min(lats))
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lat_resolution", "NC_FLOAT", 0.04)
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lat_units", "NC_CHAR", "degrees_north")
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lon_max", "NC_FLOAT", max(lons))
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lon_min", "NC_FLOAT", min(lons))
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lon_resolution", "NC_FLOAT", 0.04)
att.put.nc(outnc, "NC_GLOBAL", "geospatial_lon_units", "NC_CHAR", "degrees_east")
att.put.nc(outnc, "NC_GLOBAL", "history", "NC_CHAR", "(1.0) Initial release")
att.put.nc(outnc, "NC_GLOBAL", "id", "NC_CHAR", paste("fp.", VAR, ".", STORMNAME, ".", RES, "km.nc", sep = ""))
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
var.put.nc(outnc, 'latitude', lats)
var.put.nc(outnc, 'longitude', lons)
var.put.nc(outnc, 'quantiles', QUANTILES)
var.put.nc(outnc, 'meangust', MatrixMean)
var.put.nc(outnc, 'gust_credible_interval', gustquantile)
close.nc(outnc)
