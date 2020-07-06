#
# IKI Bangladesh (MIOASI): S6 Exceedence Probabilities
#
# Make netcdf of estimates of mean, and low and upper credible intervals of mean gust speed.
# This process fits a GAM to the ensemble data, estimates a posterior distribution,
# and samples them to find the mean and credible interval based on a set of quantiles.
#
#
# Author: HS
# Created: Jan 2020
# QA: TE 17/3/20

library(RNetCDF)
library(abind)
library(doParallel)
registerDoParallel(cores = 10)

# Define Bangladesh tropical cyclon categories in knots converted to m/s
# WINDTHRESHOLDS <- c(17, 22, 28, 34, 48, 64, 120) * 0.514444
WINDTHRESHOLDS <- seq(0, 100, 1)



INDIR <- ""
VAR <- 'fg.T1Hmax'
RES <- '4p4'
# Load base netcdf to get lat/lon variables
filename <- paste('fp.', VAR, '.AILA.', RES, 'km.nc', sep = '')
stormsnc <- open.nc(paste(INDIR, filename, sep = '/'))

# Extract netcdf variables
lat <- var.get.nc(stormsnc, 'latitude')
lon <- var.get.nc(stormsnc, 'longitude')

nLon <- length(lon)
nLat <- length(lat)
nCells <- nLon*nLat

# Make lon-lat-ens grid
lonlat <- expand.grid(lon = lon, lat = lat)

# Load predictions data from RData object
predsfile <- "grand_fp_preds.Rdata"
load(predsfile)

# Calculate posterior stats from preds
# Use the simulated values at each gridpoint to calculate exceedend probabilities
# Do in parallel with %dopar% and write out to list qlist
plist <- foreach(i = seq_along(WINDTHRESHOLDS)) %dopar% {
  # Make bool array of T/F values based on wind thresholds
  bool <- preds >= WINDTHRESHOLDS[i]
  # Find threshold probabilities based on total True values / total number of preds -> mean
  # N.B. Here the mean is the Monte Carlo approximation to the integral -- so unrelated to the mean of any distribution
  p <- apply(bool, 2, mean)
  matrix(p, nrow = nLon, ncol = nLat)
}
# Convert list to 3D matrix
exceedm <- abind(plist, along=0)

# Write out netCDF file
outfile <- paste('pgrand.exceedance_tmp.', VAR, '.', RES, 'km.nc', sep = '')
outnc <- create.nc(outfile, large=T)
# Dimensions
dim.def.nc(outnc, 'latitude', nLat)
dim.def.nc(outnc, 'longitude', nLon)
dim.def.nc(outnc, 'gust_threshold', length(WINDTHRESHOLDS))
# Variables
var.def.nc(outnc, 'latitude', 'NC_FLOAT', 'latitude')
var.def.nc(outnc, 'longitude', 'NC_FLOAT', 'longitude')
var.def.nc(outnc, 'gust_threshold', 'NC_INT', 'gust_threshold')
var.def.nc(outnc, 'latitude_longitude', 'NC_CHAR', NA)

var.def.nc(outnc, 'probability_of_exceedance', 'NC_FLOAT', c('gust_threshold', 'longitude', 'latitude'))

# Define some attributes
att.put.nc(outnc, "latitude", "axis", "NC_CHAR", "Y")
att.put.nc(outnc, "latitude", "units", "NC_CHAR", "degrees_north")
att.put.nc(outnc, "latitude", "standard_name", "NC_CHAR", "latitude")

att.put.nc(outnc, "longitude", "axis", "NC_CHAR", "X")
att.put.nc(outnc, "longitude", "units", "NC_CHAR", "degrees_east")
att.put.nc(outnc, "longitude", "standard_name", "NC_CHAR", "longitude")

att.put.nc(outnc, "gust_threshold", "long_name", "NC_CHAR", "Gust speed classification threshold")
att.put.nc(outnc, "gust_threshold", "units", "NC_CHAR", "m s-1")

att.put.nc(outnc, "probability_of_exceedance", "long_name", "NC_CHAR", "posterior wind_speed_of_gust probability_of_exceedance")
att.put.nc(outnc, "probability_of_exceedance", "units", "NC_INT", 1)
att.put.nc(outnc, "probability_of_exceedance", "grid_mapping", "NC_CHAR", "latitude_longitude")

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
att.put.nc(outnc, "NC_GLOBAL", "id", "NC_CHAR", paste("fpgrandw.", VAR, ".", RES, "km.nc", sep=""))
att.put.nc(outnc, "NC_GLOBAL", "institution", "NC_CHAR", "Met Office, UK")
att.put.nc(outnc, "NC_GLOBAL", "licence", "NC_CHAR", "Creative Commons Attribution 4.0 International (CC BY 4.0)")
att.put.nc(outnc, "NC_GLOBAL", "product_version", "NC_CHAR", "v1.0")
att.put.nc(outnc, "NC_GLOBAL", "project", "NC_CHAR", "Oasis Platform for Climate and Catastrophe Risk Assessment – Asia")
att.put.nc(outnc, "NC_GLOBAL", "institution", "NC_CHAR", "Met Office, UK")
att.put.nc(outnc, "NC_GLOBAL", "keywords", "NC_CHAR", "Bangladesh, footprint, quantiles, Met Office")
att.put.nc(outnc, "NC_GLOBAL", "source", "NC_CHAR", "Met Office UM RA2T CON")
att.put.nc(outnc, "NC_GLOBAL", "spatial_resolution", "NC_CHAR", "4.4km")
att.put.nc(outnc, "NC_GLOBAL", "summary", "NC_CHAR", "Tropical cyclone gust speed threshold exceedence probabilities over Bangladesh")
att.put.nc(outnc, "NC_GLOBAL", "keywords", "NC_CHAR", "Bangladesh, footprint, posterior, exceedence probabilities, Met Office")

# Write in data
var.put.nc(outnc, 'latitude', lat)
var.put.nc(outnc, 'longitude', lon)
var.put.nc(outnc, 'gust_threshold', WINDTHRESHOLDS)
var.put.nc(outnc, 'probability_of_exceedance', exceedm)
close.nc(outnc)