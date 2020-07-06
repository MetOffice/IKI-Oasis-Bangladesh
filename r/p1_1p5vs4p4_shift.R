#
# IKI Bangladesh (MIOASI): P1 1p5 vs 4p4 Shift Function
#
# This script compares the 4.4km and 1.5km model output, part of the IKI Bangladesh project. 
# As the data exist at different resolution, we compare distribution rather than absolute values.
#
# Author: HS

# Load libraries
library(RNetCDF)
library(reshape2)
library(fields)
library(evd)
library(plyr)
library(tibble)
library(rogme)
library(ggplot2)


# Load the data
#STORMNAME <- 'ROANU'
VAR <- 'fg.T1Hmax'
NCNAME <- 'wind_speed_of_gust'

# VAR <- 'psl.T1Hmin'
# NCNAME <- 'air_pressure_at_sea_level'

INDIR <- ''
# Set lat-lon limits c(min, max)
LAT <- c(20.503502, 25.3)
LON <- c(87.8, 92.942)

# Load netcdf
# filename <- paste('fpens.', VAR, '.', STORMNAME, '.4p4km.nc', sep = '')
filename <- paste('fpens.', VAR, '.grand.4p4km.nc', sep = '')
stormsnc <- open.nc(paste(INDIR, filename, sep = '/'))

# Extract netcdf variables
gust4p4 <- var.get.nc(stormsnc, NCNAME) # an array of (810, 790, 99)
lat <- var.get.nc(stormsnc, 'latitude')
lon <- var.get.nc(stormsnc, 'longitude')

# Make lon-lat-storm grid
lonlatst <- expand.grid(lon, lat, rep(1:13, each=9))

# Make gust speed df
gustdf <- melt(gust4p4, varnames = c('lon', 'lat', 'st'), value.name = 'gust4p4')
gustdf[1:3] <- lonlatst

# Extract area of interest to reduce data size
gustdf <- subset(gustdf, lon >= LON[1] & lon <= LON[2] & lat >= LAT[1] & lat <= LAT[2])

gc()

# Load netcdf 1.5km
# filename <- paste('fpens.', VAR, '.', STORMNAME, '.1p5km.nc', sep = '')
filename <- paste('fpens.', VAR, '.grand.1p5km.nc', sep = '')
stormsnc <- open.nc(paste(INDIR, filename, sep = '/'))

# Extract netcdf variables
gust1p5 <- var.get.nc(stormsnc, NCNAME)
lat <- var.get.nc(stormsnc, 'latitude')
lon <- var.get.nc(stormsnc, 'longitude')

# Make lon-lat-ens grid
lonlatst <- expand.grid(lon, lat, rep(1:13, each=9))

# Make gust speed df
gustdf2 <- melt(gust1p5, varnames = c('lon', 'lat', 'st'), value.name = 'gust1p5')
gustdf2[1:3] <- lonlatst

# Extract area of interest to reduce data size
gustdf2 <- subset(gustdf2, lon >= LON[1] & lon <= LON[2] & lat >= LAT[1] & lat <= LAT[2])

gc()



# Plot shift function
# ens = 5
# gust4p4 <- subset(gustdf, ens==ens)$gust4p4
# gust1p5 <- subset(gustdf2, ens==ens)$gust1p5
# df <- mkt2(gust4p4, gust1p5)
# sf <- shifthd_pbci(data = df, formula = obs ~ gr, nboot = 1000, q = seq(.05, .95, .05))
#
# #> plot shift function
# psf <- plot_sf(sf, plot_theme = 1)
# ggsave('/project/ciid/projects/IKI/1p4vs4p4.png', psf[[1]])


# Find heirarchical shift function, comparing 2 resolutions across 12 storms for as manay grid point samples as we have
gusttbl <- tibble(gust = c(gustdf2$gust1p5*1.94384, gustdf$gust4p4*1.94384),
             st = factor(c(gustdf2$st, gustdf$st)),
             res = factor(c(rep("1.5km",length(gustdf2$gust1p5)), rep("4.4km", length(gustdf$gust4p4))))
             )


h <- hsf_pb(gusttbl, gust ~ res + st, qseq=c(.01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .99))
save(h, file = paste("hsfplotting_", VAR, ".RData", sep = ''))

hsfplot <- plot_hsf_pb(h, ind_line_size = 0)
hsfplot$labels$y <- "Difference (knots)"
save(h, hsfplot, file = paste("hsfplotting_", VAR, ".RData", sep = ''))
ggsave(paste("1p4vs4p4_hsf_", VAR, ".png", sep = ''), hsfplot)

plot_hsf_pb_dist <- function(data,
                             null_value = 0,
                             point_interv = "mode_hdi",
                             interval_width = c(0.5, 0.9),
                             fill_colour = "orange",
                             int_colour = "black"
){
  # check input is a list of data frames
  if(!is.list(data)){
    stop("data must be a list returned by hsf_pb()")
  }
  if(length(data)!=8){
    stop("data must have length 8, as returned by hsf_pb()")
  }

  df <- tibble(boot_samp = as.vector(data$bootstrap_samples),
               quantile = rep(data$quantiles, each = data$nboot))
  p <- ggplot(df, aes(x = boot_samp, y = quantile)) +
    theme_classic() +
    tidybayes::geom_halfeyeh(
      fill = fill_colour,
      color = int_colour,
      .point_interval = point_interv,
      .width = interval_width) +
    geom_vline(xintercept = null_value) +
    scale_y_continuous(breaks = data$quantiles) +
    theme(plot.title = element_text(size=22),
          axis.title.x = element_text(size = 18),
          axis.text = element_text(size = 16, colour = "black"),
          axis.title.y = element_text(size = 18)) +
    labs(y = "Quantiles", x = "Bootstrap differences") +
    ggtitle(data$comparison) +
    coord_flip()

  suppressMessages(p)
}

dist <- plot_hsf_pb_dist(h, point_interv = 'median')
dist$labels$y <- "1.5km percentiles"
dist$labels$x <- "Difference (knots)"
dist + scale_y_discrete(breaks=c("0.01", "0.05", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "0.95", "0.99"),
                        labels=c("1", "5", "10", "20", "30", "40", "50", "60", "70", "80", "90", "95", "99"))
