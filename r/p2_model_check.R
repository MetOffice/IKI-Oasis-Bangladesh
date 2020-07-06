#
# IKI Bangladesh (MIOASI): P2 Model Check plots
#
# Plot some basic model checking plots for each GAM model
#
# Author: HS

library(mgcViz)
library(doParallel)
registerDoParallel(cores = 6)

INDIR = ""
OUTDIR = ""
STORMNAMES <- c("AILA", "AKASH", "BOB01", "BOB07", "BULBUL", "FANI", "MORA", "RASHMI", "ROANU", "SIDR", "TC01B", "VIYARU")

foreach(i = seq_along(STORMNAMES)) %dopar% {
  load(paste(INDIR, "fp.model.fg.T1Hmax.", STORMNAMES[i], ".4p4km.Rdata", sep=""))
  b <- getViz(model, nsim = 100, post = TRUE, unconditional = TRUE)
  # b contains nsim vectors of n=100 responses simulated from the predictive posterior
  # Compare the empirical distribution of the observed residuals with that of the simulated ones
  distplt <- check0D(b) + l_dens1D()
  # Plot equivalent as a QQ plot
  qqplt <- qq(b, method = "simul1", a.qqpoi = list("shape" = 1), a.ablin = list("linetype" = 2), ngr=50, xlim=c(-3.1,3.1), ylim=c(-2,2), worm=TRUE)
  qqplt <- qqplt + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # Optionally, check kutosis
  # check0D(b, trans = function(.y) mean((.y - mean(.y))^4)) + l_dens1D() + l_vline()
  # or skewness
  # check0D(b, trans = function(.y) mean((.y - mean(.y))^3)) + l_dens1D() + l_vline()
  # or variance
  # check0D(b, trans = function(.y) mean((.y - mean(.y))^2)) + l_dens1D() + l_vline()
  ggsave(paste(OUTDIR, "qqplot_", STORMNAMES[i], ".ps"), plot = qqplt$ggObj, width = 3, height = 3, units = "in")
  ggsave(paste(OUTDIR, "distplot_", STORMNAMES[i], ".png"), plot = distplt$ggObj)
  rm(model)
  gc()
}

