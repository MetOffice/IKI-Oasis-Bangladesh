#
# IKI Bangladesh (MIOASI): S8 Find gust conversion factor 
#
# Find gust conversion factor for converting 30 sec RA2 model gust
# to 3 sec gust, based on WMO guidance in:
# "Guidelines for Converting Between Various Wind Averaging Periods in Tropical Cyclone Conditions", 
# WMO Tech. Rep. No. 1555, 2010.
#
# Author: HS

library(ggplot2)
theme_set(theme_bw())
# Gust Speed Conversion Factors
refperiod <- c(3600, 600, 180, 120, 60, 3)
# Prediction Period
predperiod <- data.frame(refperiod = seq(30, 3600, 10))
# Known gust factors
il <- data.frame(refperiod, factor = c(1.75, 1.66, 1.58, 1.55, 1.49, 1))  # In land
ol <- data.frame(refperiod, factor = c(1.6, 1.52, 1.44, 1.42, 1.36, 1))  # Off land
os <- data.frame(refperiod, factor = c(1.45, 1.38, 1.31, 1.28, 1.23, 1))  # Off Sea
as <- data.frame(refperiod, factor = c(1.4, 1.23, 1.17, 1.15, 1.11, 1))  # At Sea

# fit GLM to log of refperiod
gf <- as
mod <- glm(factor ~ log(refperiod), data=gf)
pred <- predict(mod, newdata = predperiod, type = 'response', se.fit=T)
se <- pred$se.fit
sei <- qt(0.025, df = df.residual(mod), lower.tail = F) # 95% interval, 2.5% in each tail

ndata <- with(gf, predperiod)
ndata <- add_column(ndata, fit = pred$fit, hi = pred$fit + (sei * se), lo = pred$fit - (sei * se))
# Plot
plt <- ggplot(ndata, aes(x = refperiod, y = fit)) +
  geom_line() + geom_point(aes(y = factor), data = il) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.1) +
  scale_x_log10(limits=c(20, 4000)) +
  labs(x = 'Reference Time (sec)', y = 'Gust Factor')