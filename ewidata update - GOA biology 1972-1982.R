library(tidyverse)
library(ncdf4)
library(chron)
library(fields)
library(maps)
library(mapdata)
library(zoo)
library(gtools)


GOA <- read.csv("updated GOA biology data.csv")

GOA <- GOA %>%
  filter(year %in% 1972:1982, code != "SE.HER")

# and check one more time
ggplot(GOA, aes(year, value)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  facet_wrap(~code, scales="free_y")

##########################################
# now fit models!
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
max_trends = 1 # changing to 1!
mcmc_chains = 3
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

sub_data = GOA

# set the name for the model output!
name <- "GOA_biol_1972-1982"

# reshape data
melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

# do the trend search, and save the table of model selection, along with the best model. By default, this isn't comparing student-t
# versus normal models, but just estimating the student-t df parameter
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

rotated = rotate_trends(dfa_summary$best_model)

# run HMM test!
y = apply(rotated$trends, c(2,3), mean)
sd_y = apply(rotated$trends, c(2,3), sd)

# run each trend sequentially through the regimes
# code to identify (1) number of changes and (2)
# where those change points occur.


set.seed(99)
f1 = find_regimes(y = y[1,], sds = sd_y[1,], max_regimes = 2)
print(f1$table) # this shows 2-regime model is best
plot_regime_model(f1$best_model)
