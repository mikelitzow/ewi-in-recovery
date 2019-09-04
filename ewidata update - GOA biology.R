library(tidyverse)
library(ncdf4)
library(chron)
library(fields)
library(maps)
library(mapdata)
library(zoo)

# load updated ewidata!
ewidata <- read.csv("ewidata.csv", row.names=1)

codes <- unique(ewidata$code)

# limit to GOA
GOA <- ewidata[ewidata$system %in% c("WGOA", "EGOA"),]
unique(GOA$code)

# drop CI pinks as they aren't present in odd and even years
GOA <- filter(GOA, code != "CI.PI")

# drop climate
drop <- grep("AKCLIM", GOA$code)
GOA <- GOA[-drop,]

unique(GOA$code)

# drop dogfish (too many NAs)
drop <- grep("DOGF", GOA$code)
GOA <- GOA[-drop,]

# drop birds and icy strait and chatnam and seward line and CPR
drop <- c(grep("PROD.", GOA$code), grep("PHEN.", GOA$code), grep("ICY", GOA$code), grep("CHAT", GOA$code), grep("SEWARD", GOA$code), 
          grep("mesozoo", GOA$code), grep("diatoms", GOA$code), grep("copepod", GOA$code))
GOA <- GOA[-drop,]

# herring time series are repeats, remove
drop <- grep("SE.HER.ROE", GOA$code)
GOA <- GOA[-drop,]

# and drop the ichthyo data
drop <- grep("ICH", GOA$code)
GOA <- GOA[-drop,]

# qnd, finally, drop PWS pinks as they are so hatchery-dominated
GOA <- GOA %>%
  filter(code != "PWS.PI")

# and start in 1972
GOA <- GOA[GOA$year>=1972,]

unique(GOA$code)

ggplot(GOA, aes(year, value)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  facet_wrap(~code, scales="free_y") 

# update with latest salmon values
update <- read.csv("salmon catch updates.csv")
update$system <- update$subtype <- NA

GOA <- rbind(GOA, update)

# now scale odd and even pink catches for each area

# first, remove pinks from GOA

pinks <- grep(".PI", GOA$code)

pink.dat <- GOA[pinks,]
GOA <- GOA[-pinks,]

# check
ggplot(pink.dat, aes(year, value)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  facet_wrap(~code, scales="free_y")

# scale area-by-area
SE.O <- pink.dat %>%
  filter(code=="SE.PI", odd(year))
SE.O$value <- scale(SE.O$value)

SE.E <- pink.dat %>%
  filter(code=="SE.PI", even(year))
SE.E$value <- scale(SE.E$value)

CH.O <- pink.dat %>%
  filter(code=="CH.PI", odd(year))
CH.O$value <- scale(CH.O$value)

CH.E <- pink.dat %>%
  filter(code=="CH.PI", even(year))
CH.E$value <- scale(CH.E$value)

KOD.O <- pink.dat %>%
  filter(code=="KOD.PI", odd(year))
KOD.O$value <- scale(KOD.O$value)

KOD.E <- pink.dat %>%
  filter(code=="KOD.PI", even(year))
KOD.E$value <- scale(KOD.E$value)

SP.O <- pink.dat %>%
  filter(code=="SP.PI", odd(year))
SP.O$value <- scale(SP.O$value)

SP.E <- pink.dat %>%
  filter(code=="SP.PI", even(year))
SP.E$value <- scale(SP.E$value)

scaled.pink <- rbind(SE.O, SE.E, CH.O, CH.E, KOD.O, KOD.E, SP.O, SP.E)

# check
ggplot(scaled.pink, aes(year, value)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  facet_wrap(~code, scales="free_y")

# put back into GOA
GOA <- rbind(GOA, scaled.pink)

# and save
write.csv(GOA, "updated GOA biology data.csv", row.names = F)

# and check one more time
ggplot(GOA, aes(year, value)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  facet_wrap(~code, scales="free_y")

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
name <- "GOA_biol_9.4.19"

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


