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

# load updated ewidata!
ewidata <- read.csv("ewidata.csv", row.names=1)

# now load sewardline updates

sew <- read.csv("sewardline.updates.csv")

head(sew)
cor(sew[,2:9])

# and process
# I think for now we should just use the biomass variable rather than biomass and abundance
sew <- sew %>%
  select(1:5) %>%
  gather(key="code", value="value", -year)

ggplot(sew, aes(log(value))) +
  theme_bw() +
  geom_histogram(bins=12) +
  facet_wrap(~code, scales="free")

# pull old Sewardline data and replace
drop <- grep("SEWARD", ewidata$code)
ewidata <- ewidata[-drop,]

# and replace with new
sew$system <- "WGOA"
sew$subtype <- NA

# and add "SEWARD_" to codes
sew$code <- paste("SEWARD_", sew$code, sep="")
ewidata <- rbind(ewidata, sew)

# now add primary productivity data
pp <- read.csv("GOA-primary-production.csv")

pp <- pp %>%
  gather(key="code", value="value", -year)

ggplot(pp, aes(value)) +
  theme_bw() +
  geom_histogram(bins=8) +
  facet_wrap(~code, scales="free")

# I think we can live with that distribution 
# - although I should ask Eric/Sean if the bimodal distribution in duration is a problem for the DFA models

pp$code <- paste("PP_", pp$code, sep="")

pp$system <- "EGOA"
pp$subtype <- NA

ewidata <- rbind(ewidata, pp)
# 
# # now updated CPR data
# 
# filter(ewidata, code %in% c("copepod", "mesozoo", "diatoms"))

# GoA biology

codes <- unique(ewidata$code)

# limit to GOA
GOA <- ewidata[ewidata$system %in% c("WGOA", "EGOA"),]
unique(GOA$code)

# drop salmon
drop <- grep("PI", GOA$code)
GOA <- GOA[-drop,]

drop <- grep("CO", GOA$code)
GOA <- GOA[-drop,]

# drop climate
drop <- grep("AKCLIM", GOA$code)
GOA <- GOA[-drop,]

unique(GOA$code)

# drop Pavlof and Chiniak and herring and kittiwakes
drop <- grep("PAV", GOA$code)
GOA <- GOA[-drop,]

drop <- grep("CHI", GOA$code)
GOA <- GOA[-drop,]

drop <- grep("HER", GOA$code)
GOA <- GOA[-drop,]

drop <- grep("BLKI", GOA$code)
GOA <- GOA[-drop,]

# drop total icy zoops!
drop <- grep("ICY.ZOOP.TOTDEN", GOA$code)
GOA <- GOA[-drop,]

# and drop ichthyoplankton
drop <- grep("ICH", GOA$code)
GOA <- GOA[-drop,]

# and try dropping PP duration, as this is bimodal and might give trouble with the model 
# returning bimodal loadings
drop <- grep("_dur", GOA$code)
GOA <- GOA[-drop,]

unique(GOA$code)

# now check for outliers!
ggplot(GOA, aes(x=year, y = value))  + 
  theme_bw() + 
  geom_line() +
  facet_wrap(~code, scales = "free_y")  + xlab("Year") +
  ylab("Value")

# log transform ICY and Seward
trns <- c(grep("ICY", GOA$code), grep("SEWARD", GOA$code))
GOA$value[trns] <- log(GOA$value[trns])

# and check distributions
ggplot(GOA, aes(value)) + 
  theme_bw() + 
  geom_histogram() + facet_wrap(~code, scales="free")

sub_data = GOA

# set the name for the model output!
name <- "GOA_plankton_11.18.19.one.trends"


# and now run DFA!
mcmc_iter = 4000
max_trends = 1
mcmc_chains = 3
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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
