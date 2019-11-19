library(tidyverse)
library(ncdf4)
library(chron)
library(fields)
library(maps)
library(mapdata)
library(zoo)
library(gtools)

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

# log transform
sew$value <- log(sew$value)

# pull old Sewardline data and replace
drop <- grep("SEWARD", ewidata$code)
ewidata <- ewidata[-drop,]

# and replace with new
sew$system <- sew$subtype <- NA

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

pp$system <- pp$subtype <- NA

ewidata <- rbind(ewidata, pp)

# now updated CPR data

filter(ewidata, code %in% c("copepod", "mesozoo", "diatoms"))
