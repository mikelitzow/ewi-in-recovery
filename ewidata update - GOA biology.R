library(tidyverse)
library(ncdf4)
library(chron)
library(fields)
library(maps)
library(mapdata)
library(zoo)

dat <- read.csv("ewidata.csv", row.names = 1)

levels(dat$system)

# restrict to GOA 

goa.dat <- dat %>%
  filter(system %in% c("EGOA", "WGOA"))

# drop climate data!
drop <- grep("AKC", goa.dat$code)
goa.dat <- goa.dat[-drop,]

ggplot(goa.dat, aes(year, value)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  facet_wrap(~code, scales="free_y") 

# examine salmon timeseries alone

salmon <- c(grep(".PI", goa.dat$code), grep(".CO", goa.dat$code))

salm.dat <- goa.dat[salmon,]

ggplot(salm.dat, aes(year, value)) +
  theme_bw() +
  geom_line() +
  geom_point() +
  facet_wrap(~code, scales="free_y") 

# update