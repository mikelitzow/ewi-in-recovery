library(Rcpp)
library(ewidata)
library(knitr)
library(reshape2)
library(rstan)
library(tibble)
library(dplyr)
library(gtools)
library(bayesdfa)
library(assertthat)
library(rstanarm)
library(tidyverse)

# plot the 25-yr rolling window correlations between individual time series
# and shared DFA trends

# load the climate model
GOA.clim <- readRDS("GOA_clim_1_trend_original_data.rds") # read in GOA climate model

# recover the years and names for plots!
goa.clim <- read.csv("updated goa climate data.csv")

sub_data <- goa.clim %>%
  select(-ssh, -wind.stress) %>%
  gather(key="code", value="value", -year)

melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

# set new names
new.names <- c("East.spring.SST", "East.winter.SST", "SLP.gradient", "Papa.advection", "GAK1.salinity", "West.spring.SST", "West.winter.SST", "Downwell.54.134", "Downwell.57.137", "Downwell.60.146", "Downwell.60.149")

# look at rolling correlations between fitted and observed values, as in Ecology paper
# using 25-yr rolling windows!

# summarize various aspects of the model
modelfit <- GOA.clim$best_model
pred <- predicted(modelfit)
n_ts <- dim(modelfit$data)[1]
n_years <- dim(modelfit$data)[2]

# put them in a data frame with an inscrutable name
res.df <- data.frame(ID = rep(seq_len(n_ts),
                              n_years),
                     Time = sort(rep(seq_len(n_years), n_ts)), 
                     mean = c(t(apply(pred, c(3, 4), mean))), 
                     y = c(modelfit$data)) 

res.df$ID <- as.character(new.names)

res.df$year <-  rep(1950:2019, each=n_ts)   

# make a data frame to save results
dfa.cor <- data.frame()   

# now loop through and calculate the correlations
for(i in 1:length(new.names)){    
  # i <- 1 
  temp <- res.df %>%
    filter(ID==new.names[i])
  
  for(ii in 13:58){
  #  ii <- 13
    win <- temp[(ii-12):(ii+12),]
    temp.cor <- data.frame(ID=unique(win$ID),
                           year=win$year[13], 
                           window=25,
                           mean.cor=cor(win$y, win$mean, use="p"),
                           n.cor=sum(!is.na(win$y)))
    dfa.cor <- rbind(dfa.cor, temp.cor)
  }  }

# remove salinity correlations with < 10 obs
dfa.cor <- dfa.cor %>%
  filter(n.cor>=10)

# and plot 
# set palette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(dfa.cor, aes(year, mean.cor)) + 
  theme_bw() + theme(axis.title.x = element_blank()) +
  geom_line(color=cb[3]) + facet_wrap(~ID, 
                                      scales = "free_y") + xlab("Year") + ylab("Correlation") 

# and now the same for biology

# load in the model
GOA.biol <- readRDS("GOA_biol_11.4.19.rds") 

GOA <- read.csv("updated GOA biology data.csv")

# recover the years and names for plots!
sub_data = GOA

melted = melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names = Y$code
Y = as.matrix(Y[,-which(names(Y) == "code")])

# pretty names for plotting
new.names <- c("Chignik.coho", "Chignik.pink", "Chiniak.eulachon", "Chiniak.jellyfish", 
               "Chiniak.P.cod", "Chiniak.shrimp", "Cook.Inlet.coho", "Kodiak.coho", 
               "Kodiak.pink", "Pavlof.capelin", "Pavlof.eulachon", "Pavlof.jellyfish",
               "Pavlof.P.cod", "Pavlof.shrimp", "Prince.William.coho", "Southeast.coho", 
               "Southeast.herring", "Southeast.pink", "South.Peninsula.coho", "South.Peninsula.pink")

check <- data.frame(names=names, newnames=new.names)
check # looks right

# proceed as for climate
# summarize various aspects of the model
modelfit <- GOA.biol$best_model
pred <- predicted(modelfit)
n_ts <- dim(modelfit$data)[1]
n_years <- dim(modelfit$data)[2]

# put them in a data frame with an inscrutable name
res.df <- data.frame(ID = rep(seq_len(n_ts),
                              n_years),
                     Time = sort(rep(seq_len(n_years), n_ts)), 
                     mean = c(t(apply(pred, c(3, 4), mean))), 
                     y = c(modelfit$data)) 

res.df$ID <- as.character(new.names)

res.df$year <-  rep(1972:2019, each=n_ts)   

# make a data frame to save results
dfa.cor <- data.frame()   

# now loop through and calculate the correlations
for(i in 1:length(new.names)){    
  # i <- 1 
  temp <- res.df %>%
    filter(ID==new.names[i])
  
  for(ii in 13:36){
    #  ii <- 13
    win <- temp[(ii-12):(ii+12),]
    temp.cor <- data.frame(ID=unique(win$ID),
                           year=win$year[13], 
                           window=25,
                           mean.cor=cor(win$y, win$mean, use="p"),
                           n.cor=sum(!is.na(win$y)))
    dfa.cor <- rbind(dfa.cor, temp.cor)
  }  }


# and plot 

ggplot(dfa.cor, aes(year, mean.cor)) + 
  theme_bw() + theme(axis.title.x = element_blank()) +
  geom_line(color=cb[3]) + facet_wrap(~ID, 
                                      scales = "free_y") + xlab("Year") + ylab("Correlation") 
