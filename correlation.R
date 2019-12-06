# library(Rcpp)
# library(ewidata)
# library(knitr)
library(reshape2)
library(rstan)
# library(tibble)
library(dplyr)
# library(gtools)
library(bayesdfa)
# library(assertthat)
# library(rstanarm)
library(tidyverse)

library(future)
plan(multiprocess)

# plot the 25-yr rolling window correlations between individual time series
# and shared DFA trends

# load the climate model
GOA.clim <- readRDS("GOA_clim_1_trend_original_data.rds") # read in GOA climate model

# recover the years and names for plots!
goa.clim <- read.csv("updated goa climate data.csv")

sub_data <- goa.clim %>%
  select(-ssh, -wind.stress) %>%
  gather(key = "code", value = "value", -year)

melted <- melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names <- Y$code
Y <- as.matrix(Y[, -which(names(Y) == "code")])

# set new names
new.names <- c("East.spring.SST", "East.winter.SST", "SLP.gradient",
  "Papa.advection", "GAK1.salinity", "West.spring.SST", "West.winter.SST",
  "Downwell.54.134", "Downwell.57.137", "Downwell.60.146", "Downwell.60.149")

# look at rolling correlations between fitted and observed values, as in Ecology paper
# using 25-yr rolling windows!

# summarize various aspects of the model
modelfit <- GOA.clim$best_model
pred <- predicted(modelfit)

stan_correlation_m <- stan_model("correlation.stan")
fit_cor_stan <- function(.x) {
  message(paste0(
    "Variable: ", .x$variable[1],
    ", Chain: ", .x$chain[[1]], ", Iteration: ",
    .x$iteration[[1]]
  ))
  stan_dat <- list(x = cbind(.x$trend_value, .x$data_value), N = nrow(.x))
  m <- sampling(stan_correlation_m,
    data = stan_dat,
    iter = 500, chains = 1, verbose = FALSE, refresh = 0,
    seed = 1947823,
  )
  rho_post <- rstan::extract(m)$rho
  tibble(
    rho_post = rho_post, year = .x$year[1],
    variable = .x$variable[1], orig_iteration = .x$iteration[1],
    chain = .x$chain[1],
    ess = summary(m)$summary[, "n_eff"][["rho"]]
  )
}

pred_long <- reshape2::melt(pred) %>% as_tibble()
names(pred_long) <- c("iteration", "chain", "year", "variable", "trend_value")
set.seed(83874)
pred_long <- pred_long %>%
  group_by(variable, chain) %>%
  filter(iteration %in% sample(unique(pred_long$iteration), 20)) %>%
  ungroup() %>%
  mutate(year = year + 1949)
orig_data <- reshape2::melt(modelfit$data) %>%
  as_tibble() %>%
  rename(variable = Var1, year = Var2, data_value = value)
orig_and_pred <- inner_join(pred_long, orig_data,
  by = c("year", "variable")
) %>%
  filter(!is.na(data_value))

fit_cor_window <- function(begin, middle, end) {
  filter(orig_and_pred, year %in% seq(begin, end)) %>%
    group_split(variable, chain, iteration) %>%
    future.apply::future_lapply(fit_cor_stan) %>%
    dplyr::bind_rows() %>%
    mutate(middle = middle, begin = begin, end = end)
}

middle <- seq(1950 + 12, 1950 + 57, 5)
end <- middle + 12
begin <- middle - 12
stopifnot(min(orig_and_pred$year) == min(begin))
stopifnot(max(orig_and_pred$year) == max(end))
years_df <- data.frame(begin = begin, middle = middle, end = end)
system.time({
  out <- purrr::pmap_dfr(years_df, fit_cor_window)
})
saveRDS(out, "stan-cor-window-env.rds")
out <- readRDS("stan-cor-window-env.rds")

out <- left_join(out, tibble(
  var_name = gsub("\\.", " ", new.names),
  variable = 1:length(new.names)
))
out$var_name <- factor(out$var_name, levels = gsub("\\.", " ", new.names))

out %>%
  group_by(var_name, middle) %>%
  summarise(
    lwr = quantile(rho_post, probs = 0.1),
    lwr2 = quantile(rho_post, probs = 0.25),
    upr = quantile(rho_post, probs = 0.9),
    upr2 = quantile(rho_post, probs = 0.75),
    med = quantile(rho_post, probs = 0.5)
  ) %>%
  # filter(!(middle < 1970 & var_name == "GAK1 salinity")) %>%
  ggplot(aes(x = middle, y = med, ymin = lwr, ymax = upr)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.2) +
  geom_line(colour = "#56B4E9", lwd = .7) +
  geom_ribbon(
    alpha = 0.4, fill = "#56B4E9",
    mapping = aes(ymin = lwr2, ymax = upr2)
  ) +
  geom_ribbon(alpha = 0.25, fill = "#56B4E9") +
  facet_wrap(~var_name, scales = "fixed", ncol = 3) +
  theme_bw() +
  ylab("Correlation") +
  xlab("Year (middle of 25-year window)") +
  theme(axis.title.x = element_blank()) +
  coord_cartesian(ylim = c(-0.25, 1), expand = FALSE)
dir.create("figs", showWarnings = FALSE)
ggsave("figs/stan-cor-window-env.png", width = 5.5, height = 5.5)

# and now the same for biology ------------------------
# load in the model
GOA.biol <- readRDS("GOA_biol_11.4.19.rds")
GOA <- read.csv("updated GOA biology data.csv")

# recover the years and names for plots!
sub_data <- GOA
melted <- melt(sub_data[, c("code", "year", "value")], id.vars = c("code", "year"))
Y <- dcast(melted, code ~ year)
names <- Y$code
Y <- as.matrix(Y[, -which(names(Y) == "code")])

# pretty names for plotting
new.names <- c(
  "Chignik.coho", "Chignik.pink", "Chiniak.eulachon", "Chiniak.jellyfish",
  "Chiniak.P.cod", "Chiniak.shrimp", "Cook.Inlet.coho", "Kodiak.coho",
  "Kodiak.pink", "Pavlof.capelin", "Pavlof.eulachon", "Pavlof.jellyfish",
  "Pavlof.P.cod", "Pavlof.shrimp", "Prince.William.coho", "Southeast.coho",
  "Southeast.herring", "Southeast.pink", "South.Peninsula.coho", "South.Peninsula.pink"
)

check <- data.frame(names = names, newnames = new.names)
check # looks right

# # proceed as for climate
# # summarize various aspects of the model
modelfit <- GOA.biol$best_model
pred <- predicted(modelfit)

pred_long <- reshape2::melt(pred) %>% as_tibble()
names(pred_long) <- c("iteration", "chain", "year", "variable", "trend_value")
set.seed(83874)
pred_long <- pred_long %>%
  group_by(variable, chain) %>%
  filter(iteration %in% sample(unique(pred_long$iteration), 20)) %>%
  ungroup() %>%
  mutate(year = year + min(orig_data$year) - 1)
stopifnot(range(pred_long$year) == c(1972, 2019))
orig_data <- reshape2::melt(modelfit$data) %>%
  as_tibble() %>%
  rename(variable = Var1, year = Var2, data_value = value)
orig_and_pred <- inner_join(pred_long, orig_data,
  by = c("year", "variable")
) %>%
  filter(!is.na(data_value))

middle <- c(seq(1972 + 12, 1972 + 35, 3), 2007)
end <- middle + 12
begin <- middle - 12
stopifnot(min(orig_and_pred$year) == min(begin))
stopifnot(max(orig_and_pred$year) == max(end))
years_df <- data.frame(begin = begin, middle = middle, end = end)
system.time({
  out2 <- purrr::pmap_dfr(years_df, fit_cor_window)
})
saveRDS(out2, "stan-cor-window-biol.rds")
out2 <- readRDS("stan-cor-window-biol.rds")

out2 <- left_join(out2, tibble(var_name = gsub("\\.", " ", new.names), variable = 1:length(new.names)))
out2$var_name <- factor(out2$var_name, levels = gsub("\\.", " ", new.names))

out2 %>%
  group_by(var_name, middle) %>%
  summarise(
    lwr = quantile(rho_post, probs = 0.1),
    lwr2 = quantile(rho_post, probs = 0.25),
    upr = quantile(rho_post, probs = 0.9),
    upr2 = quantile(rho_post, probs = 0.75),
    med = quantile(rho_post, probs = 0.5)
  ) %>%
  ggplot(aes(x = middle, y = med, ymin = lwr, ymax = upr)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.2) +
  geom_line(colour = "#56B4E9", lwd = .7) +
  geom_ribbon(
    alpha = 0.4, fill = "#56B4E9",
    mapping = aes(ymin = lwr2, ymax = upr2)
  ) +
  geom_ribbon(alpha = 0.25, fill = "#56B4E9") +
  facet_wrap(~var_name, scales = "fixed", ncol = 4) +
  theme_bw() +
  ylab("Correlation") +
  xlab("Year (middle of 25-year window)") +
  theme(axis.title.x = element_blank()) +
  coord_cartesian(ylim = c(-0.25, 1), expand = FALSE)
ggsave("figs/stan-cor-window-biol.png", width = 5.5, height = 5.5)
