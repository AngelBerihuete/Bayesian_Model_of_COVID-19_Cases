library(tidyverse)
library(tsibble)
library(lubridate)
library(fabletools)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# Reading the data from 
# https://cnecovid.isciii.es/covid19/#documentaci%C3%B3n-y-datos
# --------------------------------------------------------------
#

data <- read_csv("casos_diagnostico_provincia.csv")

# To show the data
# ----------------
#
# 1. Filtering the data for the wonderful Cadiz province.
# 2. Create a new variable with all those who have a positive test on Polymerase
#    Chain Reaction (PCR) plus those positive in a rapid antibody test made in 
#    laboratory.
# 3. Coerce to a tsibble object.
# 4. Plot cases time series.

data %>%
  filter(provincia_iso=="CA") %>%
  mutate(cases = num_casos_prueba_pcr +  num_casos_prueba_test_ac) %>%
  as_tsibble() %>%
  autoplot(cases)

# Lets analyze the first wave
# ---------------------------
#
# 1. Select the dates

init_date <- ymd("2020-03-02")
end_date <- init_date + weeks(9)

# Showing the first wave
data %>%
  filter(init_date <= fecha & fecha <= end_date ) %>%
  filter(provincia_iso=="CA") %>%
  mutate(cases = num_casos_prueba_pcr +  num_casos_prueba_test_ac) %>%
  as_tsibble() %>%
  autoplot(cases)

# 2. Select the data and change the variable name "fecha" (spanish) by "dates"

data %>%
  rename(dates = fecha) %>%
  filter(init_date <= dates & dates <= end_date ) %>%
  filter(provincia_iso=="CA") %>%
  mutate(cases = num_casos_prueba_pcr +  num_casos_prueba_test_ac) %>%
  as_tsibble() -> data_first_wave


# 3. Setting up the data to feed the Bayesian model

P <- 1240155/100000 # 1240155 inhabitants in the province of Cadiz

t <- rep(data_first_wave$dates, data_first_wave$cases)
t <- sort(as.numeric(t))
y <- t - t[1]
y <- y[-1]

capital_t <- y[length(y)]
N <- length(y)

covid19_data <- list(N = N, P = P, y = y, capital_t = capital_t)


# 5. Initializing the chains

initf <- function(chain_id = 1) {
  list(
    a = runif(1, 5000, 20000),
    M0 = chain_id+1,
    c = runif(1, 0.01, 0.2)
  )
}

n_chains <- 4
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

# 6. Fitting the model

fit_model <- stan(
  file = 'COVID-19_cases.stan',
  data = covid19_data,
  chains = n_chains,
  init = init_ll,
  iter = 5000,
  control = list(adapt_delta = 0.95),
  verbose = FALSE,
  open_progress = FALSE)

# Checking the inference and plotting posterior distributions
# -----------------------------------------------------------
#
# 1. Summary of the inference

fit_model

# 2. Parameters pairs plot

pairs(fit_model, pars = c("a", "b", "c", "lp__"))

# Checkout shinystan package for a more detailed and fancy summaries.


# Finally some equations in the paper
# -----------------------------------
#
# The expected number of new cases of COVID-19 in future time intervals of the
# form (T+h1,T+h2), h1 < h2

ET <- function (a, b, c, h1, h2, capital_t) { # equation 9 in the paper
  a*(exp(-b*exp(-c*(capital_t+h2)))-exp(-b*exp(-c*(capital_t+h1))))
}

# The expected number of new cases of COVID-19 in future days from T (capital_t)
# i.e., h1=0

DT <- function (a, b, c, h2, capital_t) {
  a*(exp(-b*exp(-c*(capital_t+h2)))-exp(-b*exp(-c*(capital_t+(h2-1)))))
}

# A band of Gompertz curves given by an i.i.d. random sample of size $500$ of 
# $g(t \mid \pi_{\bf x})$ as in Fig. 3 in the paper.

gompertz <- function(a, b, c, initdate = init_date, enddate=end_date){
  t <- 0:(as.numeric(enddate - (initdate + days(1))))
  return (a*exp(-b*exp(-c*t)))
}

parameters <- tibble(as.data.frame(extract(
  fit_model, pars=c("a", "b", "c"), permuted = TRUE)))

parameters %>% rowwise() %>%
  transmute(gptz_f = list(gompertz(a, b, c)))  %>% 
  flatten_dfc() %>% 
  select(sample((1:dim(parameters)[1]), 500)) %>%
  as_tibble() %>%
  mutate(time = init_date + days(1:(as.numeric(end_date - init_date)))) %>%
  gather(key = "curve", value = "value", -time) -> df


ggplot(data = df, aes(time, value)) + geom_line(aes(group=curve), color="grey")
