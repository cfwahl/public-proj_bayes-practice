
pacman::p_load(runjags,
               tidyverse)

# simulated data ----------------------------------------------------------

source("simulation_data.R")

data <- as_tibble(data)
write_csv(data, "sample_data.csv")


# jags --------------------------------------------------------------------

data <- read_csv("sample_data.csv")

## data ####
d_jags <- list(C = data$C,
               Year = data$year,
               N = data$n)

## parameters ####
para <- c("alpha",
          "beta1",
          "beta2",
          "beta3",
          "lambda")

## model file ####
m <- read.jagsfile("model_poisson.R")

## mcmc setup ####
n_ad <- 100
n_iter <- 2.0E+4
n_thin <- max(3, ceiling(n_iter / 500))
n_burn <- ceiling(max(10, n_iter/2))
n_sample <- ceiling(n_iter / n_thin)

inits <- replicate(3,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA),
                   simplify = FALSE)

for (i in 1:3) inits[[i]]$.RNG.seed <- i


# run jags ----------------------------------------------------------------

post <- run.jags(m$model,
                 monitor = para,
                 data = d_jags,
                 n.chains = 3,
                 inits = inits,
                 method = "parallel",
                 burnin = n_burn,
                 sample = n_sample,
                 adapt = n_ad,
                 thin = n_thin,
                 n.sims = 3,
                 module = "glm")

mcmc_summary <- MCMCvis::MCMCsummary(post$mcmc)