
pacman::p_load(runjags,
               tidyverse)

# simulated data ----------------------------------------------------------

X <- rnorm(n = 100,
           mean = 0,
           sd = 1)

mat <- model.matrix(~X)

v_b <- c(0.1, 2.3)
Y <- rnorm(n = 100,
           mean = mat %*% v_b,
           sd = 1)


# jags --------------------------------------------------------------------

## data ####
d_jags <- list(Nsample = length(Y),
               Y = Y,
               X = X)

## parameters ####
para <- c("b",
          "sigma")

## model file ####
m <- read.jagsfile("model.R")

## mcmc setup ####

n_ad <- 100
n_iter <- 1.0E+4
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