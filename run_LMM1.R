
# Load packages
library(runjags)
library(tidyverse)
library(MCMCvis) 

# jags --------------------------------------------------------------------

## data ####
# standardized covariate length in this example
d_jags <- list(Y = iris$Sepal.Length,
               X1 = iris$Sepal.Width,
               X2 = iris$Petal.Length,
               Species = iris$Species,
               N = nrow(iris),
               Nspecies = n_distinct(iris$Species))

## parameters ####
para <- c("mu_alpha",
          "alpha",
          "beta",
          "sigma")

## model file ####
m <- read.jagsfile("model_LMM1.R")

## mcmc setup ####
n_ad <- 100 
n_iter <- 1.0E+4 #number of draws
n_thin <- max(3, ceiling(n_iter / 500)) #number of thins
n_burn <- ceiling(max(10, n_iter/2)) # number of draws to burn
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

mcmc_summary <- MCMCsummary(post$mcmc)  
mcmc_summary
