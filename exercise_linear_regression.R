
########  MULTIPULE LINEAR REGRESSION  ###############


# Load packages
library(lattice)
library(runjags)
library(tidyverse)
library(MCMCvis) 

# dataset
attach(iris)

y <- iris$Sepal.Length 
x1 <- iris$Sepal.Width 
x2 <- iris$Petal.Length 

# response: Y would be Sepal.Length
# explanatory: X would be Sepal.Width and Petal.Length

n <- length(Sepal.Length)  # Number of data points in iris


# jags --------------------------------------------------------------------

## data ####
# standardized covariate length in this example
d_jags <- list(y = y,
               x1 = x1,
               x2 = x2,
               n = n)

## parameters ####
para <- c("b",
          "sigma")

## model file ####
m <- read.jagsfile("model.test.linear.regression.R")

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

post
plot(post)

mcmc_summary <- MCMCsummary(post$mcmc)  
mcmc_summary

library(mcmcOutput)
mcmcOutput::diagPlot(post)
