
########  Linear Mixed Model  ###############

# Load packages
library(lattice)
library(runjags)
library(tidyverse)
library(MCMCvis) 

# make sure JAGS is talking to R
testjags()

setwd("C:/....")

# dataset
attach(iris)

# assign variables
y <- iris$Sepal.Length 
x1 <- iris$Sepal.Width 
x2 <- iris$Petal.Length 
x3 <- iris$Species

# response: Y would be Sepal.Length
# explanatory: X would be Sepal.Width and Petal.Length, species random variable
# use half student_t distribution for priors of sds

n.species <- 3				# Number of species
n.sample <- 50				# Number of samples in each pop
n <- n.species * n.sample 		# Total number of data points


# jags --------------------------------------------------------------------

## data ####
# 
d_jags <- list(y = y,
               x3 = as.numeric(x3), 
               x1 = x1,
               x2 = x2,
               n = n)


## parameters ####
para <- c("alpha",
          "beta",
          "mu.int",
          "sigma.int",
          "mu.slope", 
          "sigma.slope",
          "rho",
          "covariance",
          "sigma")


## model file ####
m <- read.jagsfile("model.LMM2.R")

## mcmc setup ####
n_ad <- 100 
n_iter <- 1.0E+4 #number of draws
n_thin <- max(3, ceiling(n_iter / 500)) #number of thins
n_burn <- ceiling(max(10, n_iter/2)) # number of draws to burn
n_sample <- ceiling(n_iter / n_thin)

# These were changed for this example. Added alpha, beta, sigma to Akiras code
inits <- replicate(3,
                   list(.RNG.name = "base::Mersenne-Twister",
                        .RNG.seed = NA, 
                        mu.int = rnorm(1, 0, 1),
                        sigma.int = rlnorm(1), 
                        mu.slope = rnorm(1, 0, 1),
                        sigma.slope = rlnorm(1),
                        rho = runif(1, -1, 1),
                        sigma = rlnorm(1)),
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

# convert to mcmcOutput and do nice plots of marginal probabilities
salB <- mcmcOutput(post, default='psi')
plot(salB)
