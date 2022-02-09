
# Load packages
library(runjags)
library(tidyverse)
library(MCMCvis) 

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

n <- nrow(iris)

# jags --------------------------------------------------------------------

## data ####
# 
d_jags <- list(y = iris$Sepal.Length,
               x1 = x1,
               x2 = x2,
               x3 = as.numeric(x3), 
               n = n,
               W = diag(3))


## parameters ####
para <- c("beta",
          "sigma",
          "Sigma.beta")


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

# convert to mcmcOutput and do nice plots of marginal probabilities
salB <- mcmcOutput(post, default='psi')
plot(salB)
