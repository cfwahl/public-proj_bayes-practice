#########################################################################
#
# Bayesian Population Analysis with JAGS
#
# Based on the book "Bayesian population analysis using WinBUGS - a hierarchical perspective" by Marc Kéry & Michael Schaub (2012, Academic Press)
#
#########################################################################

# Michael Schaub, 24.1.2012, revised 9.7.2012, revised 10.12.2013

# This is the BPA code "translated" into JAGS. In many cases there are no differences at all compared to the BUGS code, but sometimes  small changes either in the model code or in the function to provide initial values were necessary. When such changes occur, they are most often commented in the code.

#########################################################################

# Load packages

library(lattice)
library(coda)
library(R2WinBUGS)
library(runjags)
library(tidyverse)
#library(R2jags)

# make sure JAGS is talking to R
testjags()

# setwd("C:/....")     # Optional


##########################################################################
# 
# 1. Introduction
# 
##########################################################################

# 1.3. The binomial distribution as a canonical description of the observation process
N <- 16                 # Population size of sparrows in the yard
p <- 0.4                # Individual detection probability
rbinom(n = 1, size = N, prob = p)
rbinom(n = 1, size = N, prob = p)
rbinom(n = 1, size = N, prob = p)
C <- rbinom(n = 1000000, size = N, prob = p)
mean(C)
var(C)
sd(C)
hist(C, breaks = 50, col = "gray", main = "", xlab = "Sparrow count", 
     las = 1, freq = FALSE)


# 1.4. Structure and overview of the contents of this book
# Population values for mean and standard deviation of individual lengths
mu <- 65            # Population mean
sigma <- 5          # Population SD

# Draw a single sample of 10 and summarize it
x <- rnorm(n = 10, mean = mu, sd = sigma)
reps <- 10^6
sample.means <- rep(NA, reps)
for (i in 1:reps){
  sample.means[i] <- mean(rnorm(n = 10, mean = mu, sd = sigma))
}

# Produce figure
par(mfrow = c(1, 2), las = 1)
hist(x, col = "gray", main = "", xlab = "Body length (cm)", las = 1)
abline(v = mu, lwd = 3, col = "red")
abline(v = mean(x), lwd = 3, col = "blue")
hist(sample.means, col = "gray", main = "", xlab = "Body length (cm)", nclass = 50, freq = FALSE, las = 1)
abline(v = mu, lwd = 5, col = "red")
abline(v = mean(sample.means), lwd = 5, col = "blue", lty = 2)
#SE of the mean, also SD of the distribution 
sd(sample.means)

#########################################################################
# 
# 3. Introduction to the generalized linear model (GLM): The simplest model for count data
#
##########################################################################

# 3.2. Statistical models: Response = Signal + Noise
# 3.2.1. The noise component
plot(density(rbeta(n = 10000, shape1 = 2, shape2 = 4)))
hist(rbeta(10^6, 2, 4), nclass = 100, col = "gray")
     
# 3.2.2. The signal component
# Define and plot data
y <- c(25, 14, 68, 79, 64, 139, 49, 119, 111)
A <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3))
X <- c(1, 14, 22, 2, 9, 20, 2, 13, 22)
plot(X, y, col = c(rep("red", 3), rep("blue", 3), rep("green", 3)), xlim = c(-1, 25), ylim = c(0, 140))

# fit an ANCOVA
summary(fm <- lm(y ~ A-1 + X))

# effect or treatment contrast parameterization
# way to examine design matrix of the model
model.matrix(~ A + X)
# means parameterization 
model.matrix(~ A-1 + X)

####### 3.3. Poisson GLM in R and WinBUGS for modeling times series of counts

# 3.3.1. Generation and analysis of simulated data
data.fn <- function(n = 40, alpha = 3.5576, beta1 = -0.0912, 
                    beta2 = 0.0091, beta3 = -0.00014){
# n: Number of years
# alpha, beta1, beta2, beta3: coefficients of a 
#    cubic polynomial of count on year
  
# Generate values of time covariate
year <- 1:n
  
# Signal: Build up systematic part of the GLM
log.expected.count <- alpha + beta1 * year + beta2 * year^2 + beta3 * year^3
expected.count <- exp(log.expected.count)
  
  # Noise: generate random part of the GLM: Poisson noise around expected counts
  C <- rpois(n = n, lambda = expected.count)
  
  # Plot simulated data
  plot(year, C, type = "b", lwd = 2, col = "black", main = "", las = 1, 
       ylab = "Population size", xlab = "Year", cex.lab = 1.2, cex.axis = 1.2)
  lines(year, expected.count, type = "l", lwd = 3, col = "red")
  
  return(list(n = n, alpha = alpha, beta1 = beta1, beta2 = beta2, beta3 = beta3, 
              year = year, expected.count = expected.count, C = C))
}

# simulated population size of falcons over 40 years
data <- data.fn()

# fitting GLM frequentist 
fm <- glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data)
summary(fm)

# Specify model in BUGS language
sink("GLM_Poisson.jags")
cat("
model {

# Priors
alpha ~ dunif(-20, 20)
beta1 ~ dunif(-10, 10)
beta2 ~ dunif(-10, 10)
beta3 ~ dunif(-10, 10)

# Likelihood: Note key components of a GLM on one line each
for (i in 1:n){
   C[i] ~ dpois(lambda[i])          # 1. Distribution for random part
   log(lambda[i]) <- log.lambda[i]  # 2. Link function
   log.lambda[i] <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) + 
   beta3 * pow(year[i],3)          # 3. Linear predictor
   } #i
}
",fill = TRUE)
sink()

bugs.dir <- "/ProgramFiles/JAGS/JAGS-4.3.0/i386/bin/jags-terminal.exe"

# Bundle data
win.data <- list(C = data$C, n = length(data$C), year = data$year)

# Initial values
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda")

# MCMC settings
ni <- 2000 # number of draws
nt <- 2 # thinning rate
nb <- 1000 # burnin length
nc <- 3 # number of chains


# Call JAGS from R
out <- run.jags(data = win.data, inits = inits, parameters.to.save = params, 
            model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, 
            n.iter = ni, n.burnin = nb, working.directory = getwd())

# Bundle data
# Standardize year covariate
mean.year <- mean(data$year)             # Mean of year covariate
sd.year <- sd(data$year)                 # SD of year covariate
win.data <- list(C = data$C, n = length(data$C), 
                 year = (data$year - mean.year) / sd.year)

## model file ####
m <- read.jagsfile("model")

# Call JAGS from R (BRT < 1 min)
out <- run.jags(data = win.data, inits = inits, parameters.to.save = params, 
            model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, 
            n.iter = ni, n.burnin = nb, working.directory = getwd())

# Call JAGS from R (BRT < 1 min)
out <- run.jags(m$model, monitor = params, data = win.data, inits = inits,  
                n.chains = nc, thin = nt, method = "parallel",
                sample = ni, burnin = nb, n.sims=3, module = "glm")

print(out, dig = 3)

# New MCMC settings with essentially no burnin
ni <- 100
nt <- 1
nb <- 1

# Call JAGS from R (BRT < 1 min)
tmp <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

plot(1:40, data$C, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Population size", xlab = "Year")
R.predictions <- predict(glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data), type = "response")
lines(1:40, R.predictions, type = "l", lwd = 3, col = "green")
JAGS.predictions <- out$BUGSoutput$mean$lambda
lines(1:40, JAGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)
cbind(R.predictions, JAGS.predictions)

# An option in JAGS to see the traceplots of all parameters is the following:
traceplot(tmp)




