model {
  # variables defined in R are capitalized
  # lowercase variables are parameters or defined in this JAGS model
  
  # prior
  b[1] ~ dnorm(0, 0.01)
  b[2] ~ dnorm(0, 0.01)
  
  tau ~ dscaled.gamma(2.5, 1)
  sigma <- sqrt(1 / tau)
  
  # likelihood
  for(i in 1:Nsample) {
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- b[1] + b[2] * X[i]
  }
  
}