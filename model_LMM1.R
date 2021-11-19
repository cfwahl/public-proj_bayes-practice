model {
  
  # Priors
  for (j in 1:2){
    beta[j] ~ dnorm(0, 0.001)
  }
  
  tau ~ dscaled.gamma(2.5, 3)
  sigma <- sqrt(1 / tau)
  
  for (k in 1:Nspecies) {
    alpha[k] ~ dnorm(mu_alpha, tau_sp)
  }
  
  mu_alpha ~ dnorm(0, 0.001)
  tau_sp ~ dscaled.gamma(2.5, 3)
  sigma_sp <- sqrt(1 / tau_sp)
  
  # Likelihood
  for (i in 1:N) {
    Y[i] ~ dnorm(mu[i], tau)		# The 'residual' random variable
    mu[i] <- alpha[Species[i]] + beta[1]*X1[i] + beta[2]*X2[i]  # Expectation
  }
  
}
