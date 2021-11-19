model {
  
  # Priors
  for (j in 1:3){
    beta[j, 1:3] ~ dmnorm(beta.hat[], Tau.beta[,])
  }
  
  for(i in 1:3) {
    beta.hat[i] ~ dnorm(0, 0.001)
  }
  
  Tau.beta[1:3,1:3] ~ dwish(W[,], 4)
  Sigma.beta[1:3, 1:3] <- inverse(Tau.beta[,])
  
  tau ~ dscaled.gamma(2.5, 3)
  sigma <- sqrt(1 / tau)
  
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dnorm(mu[i], tau)		# The 'residual' random variable
    mu[i] <- beta[x3[i], 1] + beta[x3[i], 2]*x1[i] + beta[x3[i], 3]*x2[i]  # Expectation
  }
  
}
