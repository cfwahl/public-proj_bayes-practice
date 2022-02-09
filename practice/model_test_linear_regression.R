model {
  
  # Priors
  for(j in 1:3) {
    b[j] ~ dnorm(0, 0.001)
  }	

  sigma ~ dunif(0, 100)			# Residual standard deviation
  tau <- 1 / ( sigma * sigma)
  
 # Likelihood
    for (i in 1:n) {
      y[i] ~ dnorm(mu[i], tau) 
      mu[i] <- b[1] + b[2] * x1[i] + b[3] *x2[i]
    }
  }
