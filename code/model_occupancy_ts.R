
model {
  
  # Priors and linear models
  alpha <- logit(mean.theta)                      # Intercept on logit link scale
  mean.theta ~ dunif(0, 1)                        # Intercept on prob. scale
  beta1 ~ dnorm(0, 0.001)                        # Slope on logit link scale
  beta2 ~ dnorm(0, 0.001) 
  
  # Likelihood of the Bernoulli GLM
  for (i in 1:n){
    y[i] ~ dbern(theta[i])                        # Stochastic part of response
    logit(theta[i]) <- alpha + beta1 * agr[i] + beta2 * elv[i]    # Link function and lin. pred.
    }
}