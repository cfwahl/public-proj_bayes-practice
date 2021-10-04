model {
  
  # Priors
  alpha ~ dnorm(0, 0.01)
  beta1 ~ dnorm(0, 0.01)
  beta2 ~ dnorm(0, 0.01)
  beta3 ~ dnorm(0, 0.01)
  
  # Likelihood: Note key components of a GLM on one line each
  for (i in 1:N){
    C[i] ~ dpois(lambda[i])          # 1. Distribution for random part
    log(lambda[i]) <- log.lambda[i]  # 2. Link function
    log.lambda[i] <- alpha + beta1 * Year[i] + beta2 * pow(Year[i],2) + beta3 * pow(Year[i],3)          # 3. Linear predictor
  } #i
}
