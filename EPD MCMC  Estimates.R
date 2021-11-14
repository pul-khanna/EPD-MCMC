# Clearing Old Variables from Environment
rm(list = ls(all = TRUE))

library(mcmc)

BayesAl <- function(N, n, alpha, lambda) {
  
  # Adjusting Starting values for MCMC to shorten converging times 
  alpha_ins <- alpha - 0.1
  lambda_ins <- lambda - 0.1
  burnIn <- 5000
  
  # We chose parameters of prior distribution large in such a way the mean of the prior distribution
  # is equal to the parameter of the EPD
  b1 <- alpha / 0.01
  a1 <- alpha * b1
  b2 <- lambda / 0.01
  a2 <- lambda * b2
  
  # Mean value 
  Mean_A <- rep(0, N)
  Mean_L <- rep(0, N)
  
  # Standard Errors
  SE_A <- rep(0, N)
  SE_L <- rep(0, N)
  
  # Equal Density and the Coverage Probability
  EPD_A <- array(dim = c(N, 3))
  EPD_L <- array(dim = c(N, 3))
  LED_A <- rep(0, N)
  RED_A <- rep(0, N)
  CPED_A <- rep(0, N)
  LED_L <- rep(0, N)
  RED_L <- rep(0, N)
  CPED_L <- rep(0, N)
  
  # Highest Posterior Density and the Coverage Probability
  HPD_A <- array(dim = c(N, 3))
  HPD_L <- array(dim = c(N, 3))
  LHD_A <- rep(0, N)
  RHD_A <- rep(0, N)
  CPHD_A <- rep(0, N)
  LHD_L <- rep(0, N)
  RHD_L <- rep(0, N)
  CPHD_L <- rep(0, N)
  
  for(i in 1:N) {
    
    # Generate EPD values Using IPT
    U <- runif(n, min = 0, max = 1)
    X <- (1 / lambda) * (log(1 - log(1 - U))) ^ (1 / alpha)
    
    # Function for mcmc Package
    lupost <- function(param) {
    
      stopifnot(is.numeric(param))
      stopifnot(is.finite(param))
      stopifnot(length(param) == 2)
      
      # Absolute values to prevent Nan as log of negative numbers not possible    
      a <- abs(param[1])
      l <- abs(param[2])
      
      # Log-Likelihood 
      logL <- sum ( log(a) +
                    (a * log(l)) +
                    ((a - 1) * log(X)) +
                    ((l * X) ^ a) +
                    (1) - 
                    (exp((l * X) ^ a)) )
      
      # Prior Distribution - Gamma
      aprior <- (a1 * log(b1)) +
        (a * (-b1)) +
        ((a1 - 1) * log(a)) -
        lgamma(a1)
      lprior <- (a2 * log(b2)) +
        (l * (-b2)) +
        ((a2 - 1) * log(l)) -
        (lgamma(a2))
      # Posterior Distribution
      return (logL + aprior + lprior)
    }
    
    # Metropolis Algorithm 
    out = metrop(lupost, initial = c(alpha_ins, lambda_ins), nbatch = 5000, blen = 1)
    
    print(out$accept)
    #set scale so that acceptance rate near about 20%
    out <- metrop(out, scale = 0.21)
    print(out$accept)
    out <- metrop(out, scale = 0.095)
    out$accept
    
    out <- metrop(out, nbatch = 20000, blen = 1)
    
    # Removing BurnIn Observations
    phi <- abs(out$batch[-(1:burnIn), ])
    phi_sort <- apply(phi, 2, sort)

    # Mean Estimates
    Mean_A[i] <- mean(phi[, 1])
    Mean_L[i] <- mean(phi[, 2])
    
    # Standard Error Estimates
    SE <- apply(phi, 2, sd)
    SE_A[i] <- SE[1]
    SE_L[i] <- SE[2]
    
    # Equal-Tailed Credible Interval
    ETCI <- function(phi_sort) {
      CI <- (quantile(phi_sort, probs = c(0.025, 0.975)))
      CI_width = CI[2] - CI[1]
      CPET <- (CI[1] < alpha & CI[2] > alpha)
      return (c(CI[1], CI[2], CPET))
    }
    
    EPD_A[i,] <- ETCI(phi_sort[,1])
    EPD_L[i,] <- ETCI(phi_sort[,2])
    
    # Highest Posterior Credible Interval
    HPDCI <- function(phi_sort , credMass = 0.95) {
      ciIdxInc <- floor(credMass * length(phi_sort))
      nCIs <- length(phi_sort) - ciIdxInc
      ci_width <- rep(0 , nCIs)
      for (i in 1:nCIs) {
        ci_width[i] <- phi_sort[i + ciIdxInc] - phi_sort[i]
      }
      HPDmin <- phi_sort[which.min(ci_width)]
      HPDmax <- phi_sort[which.min(ci_width) + ciIdxInc]
      CPHPD <- (HPDmin < alpha & HPDmax > alpha)
      return (c(HPDmin, HPDmax, CPHPD))
    }
    
    HPD_A[i,] <- HPDCI(phi_sort[,1])
    HPD_L[i,] <- HPDCI(phi_sort[,2])
  
  }
  
  Mean_A1 <- mean(Mean_A, na.rm = T)
  Mean_L1 <- mean(Mean_L, na.rm = T)
  
  Bias_A1 <- mean(abs(Mean_A[i] - alpha))
  Bias_L1 <- mean(abs(Mean_L[i] - lambda))
  
  MSE_A1 <- mean((Mean_A[i] - alpha) ^ 2)
  MSE_L1 <- mean((Mean_L[i] - lambda) ^ 2)
  
  SE_A1 <- mean(SE_A, na.rm = T)
  SE_L1 <- mean(SE_L, na.rm = T)
  
  LED_A1 <- mean(EPD_A[, 1], na.rm = T)
  RED_A1 <- mean(EPD_A[, 2], na.rm = T)
  WED_A1 <- RED_A1 - LED_A1
  CPED_A1 <- mean(EPD_A[, 3], na.rm = T)
  LED_L1 <- mean(EPD_L[, 1], na.rm = T)
  RED_L1 <- mean(EPD_L[, 2], na.rm = T)
  WED_L1 <- RED_L1 - LED_L1
  CPED_L1 <- mean(EPD_L[, 3], na.rm = T)
  
  LHD_A1 <- mean(HPD_A[, 1], na.rm = T)
  RHD_A1 <- mean(HPD_A[, 2], na.rm = T)
  WHD_A1 <- RHD_A1 - LHD_A1
  CPHD_A1 <- mean(HPD_A[, 3], na.rm = T)
  LHD_L1 <- mean(HPD_L[, 1], na.rm = T)
  RHD_L1 <- mean(HPD_L[, 2], na.rm = T)
  WHD_L1 <- RHD_L1 - LHD_L1
  CPHD_L1 <- mean(HPD_L[, 3], na.rm = T)
  
  A <- data.frame(n, alpha, Mean_A1, Bias_A1, MSE_A1, SE_A1, LED_A1, RED_A1, WED_A1, CPED_A1,
                  LHD_A1, RHD_A1, WHD_A1, CPHD_A1,
                lambda, Mean_L1, Bias_L1, MSE_L1, SE_L1, LED_L1, RED_L1, WED_L1, CPED_L1, 
                LHD_L1, RHD_L1, WHD_L1, CPHD_L1)
  return (A)
    
}
  
R1 <- BayesAl(n = 20, N = 100, alpha = 0.5, lambda = 0.2)
R2 <- BayesAl(n = 50, N = 100, alpha = 0.5, lambda = 0.2)
R3 <- BayesAl(n = 20, N = 100, alpha = 1.5, lambda = 3)
R4 <- BayesAl(n = 50, N = 100, alpha = 1.5, lambda = 3)

DD1 <- rbind(R1, R2, R3, R4)
DD1           

write.csv(x = DD1, file = "Exponential Power Distribution MCMC Estimates.csv") 
