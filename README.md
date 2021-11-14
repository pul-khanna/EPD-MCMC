# Exponential Power Distribution

This is a sample code from the main research project on Exponential Power Distribution. It is a sample demonstration on how to estimate the parameters (alpha and lambda) of an EPD using Metropolis Hastings Algorithm.

We will be performing the following operations in our code:-

- Here we are first generating a sample of size 'n' from an Exponential Power Distribution using Inverse Transform Sampling (ITS).
- Then we will generate sample observations from posterior distributions using Markov Chain Monte Carlo methodology - (Metropolis Hastings algorithm in this case).
- We will obtain Bayes Estimates under SELF-loss function. 
- For Bayesian Credible Intervals (BCI), we are computing Equal-tailed Credible Interval (ETCI) and Highest Posterior Density Credible Interval (HPDCI) and we are     also calculating their respective widths and Coverage Probabilities (CP)  
- Also, we will be calcuating Average Absolute Bias (AB) and Mean Square Error (MSE) by repeating the steps above.

Note:- We will be using  mcmc library to perform M-H algorithm and we will be including burn-in as well and maintaining an acceptance sampling rate of 20%
