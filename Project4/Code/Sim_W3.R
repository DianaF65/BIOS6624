####### Simulation Worksheet 2

# Load libraries
library(hdrm)


# 1)	Research: find R packages that can perform variable selection by AIC and BIC.
# Step from R and stepAIC


# 2)	Use the code from the lecture to explore bias caused by the exclusion of a 
# variable (x2) correlated with the outcome (y) and a focus variable (x1). 
# a.	What simulation setting can stay fixed?
# B, mean of cov, and evar

# Original settings: 
# res
# 1      2 
# 0.4121 0.5879 
# B, cov my, and evar can change - basically anything besides cov var.

# b.	What setting needs to vary?
#  The cov var matrix
#   
# 3)	Build a simple function to produce a simulated data set that replicates a 
# collider effect. Confirm your theories about the effect on the exposure when the
# collider is included and excluded from the model.
# 
# Includig the collider in the model would show a false association that is not
# truly from the explanatory variable of interest
# 
# Adapting sim dat from Backward variable selecion rmd file from class
sim_dat2 <- function( n = 189 # size of each simulation dataset
                     , B = c(-1.0871, 0.7041, 0)# true model parameters
                     , cov_mu = c(0, 0)  # means of x variables
                     , cov_var = matrix(c(1, 0.8, 0,)) # covariance matrix of X variables
                     , e_var = 1
) {
  
  # add generation of multiple, possibly correlated variables using the multivariate normal distribution 
  x <- MASS::mvrnorm(n=n, mu = cov_mu, Sigma = cov_var)
  
  x_matrix <- cbind(intercept = rep(1, n), x)
  
  v = as.vector(x_matrix %*% B)
  y = sapply(v, function(x){rnorm(1, mean = x, sd = sqrt(e_var))})
  
  simdat <- as.data.frame(cbind(y, x))
  vnames <- c('y',paste0('x',1:length(cov_mu)))
  names(simdat) <- vnames
  
  return(simdat2)
}
dat2 = sim_dat2()


