################################################################################
###                               Project 4 Start                            ###
################################################################################

# Code for Project4

# Libraries
library(hdrm)
library(MASS) # AIC/BIC variable selection
library(olsrr) # Backwards variable selection

 
################################################################################
###                                Project Notes                             ###
################################################################################

## Correlation among predictors in simulation model: 
# - Exchangeable with settings 0, 0.35, and 0.7

##### Analysis Plan Overview
## Use the hrdm package to generate the simulation data
## Scenario 1: 
# N = 250, with 20 vars, 5 PEVs with betas in (0.5/3, 1/3, 1.5/3, 2.0/3, and 2.5/3)
# Case 1a: All X indepdendent from each other
# Case 1b: Correlations between Xs: 
#   - Exchangeable with settings 0, 0.35, and 0.7

## Scenario 2: 
# N = 500, with 20 vars, 5 PEVs with betas in (0.5/3, 1/3, 1.5/3, 2.0/3, and 2.5/3)
# Case 1a: All X indepdendent from each other
# Case 1b: Correlations between Xs: 
#   - Exchangeable with settings 0, 0.35, and 0.7

### Set up shell of functions needed for simulation

# Simulate the data using the hdrm package

################################################################################
###                                Backwards Selection                       ###
################################################################################


################################ Simulation 1 Data #############################

# Define betas
sim_betas <- c(0.5/3, 1/3, 1.5/3, 2/3, 2.5/3,
               rep(0, 15)) 

### Functions for model selection techniques


# Using Backward variable and correlation.Rmd from BIOS 6624 class as guide

# Set seed
set.seed(646)

# Specify correlation coefficient values
rho_vec <- c(0, 0.35, 0.7)

# Use function from olsrr R package
sim1 <- function(N = 10,
                 # Subjects in each sim
                 n = 250,
                 # Betas defined
                 B = sim_betas,
                 # Number of vars 1:n_vars to model
                 n_vars,
                 # Significance level
                 alpha = 0.05) {
  # Get the number of vars
  if(missing(n_vars)) {n_vars <- length(B)}
  # List of models 
  results <- list()
  # # Add names to list items
  # names(results) <- paste0("rho_", rho_vec)
  # For range of correlation values
  for(i in seq_along(rho_vec)) {
      # Simulate data
      sim_dat <- gen_data(n = n,
                          # Number of predictors
                          p = 20, 
                          # Non-zero predictors
                          p1 = 5,
                          # Vector of non-zero betas
                          beta = B,
                          # Correlation structure
                          corr = "exchangeable",
                          # Correlation coefficient
                          rho = rho_vec[i])
      # Extract outcome for model
      y <- sim_dat$y
      # Extract predictors
      X <- sim_dat$X
      # Change colnames
      colnames(X) <-  seq(1, 20, by = 1)
      # Create data frame with outcome and predictors
      model_dat <- data.frame(y, X)
      # Get all variable names
      vnames <- names(model_dat)
      # Fit regression model - exclude intercept
      results[[i]] <- lm(paste(vnames[1], '~', paste(vnames[2:n_vars],
                        collapse = '+')),
                        data = model_dat)
      # Apply backwards stepwise selection on covariates
      backwards_sel <- ols_step_backward_p(sim_model) # use include arg here?
      # use pvalue here?
  }
  
  # Summarize results
  # Create lists to store info
  bias <- list()
  
  for (i in 1:N) {
    # True positive - % of time are vars X1-X5 in the model
    # Get the coefficients
    
    # False positive rates
    
    # For variables remaining in model: 
    # Bias - Subtract estimates from defined betas / betas
    bias[[i]] <- (coef(results[[i]])[, 2:21] - B[1:n_vars])/ B[1:n_vars]
    # Coverage of the 95% CI
    
    # Type I error of vars selected
    
    # Type II error of vars selected
  }
  # Calculate values
  bias <- colMeans(do.call('rbind', bias))
  
  # Return results
  return(data.frame(bias = bias))
}



##################### Checks/Playing around
# Creating simulated data
set.seed(646)
sim1_data0 <- gen_data(n = 250,
                       # Number of predictors
                       p = 20, 
                       # Non-zero predictors
                       p1 = 5,
                       # Vector of non-zero betas
                       beta = sim_betas,
                       # Correlation structure
                       corr = "exchangeable",
                       # Correlation coefficient
                       rho = 0)

# Checking correlations between Xs and y
mean(cor(dat_list[[1]]$X)) # [1] 0.04970694
mean(cor(dat_list[[2]]$X)) # [1] 0.3970614
mean(cor(dat_list[[3]]$X)) # [1] 0.7335126
# Looks good

###################################### AIC #####################################


###################################### BIC #####################################


##################################### LASSO ###################################



################################## Elastic Net #################################




################################ Simulation 2 Data #############################



############################# Backwards Selection ##############################


###################################### AIC #####################################


###################################### BIC #####################################


##################################### LASSO ###################################



################################## Elastic Net #################################












