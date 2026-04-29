################################################################################
###                               Project 4 Start                            ###
################################################################################

# Code for Project4

# Libraries
library(hdrm)
library(MASS) # AIC/BIC variable selection
library(olsrr) # Backwards variable selection

 
##### Notes
## Correlation among predictors in simulation model: 
# - Exchaneable with settinggs 0, 0.35, and 0.7

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

################################ Simulation 1 Data #############################

# Define betas
sim_betas <- c(0.5/3, 1/3, 1.5/3, 2/3, 2.5/3,
               rep(0, 15))

### Functions for model selection techniques

############################### Backwards Selection ############################

# Using Backward variable and correlation.Rmd from BIOS 6624 class as guide

# Set seed
set.seed(646)

# Specify correlation coefficient values
rho_vec <- c(0, 0.35, 0.7)

# Use function from olsrr R package
sim1 <- function(n = 10000) {
  # Create list of data frame for different rhos
  dat_list <- vector("list", length(rho_vec))
  # Add names to list items
  names(dat_list) <- paste0("rho_", rho_vec)
  
  # For range of correlation values
  for(i in seq_along(rho_vec)) {
    
      # Simulate data
      dat_list[[i]] <- gen_data(n = 250,
                                # Number of predictors
                                p = 20, 
                                # Non-zero predictors
                                p1 = 5,
                                # Vector of non-zero betas
                                beta = sim_betas,
                                # Correlation structure
                                corr = "exchangeable",
                                # Correlation coefficient
                                rho = rho_vec[i])
    }
  # Extract outcome
  y <- sim1_data0$y # 
  
  # Extract predictors
  X <- sim1_data0$X
  # Change colnames
  colnames(X) <-  seq(1, 20, by = 1)
  
  # Create data frame with outcome and predictors
  model_dat <- data.frame(y, X)
  
  # Fit regression model
  sim_model <- lm(y ~ ., data = model_dat)
  
  # Apply backwards stepwise selection on covariates
  backwards_sel <- ols_step_backward_p(sim_model,
                        include = c("X1", "X2", "X3", "X4", "X5"))
  
  ### COME Back here
  # Summarize results
  
  # True positive
  
  # False positive rates
  
  # For variables remaining in model: 
  # Bias
  
  # Coverage of the 95% CI
  
  # Type I error of vars selected
  
  # Type II error of vars selected
}

#### Checks/Playing around
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












