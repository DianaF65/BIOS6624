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

### Functions for model selection techniques

############################### Backwards Selection ############################

# Using Backward variable and correlation.Rmd from BIOS 6624 class as guide

# Set seed
set.seed(646)

# Specify correlation coefficient values
rho_vec <- c(0, 0.35, 0.7)

# List of data frames for different rhos
dat_list <- list()

# Use function from olsrr R package
sim1 <- function(n = 10000) {
  # For range of correlation values
  for (i in 1:3) {
    for(j in rho_vec) {
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
                                rho = j)
    }
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
  
  # 
}

#### Checks
# Checking correlations between Xs and y
a <- data.frame(dat_list[[1]]$y, dat_list[[1]]$X)
cor(a)
b <- data.frame(dat_list[[2]]$y, dat_list[[2]]$X)
cor(b)


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












