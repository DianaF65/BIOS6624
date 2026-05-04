################################################################################
###                               Project 4 Start                            ###
################################################################################

# Code for Project4

# Libraries
library(hdrm)
library(MASS) # AIC/BIC variable selection
library(olsrr) # Backwards variable selection
library(glmnet) # LASSO

 
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
                 n_vars = 20,
                 # Significance level
                 alpha = 0.05) {
  # Get the number of vars
  if(missing(n_vars)) {n_vars <- length(B)}
  # List of model estimates for fit with lm
  model_results <- list()
  # List of selection results
  backward_sel <- list()
  # Define counter 
  counter <- 1
  # Simulate data N times
  for (i in 1:N) {
    # For range of correlation values
    for (j in seq_along(rho_vec)) {
      # Simulate data
      sim_dat <- gen_data(n = n,
                          # Number of predictors
                          p = n_vars, 
                          # Non-zero predictors
                          p1 = 5,
                          # Vector of non-zero betas
                          beta = B,
                          # Correlation structure
                          corr = "exchangeable",
                          # Correlation coefficient
                          rho = rho_vec[j])
      # Extract outcome for model
      y <- sim_dat$y
      # Extract predictors
      X <- sim_dat$X
      # Change colnames
      colnames(X) <-  paste0("X", 1:n_vars)
      # Create data frame with outcome and predictors
      model_dat <- data.frame(y = y, X)
      # Get all variable names
      vnames <- names(model_dat)
      # Fit regression model - exclude intercept
      model_results[[i]] <- lm(as.formula(paste(vnames[1], '~', 
                                  paste(vnames[2:(n_vars + 1)],
                                  collapse = '+'))),
                         data = model_dat)
      ### Backwards selection
      backwards_sel <- ols_step_backward_p(model_results[[i]]) # use include arg here?
      # Extract model
      backwards_fit <- backwards_sel$model
      # Get the selected variables and remove intercept
      selected_vars <- names(coef(backwards_fit))[-1]
      # Get the vars with defined betas
      true_vars <- paste0("X", 1:5)
      # Get the remaining vars
      null_vars <- paste0("X", 6:n_vars)
      
      ### AIC selection
      aic_sel <- step(model_results[[i]],
                      trace = 0,
                      direction = "backward",
                      k = 2) # AIC
      
      # Summarize results
      # Calculate true positives
      tp <- sum(true_vars %in% selected_vars)
      # Calculate false positives
      fp <- sum(null_vars %in% selected_vars)
      
      # Add tp, fp to resuts list
      # results[[counter]] <- final_fit
      
      # Create data frame containing results
      # results[[counter]] <- data.frame(
      #   sim = i,
      #   rho = rho_vec[j],
      #   true_positives = tp,
      #   false_positives = fp)
      # 
      # # Increase counter by 1
      # counter <- counter + 1
    }
    
  }
  
  # Summarize results
  # Create lists to store info
  bias <- list(); coverage <- list()
  
  for (i in 1:N) {
    # True positive - % of time are vars X1-X5 in the model
    # Get the true positives
    # selection_results <- do.call(rbind, results)
    
    # aggregate(true_positives ~ rho, data = results, mean)
    # aggregate(false_positives ~ rho, data = results, mean)
     
    # False positive rates
    
    # For variables remaining in model: 
    ### Backwards
    # Bias - Raw bias and not relative bias
    bias[[i]] <- coef(model_results[[i]])[-1] - B  # Remove the intercept term
    # Coverage of the 95% CI
    ci <- confint.default(model_results[[i]], level = 1 - alpha)
    coverage[[i]] <-  ci[, 1][-1] <= B & B <= ci[, 2][-1]
    
    # Type I error of vars selected
    
    # Type II error of vars selected
    
  }
  # Calculate colmeans 
  bias <- colMeans(do.call('rbind', bias))
  coverage <- colMeans(do.call('rbind', coverage))
  
  # Create a data frame of model performance
  model_performance_df <- data.frame(bias = bias, 
                                     coverage = coverage)
  
  # Return results
  return(model_performance_df)
}

# Run simulation
sim_results <- sim1()

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


# Fit model for lasso
glmnet_model <- glmnet(y , x, 
                       family = c("gaussian"),
                       alpha = 1)


################################## Elastic Net #################################












