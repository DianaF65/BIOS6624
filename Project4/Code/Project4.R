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
  
  # Lists to store model results and rho values
  model_results <- list(); rho_results <- c()
  
  # Define counter to give each model its own slot because i will not do that
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
      # Change colnames of data frame
      colnames(model_dat) <- c("y", paste0("X", 1:n_vars))
      
      # Get all variable names
      vnames <- names(model_dat)
      
      # Fit regression model - exclude intercept
      full_fit <- lm(as.formula(paste(vnames[1], '~', 
                                  paste(vnames[2:(n_vars + 1)],
                                  collapse = '+'))),
                         data = model_dat)
      
      ### Backwards selection
      backwards_sel <- ols_step_backward_p(full_fit) 
      # Extract final selected model
      backwards_fit <- backwards_sel$model
      # Store the selected variables!
      model_results[[counter]] <- backwards_fit
      
      # Store the rho used for this model fit
      rho_results[counter] <- rho_vec[j]
      
      # Update the counter
      counter <- counter + 1
      
  }
  
  # Summarize results
  # Create lists to store info
  bias <- list()
  coverage <- list()
  tp <- list()
  fp <- list()
  typeI_error <- list()
  typeII_error <- list()
  
  # Specify the "true" vars - vars that are associated with y
  true_vars <- paste0("X", 1:5)
  # Specify the null fars - vars not associated with y
  null_Vars <- paste0("X", 6:n_vars)
  
  # Add names to the vector of betas
  names(B) <- paste0("X", 1:n_vars)
  
  for (i in seq_along(model_results)) {
    # Get the variables that were selected
    selected_vars <- names(coef(model_results[[i]]))[-1] # Remove intercept
    
    # Summarize results
    ### Bias
    # Create vector for betas, with NAs for vars not selected
    bhat <- rep(NA, n_vars)
    names(bhat) <- paste0("X", 1:n_vars) # Add names to vector
    
    # Obtain coefs from model results
    selected_coefs <- coef(model_results[[i]])[-1] # Remove intercept
    # Populate bhat vector with the vars from backwards sel
    bhat[names(selected_coefs)] <- selected_coefs
    
    # Calculate bias
    bias[[i]] <- bhat - B
    
    ### Coverage
    # Get the 95% CI
    ci <- confint.default(model_results[[i]], level = 1 - alpha)
    # Remove the intercept from tis
    ci <- ci[-1, , drop = FALSE]
    
    # Create vector for coverage with NAs for vars not selected
    cover_vec <- rep(NA, n_vars)
    names(cover_vec) <- paste0("X", 1:n_vars) # Name this vector
    
    # Name CI vector with selected vars
    ci_vars <- rownames(ci)
    # Populate cover_vec with T/F if 95% CI covered true param value
    cover_vec[ci_vars] <- ci[, 1] <= B[ci_vars] & B[ci_vars] <= ci[, 2]
    # Fill coverage list with these results
    coverage[[i]] <-  cover_vec
    
    ### True Positives
    # Were X1-X5 selected
    tp_vec <- true_vars %in% selected_vars
    # Assign names of X1-X5 for this vec
    names(tp_vec) <- true_vars
    # Add to list of TPs
    tp[[i]] <- tp_vec
    
    ### False Positives
    fp_vec <- null_vars %in% selected_vars
    # Assign names of X6-X20 for this vec
    names(fp_vec) <- null_vars
    # Add to list for FPs
    fp[[i]] <- fp_vec
    
    ### Type I Error
    # Type I error - FP - Times backwards var picked X6-X20
    typeI_error[[i]] <- fp_vec
    
    ### Type II Error
    # Not selecting the true variables X1-X5
    typeII_vec <- !(true_vars %in% selected_vars)
    # Assign names to this
    names(typeII_vec) <- true_vars
    # Add to list for type II
    typeII_error[[i]] <- typeII_vec

  }
  
  # Convert the lists to data frames
  bias_df <- do.call(rbind, bias)
  coverage_df <- do.call(rbind, coverage)
  tp_df <- do.call(rbind, tp)
  fp_df <- do.call(rbind, fp)
  typeI_df <- do.call(rbind, typeI_error)
  typeII_df <- do.call(rbind, typeII_error)

  # Create data frame of model performance results
  model_performance_df <- data.frame(
    variable = paste0("X", 1:n_vars),
    bias = colMeans(bias_df, na.rm = TRUE),
    coverage = colMeans(coverage_df, na.rm = TRUE)
  )
  
  # Create data frame of variable selection results
  selection_df <- data.frame(
    # Col of variables
    variable = c(true_vars, null_vars),
    true_positive_rate = c(colMeans(tp_df), 
                           rep(NA, length(null_vars))),
    false_positive_rate = c(rep(NA, length(true_vars)), 
                            colMeans(fp_df)),
    typeI_error = c(rep(NA, length(true_vars)), 
                    colMeans(typeI_df)),
    typeII_error = c(colMeans(typeII_df), 
                     rep(NA, length(null_vars)))
  )
  
  # Return list of model performance and var selection results
  return(list(
    model_performance = model_performance_df,
    selection_performance = selection_df
  ))
  
  }
}

# Run simulation
sim_results <- sim1()
# Model performance
sim_model_performance <- sim_results$model_performance
# Variable selection performance
sim_varsel_performance <- sim_results$selection_performance

############################## Checks/Playing around ##########################
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
B <- sim_betas
y <- sim1_data0$y
X <- sim1_data0$X
# Change colnames
colnames(X) <-  paste0("X", 1:20)
# Create data frame with outcome and predictors
model_dat <- data.frame(y = y, X)
# Get all variable names
vnames <- names(model_dat)
# Fit regression model - exclude intercept
model_results <- lm(y ~ .,
                         data = model_dat)
### Backwards selection
backwards_sel <- ols_step_backward_p(model_results)
# Extract model
backwards_fit <- backwards_sel$model
# Get the selected variables and remove intercept
selected_vars <- names(coef(backwards_fit))[-1]
# Get the vars with defined betas
true_vars <- paste0("X", 1:5)
# Get the remaining vars
null_vars <- paste0("X", 6:20)
# TP
tp <- sum(true_vars %in% selected_vars)
# FP
fp <- sum(null_vars %in% selected_vars) # percentage
# CI
ci <- confint.default(model_results, level = 1 - 0.05)
# Coverage
coverage <-  ci[, 1][-1] <= B & B <= ci[, 2][-1]

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












