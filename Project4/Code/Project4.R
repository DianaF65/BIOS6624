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
###                             Simulation Functions                         ###
################################################################################

### Be a copy cat and copy Dr. Sevick's strategy 
# Three functions for three purposes: 

# 1.) Innermost nested function that will extract model results from variable 
# selection method and calculate model performance and variable selection numbers

# 2.) Outer function that will generate data, fit the lm() model, and use the
# innermost function to obtain summaries

# 3.) Outermost function that will run through the 6 different possible profiles


################################################################################
###   Function for: 
###   Extracting variable selection method results
###   Summarizing model and variable selection performance
################################################################################

# This function will extract results from a given model selection and calculate
# model performance and variable selection performance numbers
## Calculate metrics for ONE model

model_var_performance <- function(# Resulting model from variable sel method
                                  selection_results,
                                  # Define values for betas
                                  true_betas = c(0.5/3, 1/3, 1.5/3, 
                                                 2/3, 2.5/3,
                                            rep(0, 15)),
                                  # Define number of variables
                                  n_vars = 20,
                                  # Define alpha for calculations
                                  alpha = 0.05) {
  
  # Name the true betas vector for further down
  names(true_betas) <- paste0("X", 1:n_vars)
  
  # Specify the "true" vars - vars that are associated with y
  true_vars <- paste0("X", 1:5)
  # Specify the null fars - vars not associated with y
  null_vars <- paste0("X", 6:n_vars)
  
  # Get the model results
  varsel_model <- selection_results$model
  
  # Get the variables that were selected
  selected_vars <- names(coef(varsel_model))[-1] # Remove intercept
  
  # Summarize results
  # Create lists to store results
  
  
  ### Bias
  # Create vector for betas, with NAs for vars not selected
  bhat <- as.vector(rep(0, n_vars))
  names(bhat) <- paste0("X", 1:n_vars) # Add names to vector
  
  # Populate bhat vector with the vars from selection method - rest as NA
  bhat[selected_vars] <- as.numeric(varsel_model$coefficients)
  
  # Calculate bias
  bias <- bhat - true_betas
  
  ### Coverage
  # Get the 95% CI
  ci <- confint.default(varsel_model, level = 1 - alpha)
  # Remove the intercept from tis
  ci <- ci[-1, , drop = FALSE]
  
  # Create vector for coverage, NAs for vars not selected
  coverage <- rep(NA, n_vars)
  names(coverage) <- paste0("X", 1:n_vars) # Name this vector
  
  # Name CI vector with selected vars
  ci_vars <- rownames(ci)
  # Populate cover_vec with T/F if 95% CI covered true param value
  coverage[ci_vars] <- ci[, 1] <= true_betas[ci_vars] & # Greater than lower bound
     true_betas[ci_vars] <= ci[, 2] # Less than upper bound
  # Fill coverage list with these results
  # coverage <-  ifelse(is.na(cover_vec), 0, 1)
  
  ### True Positives
  # Vector of length 20
  tp_vec <- rep(NA, n_vars)
  names(tp_vec) <- paste0("X", 1:n_vars)
  # Were X1-X5 selected
  tp <- true_vars %in% selected_vars
  tp_vec[true_vars] <- tp
  
  ### False Positives
  fp_vec <- rep(NA, n_vars)
  names(fp_vec) <- paste0("X", 1:n_vars)
  # Were X6-X20 selected
  fp <- null_vars %in% selected_vars
  fp_vec[null_vars] <- fp

  
  ### Type I Error
  # Type I error - FP - Times backwards var picked X6-X20
  typeI_error <- fp_vec
  
  ### Type II Error
  typeII_vec <- rep(NA, n_vars)
  names(typeII_vec) <- paste0("X", 1:n_vars)
  # Not selecting the true variables X1-X5
  typeII <- !(true_vars %in% selected_vars)
  typeII_vec[true_vars] <- typeII
  
  # Create data frame of model performance results
  model_performance_df <- data.frame(
    variable = c(true_vars, null_vars),
    bias = bias,
    coverage = coverage
  )
  
  # Create data frame of variable selection results
  selection_df <- data.frame(
    variable = c(true_vars, null_vars),
    true_positives = tp_vec,
    false_positives = fp_vec,
    typeI_error = typeI_error,
    typeII_error = typeII_vec
  )
  
  # Return list of model performance and var selection results
  return(list(
    model_performance = model_performance_df,
    selection_performance = selection_df
  ))
  
}

# Do a test run
full_fit <- lm(y ~ ., data = model_dat)

# Backwards selection
backwards_sel <- ols_step_backward_p(full_fit)
test_eval <- model_var_performance(backwards_sel)
model_performance <- test_eval$model_performance
var_sel_performance <- test_eval$selection_performance

#### CONTINUE HERE
# LASSO
lasso_sel <- 


################################################################################
###   Function for: 
###   Generating data
###   Fitting lm() model
###   Using model_var_performance() on variable selection method
################################################################################


# This function will specify the sample size and rhos, fit the lm() model, 
# and will use the model_var_performance() function


# Simulation function that will iterate N times
simulation_sampsize_rho <- function(# Desired sample size
                                    n = 250,
                                    # Desired rho
                                    rho = 0.35 ) {
  # Simulate data
  sim_dat <- gen_data(n = n,
                      # Number of predictors
                      p = n_vars, 
                      # Non-zero predictors
                      p1 = 5,
                      # Vector of non-zero betas
                      beta = betas,
                      # Distribution
                      family = "gaussian",
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
  
  ### LASSO
  
}


# Run simulation
sim_results <- sim1(N = 100)
# Model performance
sim_model_performance <- sim_results$model_performance
# Variable selection performance
sim_varsel_performance <- sim_results$selection_performance



################################################################################
###   Function for: 
###   Iterating through the 6 different profiles
###     - 3 rhos and 2 sample sizes
################################################################################

# 


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


# Tests with calculating coverage, etc.
# Example: Vector for coverage with NA and TRUEs
a <- data.frame(A = c(rep(NA, 5), rep(TRUE, 5)),
                B = c(rep(TRUE, 5), rep(NA, 5)),
                C = c(rep(NA, 3), rep(TRUE, 3), rep(NA, 4))
                )
# Taking the colMeans: 
colMeans(a, na.rm = TRUE)
# A B C 
# 1 1 1 
# Or do we want 0s and 1s
b <- data.frame(A = c(rep(0, 5), rep(1, 5)),
                B = c(rep(1, 5), rep(0, 5)),
                C = c(rep(0, 3), rep(1, 3), rep(0, 4))
)
# Taking colMeans: 
# A   B   C 
# 0.5 0.5 0.3 

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












