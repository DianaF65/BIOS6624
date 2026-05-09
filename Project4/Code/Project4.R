################################################################################
###                               Project 4 Start                            ###
################################################################################

# Code for Project4

# Libraries
library(hdrm)
library(MASS) # AIC/BIC variable selection
library(olsrr) # Backwards variable selection
library(glmnet) # LASSO
library(tidyverse)
library(kableExtra)

 
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

#################################### For Checking ##############################

# Simulate data
sim_dat <- gen_data(n = 50,
                    # Number of predictors
                    p = 20, 
                    # Non-zero predictors
                    p1 = 5,
                    # Vector of non-zero betas
                    beta = c(0.5/3, 1/3, 1.5/3, 2/3, 2.5/3,
                             rep(0, 15)),
                    # Distribution
                    family = "gaussian",
                    # Correlation structure
                    corr = "exchangeable",
                    # Correlation coefficient
                    rho = 0.7)

# Extract outcome for model
y <- sim_dat$y
# Extract predictors
X <- sim_dat$X
# Change colnames
colnames(X) <-  paste0("X", 1:20)

# Create data frame with outcome and predictors
model_dat <- data.frame(y = y, X)

################################################################################
###   Function for: 
###   Extracting variable selection method results
###   Calculating model and variable selection performance for one model
################################################################################

# This function will extract results from a given model selection and calculate
# model performance and variable selection performance numbers
## Calculate metrics for ONE model


model_var_performance <- function(# Model from selection method
                                  selection_results,
                                  # Variables that were selected from method
                                  selected_vars,
                                  # Coefs of variables that were selected
                                  selected_coefs,
                                  # Model data for CIs (lasso,elastic net)
                                  model_dat = NULL,
                                  # Number of variables
                                  n_vars = 20, # Default value of 20
                                  # Alpha level
                                  alpha = 0.05
                                  ) {
  
  # Define values for betas
  true_betas = c(0.5/3, 1/3, 1.5/3, 2/3, 2.5/3,
                        rep(0, 15))
  
  # Name the true betas vector for further down
  names(true_betas) <- paste0("X", 1:n_vars)
  
  # Specify the "true" vars - vars that are associated with y
  true_vars <- paste0("X", 1:5)
  # Specify the null vars - vars not associated with y
  null_vars <- paste0("X", 6:n_vars)
  
  # Calculate metrics
  ### Bias
  # Create vector for betas, with NAs for vars not selected
  bhat <- rep(0, n_vars)
  names(bhat) <- paste0("X", 1:n_vars) # Add names to vector
  
  # Populate bhat vector with the vars from selection method - rest as NA
  bhat[selected_vars] <- as.numeric(selected_coefs)
  
  # Calculate bias
  bias <- bhat - true_betas
  
  ### Coverage
  # Create vector for coverage, NAs for vars not selected
  coverage <- rep(0, n_vars)
  names(coverage) <- paste0("X", 1:n_vars) # Name this vector
  
  # Determine if we are working with lasso/elastic net
  if (inherits(selection_results, "cv.glmnet")){
    
    # Check if vars were indeed selected and that modeldat is not null
    if (length(selected_vars) > 0 && !is.null(model_dat)) {
      
      
      # Refit the selected vars with lm()
      refit_model <- lm(as.formula(paste("y ~", paste(selected_vars, 
                                           collapse = " + "))), 
                        data = model_dat)
      
      # Obtain the CI from the refit model
      ci <- confint.default(refit_model, level = 1 - alpha)
      # Remove the intercept
      ci <- ci[-1, , drop = FALSE] 
      
      # Populate coverage vectorwith T/F if 95% CI covered true param value
      coverage[rownames(ci)] <- ci[, 1] <= true_betas[rownames(ci)] & 
        true_betas[rownames(ci)] <= ci[, 2]
      
    } 
    
    ### Get unconditional coverage
    # Variables not selected
    not_selected <- setdiff(paste0("X", 1:n_vars), selected_vars)
    
    # If true variable X1-X5 was not selected, coverage = 0
    coverage[intersect(not_selected, true_vars)] <- 0
    
    # If null variable X6-X20 was not selected, coverage = 1
    coverage[intersect(not_selected, null_vars)] <- 1
    
    
  } else {
      # Backwards, AIC, or BIC
      
      # Get the 95% CI
      ci <- confint.default(selection_results, level = 1 - alpha)
      # Remove the intercept from tis
      ci <- ci[-1, , drop = FALSE]
      
      # Populate cover_vec with T/F if 95% CI covered true param value
      coverage[rownames(ci)] <- ci[, 1] <= true_betas[rownames(ci)] & # Greater than lower bound
        true_betas[rownames(ci)] <= ci[, 2] # Less than upper bound
      
      # Variables not selected
      not_selected <- setdiff(paste0("X", 1:n_vars), selected_vars)
      
      # If a true variable X1-X5 was not selected, coverage = 0
      coverage[intersect(not_selected, true_vars)] <- 0
      
      # If a null variable X6-X20 was not selected, coverage = 1
      coverage[intersect(not_selected, null_vars)] <- 1
  }
  
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
  
  # Create data frame results depending on 
  
  
  # Create data frame of model performance results
  model_performance_df <- data.frame(
    variable = c(true_vars, null_vars),
    # Add the beta hats to assess samp size for final sim run
    beta_hat = bhat,
    bias = bias,
    coverage = coverage)
  
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


################################################################################
###   Function for: 
###   Handling the different variable selection methods
###   If else statements for: 
###     - (1) Backwards, (2) AIC/BIC, (3) LASSO/Elastic Net
### The 
################################################################################


model_var_summary <- function(# Resulting model from variable sel method
                                  selection_results,
                                  # Model data
                                  model_dat = NULL) {

  # Create data frames to return - will vary based on which type of variable sel
  # First check which object we are working with:
  # LASSO
  if (inherits(selection_results, "cv.glmnet")) {
    # Better than (class(selection_results) == "cv.glmnet")
    
    # List to store results for lambda.1se and lamnda.min
    lambdas_results <- list()
    
    # For loop for lambda.1se and lambda.min
    
    for (i in c("lambda.1se", "lambda.min" )) {
      
      # Get estimates for all variables 
      sol <- coef(selection_results, s = i)
      # Remove intercept and convert to vector
      coefficients <- as.numeric(sol[-1, 1])
      names(coefficients) <- rownames(sol)[-1]
      
      # Variables that were selected
      selected_vars <- names(coefficients)[coefficients != 0]
      
      # Coefficients only for selected variables
      selected_coefs <- coefficients[selected_vars]
      
      # Return the summaries
      lambdas_results[[i]] <- model_var_performance(
                                  selection_results = selection_results,
                                  selected_vars = selected_vars,
                                  selected_coefs = selected_coefs,
                                  model_dat = model_dat)
      
    }
    # Return list of results with different lamdbdas
    return(lambdas_results)
    
    # For AIC and BIC
  } else  if (inherits(selection_results, "lm")) { 
    
    # AIC/BIC from step() returns lm model
    varsel_model <- selection_results
    
    selected_vars <- names(coef(varsel_model))[-1]
    selected_coefs <- coef(varsel_model)[-1]
    
    performance_df <- model_var_performance(
      selection_results = varsel_model,
      selected_vars = selected_vars,
      selected_coefs = selected_coefs,
      model_dat = model_dat
    )
    
    return(performance_df)
    
    
  } else {
    
    # For backwards selection
    # Get the model results
    varsel_model <- selection_results$model
    
    # Get the variables that were selected
    selected_vars <- names(coef(varsel_model))[-1] # Remove intercept
    
    # Estimates for vars that were selected
    selected_coefs <- varsel_model$coefficients[-1]
    
    # Obtain model performance and variable sel metrics
    performance_df <- model_var_performance(
      selection_results = varsel_model,
      selected_vars = selected_vars,
      selected_coefs = selected_coefs,
      model_dat = model_dat)
    
    # Return df of performance
    return(performance_df)
    
  }
  
}

# Do a test run
full_fit <- lm(y ~ ., data = model_dat)

# BACKWARDS
backwards_sel <- ols_step_backward_p(full_fit)
bk_test_eval <- model_var_summary(backwards_sel)
bk_model_performance <- bk_test_eval$model_performance
bk_var_sel_performance <- bk_test_eval$selection_performance

# AIC
aic_sel <- step(full_fit,
                trace = 0,
                direction = 'backward',
                k = 2 # AIC
)
aic_test_eval <- model_var_summary(aic_sel)
aic_model_performance <- aic_test_eval$model_performance
aic_var_sel_performance <- aic_test_eval$selection_performance


# BIC
bic_sel <- step(full_fit,
                trace = 0,
                direction = 'backward',
                k = log(nrow(model_dat)) # BIC
)
bic_test_eval <- model_var_summary(bic_sel)
bic_model_performance <- bic_test_eval$model_performance
bic_var_sel_performance <- bic_test_eval$selection_performance

# LASSO
# Create predictor matrix
x <- model.matrix(y ~ ., data = model_dat)[,-1] # Remove the intercept
# Extract outcome
y <- model_dat$y
# Call to cv.glmnet
lasso_sel <- cv.glmnet(x, y, 
                       alpha = 1, # Lasso
                       standardize = TRUE)

lasso_test_eval <- model_var_summary(lasso_sel, model_dat = model_dat)
# Get 1se results
lasso_1se_model_perf <- lasso_test_eval$lambda.1se$model_performance
lasso_1se_var_perf <- lasso_test_eval$lambda.1se$selection_performance
# Get min results
lasso_min_model_perf <- lasso_test_eval$lambda.min$model_performance
lasso_min_var_perf <- lasso_test_eval$lambda.min$selection_performance

# ELASTIC NET
# Call to cv.glmnet
elastic_sel <- cv.glmnet(x, y, 
                       alpha = 0.5, # Elastic net
                       standardize = TRUE)

elastic_test_eval <- model_var_summary(elastic_sel, 
                                       model_dat = model_dat)
# Get 1se results
elastic_1se_model_perf <- elastic_test_eval$lambda.1se$model_performance
elastic_1se_var_perf <- elastic_test_eval$lambda.1se$selection_performance
# Get min results
elastic_min_model_perf <- elastic_test_eval$lambda.min$model_performance
elastic_min_var_perf <- elastic_test_eval$lambda.min$selection_performance

################################################################################
###   Function for: 
###   Generating data
###   Fitting lm() model
###   Using model_var_summary() on variable selection method
################################################################################


# This function will specify the sample size and rhos, fit the lm() model, 
# and will use the model_var_performance() function

# 

# Simulation function that will iterate N times
all_var_sel_methods <- function(# Desired sample size
                                    n = 250,
                                    # Number of variables
                                    n_vars = 20,
                                    # Desired rho
                                    rho = 0.35 ) {
  
  # List to store results from variable selection methods
  selection_list <- list()
  
  # Simulate data
  sim_dat <- gen_data(n = n,
                      # Number of predictors
                      p = n_vars, 
                      # The values for the "true" betas
                      beta = c(0.5/3, 1/3, 1.5/3, 2/3, 2.5/3,
                               rep(0, 15)),
                      # Non-zero predictors
                      p1 = 5,
                      # Distribution
                      family = "gaussian",
                      # Correlation structure
                      corr = "exchangeable",
                      # Correlation coefficient
                      rho = rho)
  
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
  
  ### BACKWARDS SELECTION
  backwards_sel <- ols_step_backward_p(full_fit) 
  # Get model metrics n shit
  # backwards_result <- model_var_summary(backwards_sel)
  selection_list[["backward"]] <- model_var_summary(backwards_sel,
                                                    model_dat = model_dat)
  
  ### AIC
  aic_sel <- step(full_fit,
                  trace = 0,
                  direction = 'backward',
                  k = 2 # AIC
  )
  selection_list[["AIC"]] <- model_var_summary(aic_sel,
                                               model_dat = model_dat)
  
  ### BIC
  bic_sel <- step(full_fit,
                  trace = 0,
                  direction = 'backward',
                  k = log(nrow(model_dat)) # BIC
  )
  selection_list[["BIC"]] <- model_var_summary(bic_sel,
                                               model_dat = model_dat)
  #### SHRINKAGE METHODS
  
  ### LASSO
  # Create predictor matrix
  x <- model.matrix(y ~ ., data = model_dat)[,-1] # Remove the intercept
  
  # Extract outcome
  y <- model_dat$y
  
  # Use the same fold id for lasso and elatic net for better comparison
  foldid <- sample(rep(1:10, length.out = nrow(x)))
  
  # Call to cv.glmnet
  lasso_sel <- cv.glmnet(x, y, 
                         alpha = 1, # Lasso
                         standardize = TRUE,
                         foldid = foldid)
  # Get model metrics
  selection_list[["lasso"]] <- model_var_summary(lasso_sel, 
                                                 model_dat = model_dat)
  
  ### ELASTIC NET
  # Call to cv.glmnet
  elastic_sel <- cv.glmnet(x, y, 
                         alpha = 0.5, # Elastic
                         standardize = TRUE,
                         foldid = foldid)
  # Get model metrics
  selection_list[["elastic_net"]] <- model_var_summary(elastic_sel,
                                                       model_dat = model_dat)

  return(selection_list)
}


# Test run this
sim1_test <- all_var_sel_methods()


################################################################################
###   Function for: 
###   Results for one profile 
################################################################################

# Create all 6 possible profiles
profiles <- expand.grid(
  rho = c(0, 0.35, 0.7),
  N = c(250, 500)
)

# Create function to iterate through the profiles
sim_one_profile <- function(nsim = 10, profile) {
  
  # Make sure we only evaluating one profile
  stopifnot(nrow(profile) == 1)
  
  # List to store results
  one_profile <- list() 
  
  # Iterate through number of sims
  for (iteration in 1:nsim) {
    
    # Temporary holder for given profile
    tmp <- all_var_sel_methods(
      n = profile[, "N"],
      rho = profile[, "rho"]
    )
    
    # Append this iteration to the list
    one_profile[[iteration]] <- tmp
    
  }
  
  # Return results for one profile
  return(one_profile)
}
  

### Test
test_one_profile <- sim_one_profile(nsim = 50, 
                                # rho = 0 and n = 250
                                profile = profiles[1, ])

# This contains results for all 5 variable selection methods for 10 sims

################################################################################
###   Function for: 
###   Summarizing results for simulation run for given profile
##      - Aggregating/combining results for each method
###   Input will be the list of profiles of interest
###     - 
################################################################################

# Within each profile, a list the length of the nsims, within these 
# the results for each variable selection method


# Function for summarizing a given profile
summarize_profiles <- function(profile_results) {
  
  # Create vectors for the 3 types of selecion methods
  methods <- c("backward", "AIC", "BIC")
  lambda_methods <- c("lasso", "elastic_net")
  lambdas <- c("lambda.1se", "lambda.min")
  
  # List to store the summary
  summary_list <- list()
  
  # I was encountering errors due to NAs with some of the calculations
  ## Create a helper function to return means without accounting for NAs
  mean_na <- function(x) {
    if(all(is.na(x))) {
      return(NA) 
    } else {
      return(mean(x, na.rm = TRUE))
    }
  }
  
  
  # Backward, AIC, BIC
  for (method in methods) {
    
    # Row bind all results for model performance
    model_perf <- do.call(
      rbind,
      # Access the model performance results in the list
      lapply(profile_results, function(x) x[[method]]$model_performance)
    )
    
    # Combine results for selection performance
    selection_perf <- do.call(
      rbind,
      # Access the selection performance results in the list
      lapply(profile_results, function(x) x[[method]]$selection_performance)
    )
    
    # Column bind the bias and coverage from model_perf df above
    model_summary <- aggregate(
      # One column for each variable
      cbind(bias, coverage) ~ variable,
      data = model_perf,
      # Find the column means
      FUN = mean_na,
      na.action = na.pass
    )
    
    # Columb nind the FP, TP, error rates from selection_perf df above
    selection_summary <- aggregate(
      # One column for each variable
      cbind(true_positives, false_positives, typeI_error, typeII_error) ~
        variable,
      data = selection_perf,
      # Find the column means
      FUN = mean_na,
      na.action = na.pass
    )
    
    # Combine everything, model performance and selection performance into one obj
    summary_list[[method]] <- list(
      model_performance = model_summary,
      selection_performance = selection_summary
    )
  }
  
  # Lasso and Elastic Net
  for (method in lambda_methods) {
    
    for (lambda in lambdas) {
      # Row bind all results for model performance
      model_perf <- do.call(
        rbind,
        lapply(profile_results, function(x) x[[method]][[lambda]]$model_performance)
      )
      
      # Row bind all results for var sel performance
      selection_perf <- do.call(
        rbind,
        lapply(profile_results, function(x) x[[method]][[lambda]]$selection_performance)
      )
      
      # Column bind the bias and coverage from model_perf df above
      model_summary <- aggregate(
        cbind(bias, coverage) ~ variable,
        data = model_perf,
        FUN = mean_na,
        na.action = na.pass
      )
      
      # Column nind the FP, TP, error rates from selection_perf df above
      selection_summary <- aggregate(
        cbind(true_positives, false_positives, typeI_error, typeII_error) ~ variable,
        data = selection_perf,
        FUN = mean_na,
        na.action = na.pass
      )
      
      # Combine everything, model performance and selection performance into one obj
      summary_list[[paste(method, lambda, sep = "_")]] <- list(
        model_performance = model_summary,
        selection_performance = selection_summary
      )
    }
  }
  
  return(summary_list)
}



################################################################################
###   Function for: 
###   Summarizing results for ALL 6 profiles
################################################################################

# Simulation results
sim_results <- function(# Simulations
                        nsim = 10,
                        # Profile
                        profiles) {
  
  # List to store results from all profiles
  all_profiles_results <- list()
  
  # Iterate through the profiles
  for (p in 1:nrow(profiles)) {
    
    # Apply sim_one_profile to the pth profile
    all_profiles_results[[p]] <- sim_one_profile(
      nsim = nsim, # Take in nsim argument
      profile = profiles[p, ] # pth profile
    )
  }
  
  # Name the results
  names(all_profiles_results) <- paste0(
    "n", profiles$N, "_rho", profiles$rho
  )
  
  return(all_profiles_results)
}

### Test
# Look at one profile
all_profile_summaries$n250_rho0.35$backward$model_performance


################################################################################
###                               Testing                                    ###
################################################################################

### 1. One simulated dataset
test_one <- all_var_sel_methods(n = 250, rho = 0.35)

### 2. One profile, nsim iterations
test_profile <- sim_one_profile(nsim = 10, profile = profiles[1, ])


############### 
### 3. Two profiles, two iterations each
### Doing a pilot run go get sampe size for sim run

test_profiles <- profiles[1:6, ]
# rho   N
# 1 0.00 250
# 2 0.35 250

# List
test <- list()

# Test profiles
pilot_results <- sim_results(
  nsim = 200,
  profiles = test_profiles
)
  
  

# 4. Summarize one profile
test_summary1 <- summarize_profiles(test_all_profiles[[1]])
# Summarize the n = 250 and rho = 0

# Create data frame to assess multiple methods
all_model_performance <- list()

# For the 2 profiles
for (p in seq_along(test_all_profiles)) {
  # Get the pth profile summary
  profile_summary <- summarize_profiles(test_all_profiles[[p]])
  # The pth profile info
  profile_info <- test_profiles[p, ]
  
  # List of all model performance for all selection methods
  model_perf_list <- list(
    Backward = profile_summary$backward$model_performance,
    AIC = profile_summary$AIC$model_performance,
    BIC = profile_summary$BIC$model_performance,
    LASSO_lambda.1se = profile_summary$lasso_lambda.1se$model_performance,
    LASSO_lambda.min = profile_summary$lasso_lambda.min$model_performance,
    Elastic_lambda.1se = profile_summary$elastic_net_lambda.1se$model_performance,
    Elastic_lambda.min = profile_summary$elastic_net_lambda.min$model_performance
  )
  
  # Create a data frame
  all_model_performance[[p]] <- do.call(
    # Combine into rows
    rbind,
    # Create data frame with the cols of interst
    lapply(names(model_perf_list), function(m) {
      data.frame(
        n = profile_info$N,
        rho = profile_info$rho,
        method = m,
        model_perf_list[[m]]
      )
    })
  )
}

all_model_performance <- do.call(rbind, all_model_performance)
row.names(all_model_performance) <- NULL

# Adjust the order of the variables
# Make variable a factor
all_model_performance$variable <- factor(all_model_performance$variable,
                                     levels = paste0("X", 1:20))
# Sort
all_model_performance <- all_model_performance[order(all_model_performance$variable,
                                            all_model_performance$method), ]

### Display something Dr. S has..
# Coverage for X1-X5 for backwards, AIC, BIC
coverage <- all_model_performance %>% 
  filter(variable %in% paste0("X", 1:5)) %>% 
  # filter(method %in% c("Backward", "AIC", "BIC")) %>% 
  group_by(n, rho, method, variable) %>% 
  summarise(prob_cover = mean(coverage) * 100)

# Df for now to look at
cov_df <- pivot_wider(coverage,
                      names_from = variable,
                      values_from = prob_cover)

# Display with kable
knitr::kable(
  coverage,
    pivot_wider(
      names_from = variable,
      values_from = prob_cover,
      format = "html"
    ),
  digits = 3
)

# Bias for X1-X5 for backwards, AIC, BIC
bias <- all_model_performance %>% 
  # filter(variable %in% paste0("X", 1:5)) %>% 
  # filter(method %in% c("Backward", "AIC", "BIC")) %>% 
  group_by(n, rho, method, variable) %>% 
  summarise(prob_bias = mean(bias))

# Df for now to look at
bias_df <- pivot_wider(bias,
                      names_from = variable,
                      values_from = prob_bias)


# Selection performance
# Create data frame to assess multiple methods
all_selection_performance <- list()

# For the 2 profiles
for (p in seq_along(test_all_profiles)) {
  # Get the pth profile summary
  profile_summary <- summarize_profiles(test_all_profiles[[p]])
  # The pth profile info
  profile_info <- test_profiles[p, ]
  
  # List of all selection performance for all selection methods
  selection_perf_list <- list(
    Backward = profile_summary$backward$selection_performance,
    AIC = profile_summary$AIC$selection_performance,
    BIC = profile_summary$BIC$selection_performance,
    LASSO_lambda.1se = profile_summary$lasso_lambda.1se$selection_performance,
    LASSO_lambda.min = profile_summary$lasso_lambda.min$selection_performance,
    Elastic_lambda.1se = profile_summary$elastic_net_lambda.1se$selection_performance,
    Elastic_lambda.min = profile_summary$elastic_net_lambda.min$selection_performance
  )
  
  # Create a data frame
  all_selection_performance[[p]] <- do.call(
    # Combine into rows
    rbind,
    # Create data frame with the cols of interst
    lapply(names(selection_perf_list), function(m) {
      data.frame(
        n = profile_info$N,
        rho = profile_info$rho,
        method = m,
        selection_perf_list[[m]]
      )
    })
  )
}

all_selection_performance <- do.call(rbind, all_selection_performance)
row.names(all_selection_performance) <- NULL

# Adjust the order of the variables
# Make variable a factor
all_selection_performance$variable <- factor(all_selection_performance$variable,
                                         levels = paste0("X", 1:20))
# Sort
all_selection_performance <- all_selection_performance[order(all_selection_performance$variable,
                                                     all_selection_performance$method), ]

# False postiives ==TYPE I ERROR 
fps <- all_selection_performance %>% 
  filter(!(variable %in% c("X1", "X2", "X3", "X4", "X5"))) %>% 
  group_by(n, rho, method, variable) %>% 
  summarise(prob_fps = mean(false_positives) * 100)

# Df for now to look at
fps_df <- pivot_wider(fps,
                      names_from = variable,
                      values_from = prob_fps)

# Type II error
errorII <- all_selection_performance %>% 
  filter(variable %in% c("X1", "X2", "X3", "X4", "X5")) %>% 
  group_by(n, rho, method, variable) %>% 
  summarise(prob_errorII = mean(typeII_error) * 100)

errorII_df <- pivot_wider(
  errorII, 
  names_from = variable,
  values_from = prob_errorII
)

# Type I error
errorI <- all_selection_performance %>% 
  filter(!(variable %in% c("X1", "X2", "X3", "X4", "X5"))) %>% 
  group_by(n, rho, method, variable) %>% 
  summarise(prob_errorI = mean(typeI_error) * 100)

errorI_df <- pivot_wider(
  errorI, 
  names_from = variable,
  values_from = prob_errorI
)


################################################################################
###                               SIMULATION RUN                             ###
################################################################################

# Load the future library for parallelization
library(future)
library(future.apply)

# Define the number of profiles 
n_profiles <- nrow(profiles)

plan(multisession, workers = 5)
# I have 8 cores on my computer, this leaves 2 for doing other work

# List to store results for all 6 profiles
sim_results <- list()

# Seed for reproducibility
set.seed(6456)

## SIMULATION RUN
run_time <- system.time({
  
  sim_results <- future_lapply(
    1:n_profiles,
    function(i) {
      sim_one_profile(
        nsim = 2500, 
        profile = profiles[i, ]
      )
    },
    future.seed = TRUE
  )
  
})

# Name the results by profile
names(sim_results) <- paste0("n", profiles$N, "_rho", profiles$rho)

# Save results immediately
saveRDS(sim_results, "simulation_results_2500.rds")

# Return to sequential plan when done
plan(sequential)

# Print run time
run_time

# Then check

length(sim_results)
length(sim_results[[1]])
names(sim_results)


## Note on nsim - I am mainly interested in evaluating bias, therefore
# I will run 2500 repititions for the simulation
# 


###### DATA FRAME CREATION FOR SIMULATION RUN RESULTS
# Create data frame to assess multiple methods
sim_model_performance <- list()

# For the 2 profiles
for (p in seq_along(sim_results) {
  # Get the pth profile summary
  profile_summary <- summarize_profiles(sim_results[[p]])
  # The pth profile info
  profile_info <- profiles[p, ]
  
  # List of all model performance for all selection methods
  model_perf_list <- list(
    Backward = profile_summary$backward$model_performance,
    AIC = profile_summary$AIC$model_performance,
    BIC = profile_summary$BIC$model_performance,
    LASSO_lambda.1se = profile_summary$lasso_lambda.1se$model_performance,
    LASSO_lambda.min = profile_summary$lasso_lambda.min$model_performance,
    Elastic_lambda.1se = profile_summary$elastic_net_lambda.1se$model_performance,
    Elastic_lambda.min = profile_summary$elastic_net_lambda.min$model_performance
  )
  
  # Create a data frame
  model_performance[[p]] <- do.call(
    # Combine into rows
    rbind,
    # Create data frame with the cols of interst
    lapply(names(model_perf_list), function(m) {
      data.frame(
        n = profile_info$N,
        rho = profile_info$rho,
        method = m,
        model_perf_list[[m]]
      )
    })
  )
}

all_model_performance <- do.call(rbind, all_model_performance)
row.names(all_model_performance) <- NULL

# Adjust the order of the variables
# Make variable a factor
siml_model_performance$variable <- factor(sim_model_performance$variable,
                                         levels = paste0("X", 1:20))
# Sort
sim_model_performance <- sim_model_performance[order(sim_model_performance$variable,
                                                     sim_model_performance$method), ]

### Display something Dr. S has..
# Coverage for X1-X5 for backwards, AIC, BIC
coverage <- sim_model_performance %>% 
  filter(variable %in% paste0("X", 1:5)) %>% 
  # filter(method %in% c("Backward", "AIC", "BIC")) %>% 
  group_by(n, rho, method, variable) %>% 
  summarise(prob_cover = mean(coverage) * 100)

# Df for now to look at
cov_df <- pivot_wider(coverage,
                      names_from = variable,
                      values_from = prob_cover)


# Bias for X1-X5 for backwards, AIC, BIC
bias <- sim_model_performance %>% 
  # filter(variable %in% paste0("X", 1:5)) %>% 
  # filter(method %in% c("Backward", "AIC", "BIC")) %>% 
  group_by(n, rho, method, variable) %>% 
  summarise(prob_bias = mean(bias))

# Df for now to look at
bias_df <- pivot_wider(bias,
                       names_from = variable,
                       values_from = prob_bias)



############################## Checks/Playing around ##########################
# Creating simulated data
set.seed(646)
sim1_data0 <- gen_data(n = 250,
                       # Number of predictors
                       p = 20, 
                       # Non-zero predictors
                       p1 = 5,
                       # Vector of non-zero betas
                       beta = c(0.5/3, 1/3, 1.5/3, 
                                2/3, 2.5/3,
                                rep(0, 15)),
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


