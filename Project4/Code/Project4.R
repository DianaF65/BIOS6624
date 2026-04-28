################################################################################
###                               Project 4 Start                            ###
################################################################################

# Code for Project4

# Libraries
library(hdrm)

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










