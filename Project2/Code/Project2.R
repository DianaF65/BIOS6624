# Code for Project0

# Libraries
library(powertools)
library(here)
library(readr)

# Read in the data
dat <- read_csv(here("Project2", "Data", "PrelimData.csv"))

### Calculating correlations for power analyses
# Correlations are used in power calculations when measuring the strength
# and direction of a relationship between two continuous variables instead of
# comparing group means

colnames(dat)
# CVLT_CNG3 - Measure of episodic memory
#   - Recall/Recognition cognitive outcome variable
# CORT_CNG3 - Measure of corticol thickness
# IL_6 - A novel Cytokine/chemokine of interest
# MCP_1 - A tested cytokine/chemokine of interest

#### Aim 1
# Correlation between IL6 and episodic memory
cor(dat$IL_6, dat$CVLT_CNG3)
# -0.2585617

# Correlation between MCP1 and episodic memory
cor(dat$MCP_1, dat$CVLT_CNG3)
# -0.3183585

# Correlation between IL6 and corticol thickness
cor(dat$IL_6, dat$CORT_CNG3)
# -0.5993366

# Correlation between MCP1 and corticol thickness
cor(dat$MCP_1, dat$CORT_CNG3)
# -0.685273

#### Correlation between cytokines/chemokines
# Correlation between IL6 and MCP1
cor(dat$IL_6, dat$MCP_1)
# 0.9335558
# - Pretty highly correlated so perhaps we can reduce the number of chemokines/
#cytokines to just 1 instead of examining all 4 which will also address the
# multiple testing issue

### OR JUST DO THIS
cor(dat)

#################################### Aim 1a ####################################

### Pre power analysis calculations
# For Aim 1, correlations range from -0.26 to -0.69
# We will use the corr.1samp() function
# NOTE: The correlation coefficient becomes an effect size when squared
# AKA the r-squared value
test1 <- corr.1samp(N = 175,
                    rho0 = 0,
                    rhoA = (0.26)^2,
                    alpha = 0.05,
                    power = NULL,
                    sides = 1)

test1

test2 <- corr.1samp(N = 175,
                    rho0 = 0,
                    rhoA = (0.5)^2,
                    alpha = 0.05,
                    power = NULL,
                    sides = 1)

test2

test3 <- corr.1samp(N = NULL,
                    rho0 = 0,
                    rhoA = (0.5)^2,
                    alpha = 0.05,
                    power = 0.8,
                    sides = 1)

test3


# Vector of correlation values
r <- seq(0.2, 0.80, by = 0.05)
# Square these values for r^squared vals
r_squared <- r^2

# Apply the corr.1samp() function to this range of effect sizes
# What are the possible values of power for the range of correlations with a set
# sample size of 175?
aim1_power <- vapply(r_squared, function(k) {
  corr.1samp(
    N = 175,
    rho0 = 0,
    rhoA = k,
    alpha = 0.05,
    power = NULL,
    sides = 1
  )}, numeric(1)
)

# What are the required sample sizes to achieve 80% power for the range of 
# correlations?
aim1_sampsizes <- vapply(r_squared, function(k) {
  corr.1samp(
    N = NULL,
    rho0 = 0,
    rhoA = k,
    alpha = 0.05,
    power = 0.8,
    sides = 1
  )}, numeric(1)
)
# Maybe this is relevant for Aim 2 more


# Create a data frame of the powers
aim1_df <- data.frame(`Correlation` = r,
                      `R-Squared` = r_squared,
                      `Power` = round(aim1_power, 3),
                      check.names = FALSE)


# This power analysis does not account for multiple predictors, baseline adjustment,
# or interactions
# Make sure we state that this function is not adjusting for baseline or other
# covariates


##################################### Aim 2 ####################################

### From doc: We will use global SUVR in our analyses, and we will consider
# a dichotomous (amyloid positive or negative) outcome as an alt. analytic strategy
# Also, Amyloid is measured at baseline. So our interaction will be
# amyloid(Y/N) x IL-6 (or MCP-1)
# This will also be a lot easier for interpretation

#1) Although we have proposed using a continuous SUVR variable for amyloid PET, we will
# examine the utility of incorporating a dichotomous variable for amyloid positivity using recent cutoffs described
# in the literature (73). 

# TheSUVR threshold of 1.10 had a sensitivity of 97% and a specificity of 100% for 
# identifying individuals with clinically significant amyloid pathology 
# (defined as moderate to frequent neuritic plaques by CERAD criteria) (73) 

# (73) reported a rho = 0.79 for correlatino between PET and postmort amyloid
# burden within one year 

# We will use a dichotomous variable of amyloid deposition
# Create an "interaction" variable for Aim 2
# Interaction between IL-6 and amyloid (Y/N)

# We have these possible correlations or effect sizes for the outcomes and 
# PEVs
# Interaction terms are known to have significantly smaller effect sizes
# than main effects
# So what are the sample sizes required to achieve 80% power for effect 
# sizes that are smaller than we have for the main effects?

# Define the new vector of correlations
aim2_r <- seq(0.05, 0.6, by = 0.05)

# Get the rsqured values
aim2_rsquared <- aim2_r^2

aim2_sampsizes <- vapply(aim2_rsquared, function(k) {
  corr.1samp(
    N = NULL,
    rho0 = 0,
    rhoA = k,
    alpha = 0.05,
    power = 0.8,
    sides = 1
  )}, numeric(1)
)

# Combine this into a data frame
aim2_df <- data.frame(`Correlation` = aim2_r,
                      `R-Squared` = aim2_rsquared,
                      `Sample Sizes` = aim2_sampsizes,
                      check.names = FALSE)
# The smallest effect size seen with the main effects is 0.04 and the sample
# size needed to obtain 80% power is 3,861 subjects
# The aim is to recruit 175 subjects and with this you can achieve power of
aim2_power <- vapply(aim2_rsquared, function(k) {
  corr.1samp(
    N = 175,
    rho0 = 0,
    rhoA = k,
    alpha = 0.05,
    power = NULL,
    sides = 1
  )}, numeric(1)
)
# 80% with an effect size of 



