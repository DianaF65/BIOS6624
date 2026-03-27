# Code for Project0

# Libraries
library(powertools)
library(here)
library(readr)
library(dplyr)

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
aim1_r <- seq(0, 0.80, by = 0.02)
# Square these values for r^squared vals
# r_squared <- r^2

# Adjust this based on what we observe with power
# Reaches 1 after 0.36
aim1_r2 <- seq(0, 0.36, by = 0.02)

# Create some more visually appealing r_squared values to report
# It's in the range of 0.04 to 0.64
# Perhaps 0.1 to 0.6 by 0.5
# aim1_rsqrd <- seq(0.02, 0.2, by = 0.02)

# Apply the corr.1samp() function to this range of effect sizes
# What are the possible values of power for the range of correlations with a set
# sample size of 175?
aim1_power <- vapply(aim1_r2, function(k) {
  corr.1samp(
    N = 175,
    rho0 = 0,
    rhoA = k,
    alpha = 0.05,
    power = NULL,
    sides = 1
  )}, numeric(1)
)

# Create a data frame of the powers
aim1_power_df <- data.frame(`Effect Size` = aim1_r2,
                      `Power` = round(aim1_power, 3),
                      check.names = FALSE)


# This power analysis does not account for multiple predictors, baseline adjustment,
# or interactions
# Make sure we state that this function is not adjusting for baseline or other
# covariates

# Transpose data frame for report
# > a <- data.frame(t(aim1_power_df))


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
# (73) reported a rho = 0.79 for correlation between PET and postmort amyloid
# burden within one year 

# We will use a dichotomous variable of amyloid deposition
# Interaction between IL-6 and amyloid (Y/N)

#### From Tasha
# ...use the dichotomous amyloid variable as a simplifying assumption for power. If 
# you use that dichotomous variable, then what 2 things will you essentially be 
# estimating in terms of the relationship between the inflammatory marker and the 
# outcome? You can test the difference between those things using a very similar 
# package to the one you are using for your other calculations.

# If we dichotomize, the interaction is interpreted as the change in the outcomes
# for the signficant amyloid depo. group compared to the group with no sig. 
# amyloid dep
# We will be testing the difference in the outcome for those with sig. AD compared
# to those with non sig. AD
# So we will use corr.2smp to get power for comparing the correlation
# coefficient of both of these groups
# We have to get the correlation coefficient for those with no sig. AD
# And those with sig. AD
# We have the correlation effects for the main effect of the cytokines/chemokines
# So we will just have to do it based on that I believe..
# We are going to suppose that the correlation coefficients for the AD
# group is larger than that of the non AD group

# How to choose n? There are 125 aMCI and 50 HC planned
# But we also have: Literature indicates that AD pathology may be present in
# functionally normal adults decades before a diagnosis of MCI or AD (42, 43)

# How much do we expect the correlations between the two groups to differ?
# Pick values from: 0.1 - 0.5

# Define the new vector of correlations
# This is for one group say AD - which we expect to be lower than AD +
aim2_r <- seq(0.1, 0.3, by = 0.1)
# Vector of correlations for AD + group
# Testing for changes of 0.1 - 0.5
aim2_r2 <- seq(0.2, 0.6, by = 0.1)

# Get the rsqaured values
# aim2_rsquared <- aim2_r^2
  
#### WAY 2
# Obtain the different combinations of rho1 and rho2
aim2_powerdf <- expand.grid(rho1 = aim2_r, rho2 = aim2_r2)

# Calculate power
aim2_powerdf$power <- mapply(function(k, i) {
  corr.2samp(
    n1 = 175,
    n.ratio = 1,
    rho1 = k,
    rho2 = i,
    alpha = 0.05,
    power = NULL,
    sides = 2
  )
}, aim2_powerdf$rho1, aim2_powerdf$rho2)

# Add a column with differences
aim2_powerdf$change <- aim2_powerdf$rho1 - aim2_powerdf$rho2

aim2_powerdf$change <- abs(round(aim2_powerdf$change, 2))

# Round power
aim2_powerdf$power <- round(aim2_powerdf$power, 4)

# Arrange by change
aim2_powerdf <- aim2_powerdf %>% 
  arrange(change)



